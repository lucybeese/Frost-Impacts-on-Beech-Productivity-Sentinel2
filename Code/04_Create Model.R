##%%%%%%%%%%%%%%%%%%##

### Frost Analysis ###

##%%%%%%%%%%%%%%%%%%##

#%%%%%%%%%%%%%%%%%%%%#

#### Create Model ####

#%%%%%%%%%%%%%%%%%%%%#

#------------------------#
##### Load Libraries #####
#------------------------#

library(data.table)
library(dplyr)
library(lubridate)
library(lme4)
library(car)

#-------------------#
##### Load Data #####
#-------------------#

# Read data from csv files
NDVI_complete <- fread("Data/classified_NDVI_all.csv")
data<- fread('Data/Results_filt.csv')

#%%%%%%%%%%%%%%%%%#

#### Aggregate ####

#%%%%%%%%%%%%%%%%%#

#----------------------------#
#####  Aggregate by Week #####
#----------------------------#

# Convert Date to IDate (lighter than character)
NDVI_complete[, Date := as.IDate(Date)]

# compute year/week columns efficiently
NDVI_complete[, Year := year(Date)]
NDVI_complete[, WeekNum := as.integer(format(Date, "%V"))]

# Count distinct Frost_Class values per ID/week
class_check <- NDVI_complete[, .(n_classes = uniqueN(Frost_Class)), 
                             by = .(ID, Year, WeekNum)]

# See how many weeks have more than one class
table(class_check$n_classes)

#-------------------------------------------#
#####  Aggregate by ID + Year + WeekNum #####
#-------------------------------------------#

weekly_NDVI <- NDVI_complete[, .(
  Mean_NDVI = mean(NDVI, na.rm = TRUE),
  Min_Date = min(Date),        # keep representative Date (start of week)
  Max_Date = max(Date)         # optional, end of week
), by = .(ID, Year, WeekNum, Frost_Class)]

head(weekly_NDVI)

# Clean up
rm(NDVI_complete)
gc()

# Save 
fwrite(weekly_NDVI, "Data/weekly_classified_NDVI.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### Add topographic info ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Things like elevation need to be accounted for
head(data)
head(weekly_NDVI)

# Join max_ndvi and category to data
data_with_ndvi <- weekly_NDVI %>%
  left_join(data, by = c("ID", "Year"))

# Check the result
head(data_with_ndvi)

# Make sure Elevation_Band is an ordered factor
data_with_ndvi$Elevation_Band <- factor(
  data_with_ndvi$Elevation_Band,
  levels = c("<800", "800-1000", "1000-1200", "1200-1400", ">1400"),
  ordered = FALSE
)

head(data_with_ndvi)

# Save
fwrite(data_with_ndvi,"Data/H2_data_allyear.csv")

#%%%%%%%%%%%%%%%%%%%#

#### H2/H3 Model ####

#%%%%%%%%%%%%%%%%%%%#

# Here I want to see if 2019 is lower mean NDVI than all other years. I also want to see if frost classified points have significantly lower 
# mean ndvi than frost classified points in 2019, but no significant difference between frost classes in all other years.
# all of this while controlling for elevation, the specific week of the year, and the point ID

# Factor
data_with_ndvi$Year <- factor(data_with_ndvi$Year)

# Model
mod_int <- lmer(
  Mean_NDVI ~ Year * Frost_Class + Elevation_Band + (1 | ID) + (1 | WeekNum),
  data = data_with_ndvi
)

# Save
saveRDS(mod_int, file = "data/LMM.rds")
summary(mod_int)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: Mean_NDVI ~ Year * Frost_Class + Elevation_Band + (1 | ID) +      (1 | WeekNum)
#    Data: data_with_ndvi
# 
# REML criterion at convergence: -21868672
# 
# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -12.8202  -0.4840   0.0495   0.5371   9.4259 
# 
# Random effects:
#  Groups   Name        Variance  Std.Dev.
#  ID       (Intercept) 0.0009834 0.03136 
#  WeekNum  (Intercept) 0.0133683 0.11562 
#  Residual             0.0024343 0.04934 
# Number of obs: 6936333, groups:  ID, 46414; WeekNum, 32
# 
# Fixed effects:
#                                 Estimate Std. Error         df  t value Pr(>|t|)    
# (Intercept)                    8.116e-01  2.045e-02  2.785e+01   39.688   <2e-16 ***
# Year2019                      -2.763e-02  1.146e-04  6.900e+06 -241.029   <2e-16 ***
# Year2020                       1.348e-02  1.148e-04  6.900e+06  117.464   <2e-16 ***
# Year2021                       1.588e-03  1.179e-04  6.900e+06   13.471   <2e-16 ***
# Year2022                       7.549e-03  1.140e-04  6.900e+06   66.198   <2e-16 ***
# Year2023                       1.526e-02  1.166e-04  6.901e+06  130.864   <2e-16 ***
# Frost_ClassNon-frost          -4.593e-03  3.400e-04  5.547e+04  -13.510   <2e-16 ***
# Elevation_Band800-1000        -1.573e-02  6.718e-04  4.618e+04  -23.408   <2e-16 ***
# Elevation_Band1000-1200       -2.728e-02  6.323e-04  4.617e+04  -43.151   <2e-16 ***
# Elevation_Band1200-1400       -3.594e-02  6.417e-04  4.615e+04  -56.003   <2e-16 ***
# Elevation_Band>1400           -4.073e-02  6.822e-04  4.619e+04  -59.704   <2e-16 ***
# Year2019:Frost_ClassNon-frost  2.705e-02  1.436e-04  6.901e+06  188.345   <2e-16 ***
# Year2020:Frost_ClassNon-frost -2.847e-03  1.436e-04  6.901e+06  -19.819   <2e-16 ***
# Year2021:Frost_ClassNon-frost -2.125e-03  1.471e-04  6.902e+06  -14.448   <2e-16 ***
# Year2022:Frost_ClassNon-frost -1.855e-03  1.428e-04  6.901e+06  -12.992   <2e-16 ***
# Year2023:Frost_ClassNon-frost  1.073e-04  1.455e-04  6.902e+06    0.737    0.461    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#---------------------------#
##### Check assumptions #####
#---------------------------#

set.seed(123)

res <- residuals(mod_int)
fit <- fitted(mod_int)
idx <- sample(seq_along(res), 20000)

# Residuals vs fitted
plot(fit[idx], res[idx],
     pch = 16, cex = 0.3,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2, col = "red")

# Residual QQ plot
qqnorm(res[idx], pch = 16, cex = 0.3)
qqline(res[idx], col = "red", lwd = 2)

# Random effects QQ plots
re_id <- ranef(mod_int)$ID[, 1]
re_week <- ranef(mod_int)$WeekNum[, 1]

re_id_s <- sample(re_id, min(20000, length(re_id)))

qqnorm(re_id_s, pch = 16, cex = 0.4,
       main = "Q-Q plot of random intercepts: ID")
qqline(re_id_s, col = "red")

qqnorm(re_week, pch = 16, cex = 1,
       main = "Q-Q plot of random intercepts: WeekNum")
qqline(re_week, col = "red")

# Multicollinearity among fixed effects
vif(lm(Mean_NDVI ~ factor(Year) * Frost_Class + Elevation_Band,
       data = data_with_ndvi))

# factor(Year)             182.388933  5        1.683061
# Frost_Class                7.336028  1        2.708510
# Elevation_Band             1.116410  4        1.013860
# factor(Year):Frost_Class 599.194083  5        1.895644