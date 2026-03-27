##%%%%%%%%%%%%%%%%%%##

### Frost Analysis ###

##%%%%%%%%%%%%%%%%%%##

#%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### HYPOTHESIS TESTING ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%#

#------------------------#
##### Load Libraries #####
#------------------------#

library(dplyr)
library(tidyr)
library(data.table)
library(car)
library(lmtest)
library(sandwich)

#-------------------#
##### Load Data #####
#-------------------#

# Read data from csv files
NDVI_wil<- fread('Data/wilcoxon_test_NDVI.csv', header = TRUE, sep = ",")
temp_wil<- fread('Data/Temp/wilcoxon_test_temp.csv', header = TRUE, sep = ",")
pheno<- fread('Data/Temp/pheno_variables.csv')
data<- fread('Data/Results_filt.csv', header = TRUE, sep = ",")
NDVI<- fread('Data/NDVI.csv', header = TRUE, sep = ",")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### H1) NDVI change during frost depends on topography and phenology ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#---------------------------------------#
##### NDVI change over frost period #####
#---------------------------------------#

# The magnitude of the NDVI change during spring frost events in 2019 depend on elevation and phenological status of the plants (SOS).
# According to Bellini et al: doy 124-130 were frost events in Cembre (May 4th-May 10th)
# However lookigng at our estimates, over 90% of points were lower than average temp from May 01, May 29th so we will use this 
# as our frost definition.

frost_start <- as.Date("2019-05-01") 
frost_end <- as.Date("2019-05-31")

# Make dataframe 
NDVI<-as.data.frame(NDVI)

# Subset the dataframe to keep year 2019
NDVI_2019 <- NDVI[substr(NDVI$Date, 1, 4) == "2019", ] # Extract the first 4 characters of the date string (the year)

# Save for future use
fwrite(NDVI_2019,"Data/NDVI_2019.csv")

# Calculate ndvi_change from mean 2 weeks before and mean 2 weeks after
# Define the specific dates
desired_dates <- as.Date(c("2019-04-24", "2019-04-17", 
                           "2019-06-05", "2019-06-12"))

# Filter the dataframe for these specific dates
NDVI_filtered <- NDVI_2019[NDVI_2019$Date %in% desired_dates, ]

# Define pre-frost and post-frost dates
pre_frost_start <- as.Date("2019-04-17")
pre_frost_end <- as.Date("2019-04-24")
post_frost_start <- as.Date("2019-06-05")
post_frost_end <- as.Date("2019-06-12")

# Filter the dataframe to keep only relevant dates
NDVI_filtered <- NDVI_filtered %>%
  filter((Date >= pre_frost_start & Date <= pre_frost_end) |
           (Date >= post_frost_start & Date <= post_frost_end))

# Add a new column to categorise pre-frost and post-frost
NDVI_filtered <- NDVI_filtered %>%
  mutate(
    frost_period = case_when(
      Date >= pre_frost_start & Date <= pre_frost_end ~ "pre_frost",
      Date >= post_frost_start & Date <= post_frost_end ~ "post_frost",
      TRUE ~ NA_character_  # Handle any dates not in the specified ranges
    )
  )

# Calculate mean NDVI for each id and frost period
mean_ndvi <- NDVI_filtered %>%
  group_by(id, frost_period) %>%
  summarise(mean_ndvi = mean(NDVI, na.rm = TRUE), .groups = 'drop')

# Pivot the table to have separate columns for pre_frost and post_frost means
mean_ndvi_wide <- mean_ndvi %>%
  pivot_wider(names_from = frost_period, values_from = mean_ndvi, names_prefix = "mean_")

# View the resulting dataframe
print(mean_ndvi_wide)

# Remove NA
mean_ndvi<- na.omit(mean_ndvi_wide)

# Calculate the mean NDVI pre and post frost
mean_ndvi$ndvi_change<- mean_ndvi$mean_post_frost- mean_ndvi$mean_pre_frost

# Subset df_2019
df_2019 <- subset(data, Year == 2019)

# Merge ndvi_change and df_2019 by id
head(df_2019)
mean_ndvi<-as.data.frame(mean_ndvi)
head(mean_ndvi)
names(df_2019)[names(df_2019) == "id"] <- "ID"
names(mean_ndvi)[names(mean_ndvi) == "id"] <- "ID"
merged_df_2019 <- merge(df_2019, mean_ndvi, by = "ID", all.x = TRUE)

# View the merged dataframe
head(merged_df_2019)

# Remove rows with NA in greenup
df_filtered <- subset(merged_df_2019, !is.na(Greenup))

summary(df_filtered$Greenup)
hist(df_filtered$Greenup)
hist(df_filtered$ndvi_change)
plot(df_filtered$Greenup,df_filtered$ndvi_change)

# Save
fwrite(df_filtered,"Data/df_filtered.csv")

#-----------------------#
##### Linear model #####
#-----------------------#

quad_ndvi<-df_filtered$mean_pre_frost^2
lin_model <- lm(ndvi_change ~ mean_pre_frost + quad_ndvi + Elevation_Band + Aspect, data = df_filtered)
summary(lin_model)
vif(lin_model)

# Center quadratic to reduce VIF 
df_filtered$mean_pre_frost_c <- df_filtered$mean_pre_frost - mean(df_filtered$mean_pre_frost, na.rm = TRUE) 
df_filtered$mean_pre_frost_c2 <- df_filtered$mean_pre_frost_c^2

# Fit the centered quadratic model
quad_model_centered <- lm(ndvi_change ~ mean_pre_frost_c + mean_pre_frost_c2 +
                            Elevation_Band + Aspect, data = df_filtered)
summary(quad_model_centered)

# Save
saveRDS(quad_model_centered, file = "data/quad_model_centered.rds")

# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.45274 -0.02267  0.01235  0.03658  0.11200 
# 
# Coefficients:
#   Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)              0.2909283  0.0012233  237.813  < 2e-16 ***
#   mean_pre_frost_c        -0.9489023  0.0040064 -236.848  < 2e-16 ***
#   mean_pre_frost_c2        0.2021578  0.0287001    7.044 1.91e-12 ***
#   Elevation_Band>1400     -0.0367510  0.0014851  -24.747  < 2e-16 ***
#   Elevation_Band1000-1200 -0.0410684  0.0013091  -31.371  < 2e-16 ***
#   Elevation_Band1200-1400 -0.0513831  0.0013488  -38.097  < 2e-16 ***
#   Elevation_Band800-1000  -0.0095689  0.0013430   -7.125 1.06e-12 ***
#   Aspect                  -0.0078232  0.0005059  -15.463  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05504 on 32212 degrees of freedom
# (11459 observations deleted due to missingness)
# Multiple R-squared:  0.6844,	Adjusted R-squared:  0.6843 
# F-statistic:  9978 on 7 and 32212 DF,  p-value: < 2.2e-16

#-------------------#
# Check Assumptions #
#-------------------#

# Multicollinearity
vif(quad_model_centered) # All values are below 3
# GVIF Df GVIF^(1/(2*Df))
# mean_pre_frost_c  1.403978  1        1.184896
# mean_pre_frost_c2 1.072956  1        1.035836
# Elevation_Band    1.447344  4        1.047301
# Aspect            1.114545  1        1.055720

# Residuals and fitted values
res <- residuals(quad_model_centered)
fit <- fitted(quad_model_centered)

plot(fit, res,
     xlab = "Fitted values",
     ylab = "Residuals",
     pch = 16, cex = 0.3)
abline(h = 0, lty = 2)

qqnorm(res, pch = 16, cex = 0.3)
qqline(res, col = "red", lwd = 2)

# Histogram
hist(res, breaks = 50, main = "Residuals", xlab = "Residuals")
hist(df_filtered$ndvi_change)

# Outlier
cd <- cooks.distance(quad_model_centered)
sort(cd, decreasing = TRUE)[1:10]

# Breusch-Pagan test for heteroskedasticity
bptest(quad_model_centered) # there is statistical evidence of heteroscedascity, use heteroscedasticity-robust standard errors
# studentized Breusch-Pagan test
# 
# data:  quad_model_centered
# BP = 1125.7, df = 7, p-value < 2.2e-16

robust_q4 <- coeftest(quad_model_centered, vcov = vcovHC(quad_model_centered, type = "HC3"))
print(robust_q4)

#  t test of coefficients:
# 
# Estimate  Std. Error   t value  Pr(>|t|)    
# (Intercept)              0.29092828  0.00091210  318.9658 < 2.2e-16 ***
#   mean_pre_frost_c        -0.94890226  0.00445793 -212.8571 < 2.2e-16 ***
#   mean_pre_frost_c2        0.20215778  0.03593626    5.6255 1.866e-08 ***
#   Elevation_Band>1400     -0.03675098  0.00119536  -30.7446 < 2.2e-16 ***
#   Elevation_Band1000-1200 -0.04106835  0.00105503  -38.9264 < 2.2e-16 ***
#   Elevation_Band1200-1400 -0.05138311  0.00111168  -46.2213 < 2.2e-16 ***
#   Elevation_Band800-1000  -0.00956893  0.00089845  -10.6505 < 2.2e-16 ***
#   Aspect                  -0.00782324  0.00050816  -15.3951 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
