##%%%%%%%%%%%%%%%%%%##

### Frost Analysis ###

##%%%%%%%%%%%%%%%%%%##

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### PRELIMINARY ANALYSIS ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#------------------------#
##### Load Libraries #####
#------------------------#

library(lmerTest)
library(DHARMa)
library(car)
library(ggplot2)

#-------------------#
##### Load Data #####
#-------------------#

data<-fread("Data/Results_filt.csv")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#### Beech phenology varies with Elevation, Aspect and Year ####
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

data$Year <- as.factor(data$Year)
str(data$Year)  # should show Factor w/ levels 2018 2019 2020 2021 2022 2023

# Conduct linear mixed effect model (point is a random effect since it is a repeated measure, allowing each point to have its own intercept)
Green <- lmerTest::lmer(Greenup ~ Year + Elevation_Band + Aspect + (1 | ID), data=data)
summary(Green)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: Greenup ~ Year + Elevation_Band + Aspect + (1 | ID)
#    Data: data
# 
# REML criterion at convergence: 1968144
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -6.7581 -0.4294 -0.0028  0.4393  5.4145 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  ID       (Intercept)  66.42    8.15   
#  Residual             113.15   10.64   
# Number of obs: 251303, groups:  ID, 47531
# 
# Fixed effects:
#                           Estimate Std. Error         df  t value Pr(>|t|)    
# (Intercept)              1.090e+02  1.789e-01  5.470e+04  609.422  < 2e-16 ***
# Year2019                -7.084e-01  7.744e-02  2.059e+05   -9.148  < 2e-16 ***
# Year2020                -1.250e+01  7.753e-02  2.059e+05 -161.173  < 2e-16 ***
# Year2021                 5.926e+00  7.802e-02  2.061e+05   75.961  < 2e-16 ***
# Year2022                -2.072e+00  7.707e-02  2.058e+05  -26.879  < 2e-16 ***
# Year2023                -2.131e-01  7.810e-02  2.062e+05   -2.728  0.00637 ** 
# Elevation_Band>1400      2.291e+01  1.998e-01  4.526e+04  114.670  < 2e-16 ***
# Elevation_Band1000-1200  9.749e+00  1.873e-01  4.533e+04   52.033  < 2e-16 ***
# Elevation_Band1200-1400  1.504e+01  1.881e-01  4.516e+04   79.920  < 2e-16 ***
# Elevation_Band800-1000   4.898e+00  1.988e-01  4.547e+04   24.638  < 2e-16 ***
# Aspect                   1.720e+00  6.814e-02  4.500e+04   25.240  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) Yr2019 Yr2020 Yr2021 Yr2022 Yr2023 E_B>14 E_B100 E_B120 E_B800
# Year2019    -0.244                                                               
# Year2020    -0.243  0.562                                                        
# Year2021    -0.241  0.558  0.558                                                 
# Year2022    -0.245  0.566  0.565  0.561                                          
# Year2023    -0.242  0.558  0.557  0.553  0.561                                   
# Elvt_B>1400 -0.813  0.001  0.000  0.001  0.001  0.000                            
# E_B1000-120 -0.863  0.000 -0.001 -0.001 -0.001  0.000  0.772                     
# E_B1200-140 -0.866  0.003  0.003  0.002  0.003  0.003  0.775  0.822              
# E_B800-1000 -0.804 -0.002 -0.003 -0.003 -0.003 -0.002  0.721  0.767  0.766       
# Aspect      -0.152 -0.007 -0.003  0.000 -0.009 -0.004  0.134  0.115  0.155  0.057

#------------------------------#
##### Checking assumptions #####
#------------------------------#

# Residual distribution
qqnorm(resid(Green))
qqline(resid(Green))
hist(resid(Green), breaks=30)

# Homoscedasticity
ggplot(data.frame(fitted=fitted(Green), resid=resid(Green)), aes(x=fitted, y=resid)) +
  geom_smooth(se=FALSE, color="red") +
  geom_point(alpha=0.02) +
  labs(x="Fitted values", y="Residuals") +
  theme_minimal()

sim_res <- simulateResiduals(Green)
plot(sim_res)

# Heteroscedasticity / dispersion
testDispersion(sim_res)       

# Overall fit
testUniformity(sim_res)       

# Multicollinearity
vif(Green) # All values are below 3

# Histogram
hist(res, breaks = 50, main = "Residuals", xlab = "Residuals")
hist(data$Greenup)

# Outliers
cd <- cooks.distance(Green)
sort(cd, decreasing = TRUE)[1:10]

