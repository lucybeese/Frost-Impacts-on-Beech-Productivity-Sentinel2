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
library(ggplot2)
library(emmeans)
library(lme4)
library(car)
library(lubridate)
library(scales)

#-------------------#
##### Load Data #####
#-------------------#

mod_int <- readRDS("Data/LMM.rds")
data_with_ndvi <- fread('Data/H2_data_allyear.csv')
NDVI_complete<-fread('Data/classified_NDVI_all.csv')
source(file="Scripts/wvioplot.r")

str(data_with_ndvi[, c("Year", "Frost_Class", "Elevation_Band")])

data_with_ndvi$Frost_Class <- factor(data_with_ndvi$Frost_Class)
data_with_ndvi$Elevation_Band <- factor(data_with_ndvi$Elevation_Band)

#---------------------#
##### Get EMMeans #####
#---------------------#

summary(mod_int)
emms <- emmeans(mod_int, ~ Year * Frost_Class)
summary(emms)

# Year Frost_Class emmean     SE  df asymp.LCL asymp.UCL
# 2018 Frost        0.788 0.0204 Inf     0.748     0.828
# 2019 Frost        0.760 0.0204 Inf     0.720     0.800
# 2020 Frost        0.801 0.0204 Inf     0.761     0.841
# 2021 Frost        0.789 0.0204 Inf     0.749     0.829
# 2022 Frost        0.795 0.0204 Inf     0.755     0.835
# 2023 Frost        0.803 0.0204 Inf     0.763     0.843
# 2018 Non-frost    0.783 0.0204 Inf     0.743     0.823
# 2019 Non-frost    0.782 0.0204 Inf     0.742     0.823
# 2020 Non-frost    0.794 0.0204 Inf     0.754     0.834
# 2021 Non-frost    0.783 0.0204 Inf     0.742     0.823
# 2022 Non-frost    0.789 0.0204 Inf     0.749     0.829
# 2023 Non-frost    0.798 0.0204 Inf     0.758     0.838
# 
# Results are averaged over the levels of: Elevation_Band 
# Degrees-of-freedom method: asymptotic 
# Confidence level used: 0.95 

# Test frost effect in each year H2 and H3 (within year contrast)
pairs(emms, by = "Year")

# Year = 2018:
#   contrast            estimate       SE  df z.ratio p.value
# Frost - (Non-frost)  0.00459 0.000340 Inf  13.510  <.0001
# 
# Year = 2019:
#   contrast            estimate       SE  df z.ratio p.value
# Frost - (Non-frost) -0.02246 0.000336 Inf -66.907  <.0001
# 
# Year = 2020:
#   contrast            estimate       SE  df z.ratio p.value
# Frost - (Non-frost)  0.00744 0.000336 Inf  22.163  <.0001
# 
# Year = 2021:
#   contrast            estimate       SE  df z.ratio p.value
# Frost - (Non-frost)  0.00672 0.000337 Inf  19.919  <.0001
# 
# Year = 2022:
#   contrast            estimate       SE  df z.ratio p.value
# Frost - (Non-frost)  0.00645 0.000335 Inf  19.229  <.0001
# 
# Year = 2023:
#   contrast            estimate       SE  df z.ratio p.value
# Frost - (Non-frost)  0.00449 0.000337 Inf  13.330  <.0001
# 
# Results are averaged over the levels of: Elevation_Band 
# Degrees-of-freedom method: asymptotic 

# Or test year effect separately for each frost class H3 (Within-frost year changes)
pairs(emms, by = "Frost_Class")

# Frost_Class = Frost:
#   contrast              estimate       SE  df  z.ratio p.value
# Year2018 - Year2019  2.763e-02 1.15e-04 Inf  241.029  <.0001
# Year2018 - Year2020 -1.348e-02 1.15e-04 Inf -117.464  <.0001
# Year2018 - Year2021 -1.588e-03 1.18e-04 Inf  -13.471  <.0001
# Year2018 - Year2022 -7.549e-03 1.14e-04 Inf  -66.198  <.0001
# Year2018 - Year2023 -1.526e-02 1.17e-04 Inf -130.864  <.0001
# Year2019 - Year2020 -4.111e-02 1.07e-04 Inf -384.684  <.0001
# Year2019 - Year2021 -2.922e-02 1.10e-04 Inf -264.661  <.0001
# Year2019 - Year2022 -3.518e-02 1.06e-04 Inf -331.674  <.0001
# Year2019 - Year2023 -4.289e-02 1.09e-04 Inf -393.477  <.0001
# Year2020 - Year2021  1.189e-02 1.11e-04 Inf  107.576  <.0001
# Year2020 - Year2022  5.932e-03 1.06e-04 Inf   55.798  <.0001
# Year2020 - Year2023 -1.776e-03 1.09e-04 Inf  -16.263  <.0001
# Year2021 - Year2022 -5.961e-03 1.10e-04 Inf  -54.358  <.0001
# Year2021 - Year2023 -1.367e-02 1.12e-04 Inf -121.560  <.0001
# Year2022 - Year2023 -7.707e-03 1.08e-04 Inf  -71.228  <.0001
# 
# Frost_Class = Non-frost:
#   contrast              estimate       SE  df  z.ratio p.value
# Year2018 - Year2019  5.812e-04 8.66e-05 Inf    6.711  <.0001
# Year2018 - Year2020 -1.063e-02 8.65e-05 Inf -122.985  <.0001
# Year2018 - Year2021  5.375e-04 8.82e-05 Inf    6.097  <.0001
# Year2018 - Year2022 -5.694e-03 8.60e-05 Inf  -66.219  <.0001
# Year2018 - Year2023 -1.536e-02 8.72e-05 Inf -176.264  <.0001
# Year2019 - Year2020 -1.121e-02 7.93e-05 Inf -141.385  <.0001
# Year2019 - Year2021 -4.368e-05 8.13e-05 Inf   -0.537  0.9947
# Year2019 - Year2022 -6.275e-03 7.88e-05 Inf  -79.597  <.0001
# Year2019 - Year2023 -1.594e-02 8.03e-05 Inf -198.506  <.0001
# Year2020 - Year2021  1.117e-02 8.12e-05 Inf  137.583  <.0001
# Year2020 - Year2022  4.940e-03 7.87e-05 Inf   62.727  <.0001
# Year2020 - Year2023 -4.730e-03 8.02e-05 Inf  -58.992  <.0001
# Year2021 - Year2022 -6.232e-03 8.06e-05 Inf  -77.315  <.0001
# Year2021 - Year2023 -1.590e-02 8.20e-05 Inf -193.932  <.0001
# Year2022 - Year2023 -9.669e-03 7.95e-05 Inf -121.616  <.0001