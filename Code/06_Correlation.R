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

#-------------------#
##### Load Data #####
#-------------------#

# Read data from csv files
data<-fread("Data/df_filtered.csv")

#-------------------------------------------------------------#
#####	Correlation between phenological metrics #####
#-------------------------------------------------------------#

data<- as.data.frame(data)
# Perform correlation on unfiltered data
colnames(data)

# Select the metrics you want to include in the correlation analysis
selected_metrics <- c("Aspect", "Elevation", "Greenup", "mean_pre_frost")  

# Subset the data frame to include only the selected metrics
subset_data <- data[, selected_metrics]

subset_data <- subset_data %>%
  rename(
    'Pre-frost NDVI' = 'mean_pre_frost',
    'Aspect(cos)' = 'Aspect',
    'Green-up' = 'Greenup'
  )

# Correlation between metrics
corr<-cor(subset_data[,unlist(lapply(subset_data, is.numeric))], use = "complete.obs",  method = "spearman")
print(corr)

# Find indices where correlation values meet the criteria
high_corr_indices <- which(corr > 0.7 | corr < -0.7, arr.ind = TRUE)

# Extract correlations greater than 0.7 or less than -0.7
high_corr_values <- corr[high_corr_indices]

# Extract row and column names for these correlations
high_corr_pairs <- rownames(corr)[high_corr_indices[, 1]] # Row names
high_corr_pairs <- cbind(high_corr_pairs, colnames(corr)[high_corr_indices[, 2]]) # Column names

# Combine names and corresponding correlation values
high_corr_results <- cbind(high_corr_pairs, high_corr_values)

# Print the results
print(high_corr_results)
