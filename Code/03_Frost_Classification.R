##%%%%%%%%%%%%%%%%%%##

### Frost Analysis ###

##%%%%%%%%%%%%%%%%%%##

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### Frost Classification ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#------------------------#
##### Load Libraries #####
#------------------------#

library(dplyr)
library(tidyr)
library(data.table)

#-------------------#
##### Load Data #####
#-------------------#

# Read data from csv files
NDVI_wil<- fread('Data/wilcoxon_test_NDVI.csv')
temp_wil<- fread('Data/Temp/wilcoxon_test_temp.csv')
NDVI <-fread('Data/NDVI.csv')
NDVI_2019 <- fread('Data/NDVI_2019.csv')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### Identify frost exposed points ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#------------------------------------------#
# Identify points exposed to frost by temp #
#------------------------------------------#

# Convert the data from wide to long format
temp_df <- temp_wil %>%
  pivot_longer(
    cols = starts_with("R_"), 
    names_to = "date", 
    values_to = "p_value"
  )

head(temp_df)

#-------------------------#
###### Process data  ######
#-------------------------#

temp_df_res <- temp_df %>%
  pivot_wider(names_from = Dir, values_from = p_value) %>%
  filter(!is.na(L) & !is.na(G)) %>%  # Remove rows where either L or G is NA
  mutate(
    G_significant = ifelse(G < 0.05, TRUE, FALSE),
    L_significant = ifelse(L < 0.05, TRUE, FALSE)
  ) %>%
  rowwise() %>%
  mutate(
    both_significant = all(c(G_significant, L_significant))
  ) %>%
  select(id, date, both_significant)

head(temp_df_res)
unique(temp_df_res$both_significant) # All Ids have only got a sig value in either G or L

# Create the new dataset
temp_df <- temp_df %>%
  pivot_wider(names_from = Dir, values_from = p_value) %>%
  filter(!is.na(L) & !is.na(G)) %>%  # Remove rows where either L or G is NA
  mutate(
    category = case_when(
      L < 0.05 & G >= 0.05 ~ "significantly lower",
      G < 0.05 & L >= 0.05 ~ "significantly higher",
      TRUE ~ "no sig dif"  # If neither is significant or both are NA
    )
  ) %>%
  select(id, date, category)

head(temp_df)
unique(temp_df$category)

temp_df <- temp_df %>%
  mutate(
    date = as.Date(sub("R_", "", date), format = "%Y%m%d")
  )

fwrite(temp_df,"Data/Long_temperature.csv")

# Now filter for May 2019
spring_2019 <- temp_df %>%
  filter(date >= as.Date("2019-05-01") & date <= as.Date("2019-05-22"))

# IDs that are consistently “significantly lower” across the whole May period are kept
frost_ids <- spring_2019 %>%
  group_by(id) %>%
  filter(all(category == "significantly lower")) %>%
  pull(id) %>% 
  unique()

frost_points<-as.data.frame(frost_ids)

# Identify points without frost - never significantly lower temp than average during may-june 
non_frost_ids <- spring_2019 %>%
  filter(category != "significantly lower") %>%
  pull(id) %>% unique()

non_frost_points<-as.data.frame(non_frost_ids)

head(frost_points)
head(non_frost_points)

#---------#
# Combine #
#---------#

# Frost IDs
frost_df <- tibble(
  ID = frost_ids,
  category = "Frost"
)

# Non-frost IDs
non_frost_df <- non_frost_points %>%
  rename(ID = non_frost_ids) %>%
  mutate(category = "Non-frost")

# Combine
all_ids <- bind_rows(frost_df, non_frost_df)

head(all_ids)
unique(all_ids$category)

length(non_frost_df$ID)
length(frost_df$ID)

#--------------------------------------------------------#
###### Identify points exposed to frost - NDVI drop ######
#--------------------------------------------------------#

# Of these points, find any where NDVI was consistently lower during frost period

# Convert NDVI data from wide to long format
ndvi_df <- NDVI_wil %>%
  pivot_longer(
    cols = starts_with("R_"), 
    names_to = "date", 
    values_to = "p_value"
  )

head(ndvi_df)

# Process data to create categories per direction
ndvi_df_res <- ndvi_df %>%
  pivot_wider(names_from = Dir, values_from = p_value) %>%
  filter(!is.na(L) & !is.na(G)) %>%
  mutate(
    # Clean up date column
    date = as.Date(gsub("R_", "", date), format = "%Y%m%d"),
    category = case_when(
      L < 0.05 & G >= 0.05 ~ "significantly lower",
      G < 0.05 & L >= 0.05 ~ "significantly higher",
      TRUE ~ "no sig dif"
    )
  ) %>%
  select(id, date, category)

head(ndvi_df_res)
unique(ndvi_df_res$category)

# Filter for 2019 spring dates (May–June)
spring_2019_ndvi <- ndvi_df_res %>%
  filter(date >= as.Date("2019-05-01") & date <= as.Date("2019-05-22"))

# IDs consistently “significantly lower” across May–June
frost_ids_ndvi <- spring_2019_ndvi %>%
  group_by(id) %>%
  filter(all(category == "significantly lower")) %>%
  pull(id) %>%
  unique()

frost_points_ndvi <- tibble(ID = frost_ids_ndvi, category = "Frost")

# IDs never significantly lower (Non-frost)
non_frost_ids_ndvi <- spring_2019_ndvi %>%
  filter(category != "significantly lower") %>%
  pull(id) %>%
  unique()

non_frost_points_ndvi <- tibble(ID = non_frost_ids_ndvi, category = "Non-frost")

# Combine frost and non-frost IDs
all_ids_ndvi <- bind_rows(frost_points_ndvi, non_frost_points_ndvi)

head(all_ids_ndvi)
head(all_ids)

#--------------------#
###### Combine  ######
#--------------------#

# Rename columns first
all_ids_temp <- all_ids %>%
  rename(Temp_class = category)

all_ids_ndvi <- all_ids_ndvi %>%
  rename(NDVI_class = category)

# Full join to keep all IDs that appear in either dataset
combined_ids <- full_join(all_ids_temp, all_ids_ndvi, by = "ID")

head(combined_ids)

# Add a column checking if they differ
combined_ids <- combined_ids %>%
  mutate(classes_differ = Temp_class != NDVI_class)

# Summary
table(combined_ids$classes_differ, useNA = "ifany")

# Inspect only mismatched IDs
diff_ids <- combined_ids %>%
  filter(classes_differ == TRUE)

head(diff_ids)

# Where temp_class is non frost but ndvi_class is frost, make a new column called Class and classify as non-frost
combined_ids <- combined_ids %>%
  mutate(Frost_Class = case_when(
    Temp_class == "Non-frost" ~ "Non-frost",
    Temp_class == "Frost" & NDVI_class == "Non-frost" ~ "Non-frost",
    TRUE ~ "Frost"
  ))

# Add this to the NDVI data
head(combined_ids)
head(NDVI_2019)

# Ensure data and all_ids have matching column names for joining
names(NDVI_2019)[names(NDVI_2019) == "id"] <- "ID"

# Filter to keep only IDs in all_ids and add the category column
ndvi_data <- NDVI_2019 %>%
  inner_join(combined_ids %>% select(ID, Frost_Class), by = "ID")

head(ndvi_data)

# Check the result
summary(ndvi_data)
table(ndvi_data$Frost_Class)

# Save
fwrite(ndvi_data,"Data/ndvi_2019_classified.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

####  Constrain data of other years to classified points ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Keep only the ids that have been classified
setDT(NDVI)
setDT(ndvi_data)

# Filter NDVI to only the IDs present in ndvi_data
NDVI_all <- NDVI[id %in% ndvi_data$ID]

head(NDVI_all)
glimpse(NDVI_all)

# Rename
NDVI_all <- NDVI_all %>% rename(ID = id)

fwrite(NDVI_all,"Data/filtered_ndvis.csv")
rm(NDVI)

# Join them together
head(NDVI_all)
head(ndvi_data)

NDVI_complete <- NDVI_all %>%
  left_join(
    ndvi_data %>%
      distinct(ID, Frost_Class),
    by = "ID"
  )

summary(NDVI_complete)

# Save
fwrite(NDVI_complete,"Data/classified_NDVI_all.csv")
