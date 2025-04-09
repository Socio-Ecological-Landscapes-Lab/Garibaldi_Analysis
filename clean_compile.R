# Load necessary library
library(readxl)
library(dplyr)
library(openxlsx)

# Load data into frame 
data <- read_excel("/Users/lizmccleary/Desktop/Garibaldi_R/2024_point_frame_raw.xlsx")

plot_naming <- function(data_column, orig_name1, new_name1){
  data_column <- as.character(data_column)
  data_column[grep(orig_name1, data_column)] <- new_name1
  data_column
}

# Declare species column as column of interest
colnames(data) 

unique(data$SPP)
unique(data$SPP_sub)

# Replace species names that start with "carmini" with "poaspp"
data$SPP <- plot_naming(data$SPP, "carmini.*", "poaspp")
data$SPP <- plot_naming(data$SPP, "soli", "eriper")
data$SPP <- plot_naming(data$SPP, "carmed", "carnig")

data$SPP_sub <- plot_naming(data$SPP_sub, "carmini.*", "poaspp")
data$SPP_sub <- plot_naming(data$SPP_sub, "soli", "eriper")
data$SPP_sub <- plot_naming(data$SPP_sub, "carmed", "carnig")


data$SPP <- data$SPP_sub

# Write changes to the dataset
# write.xlsx(data, file = "/Users/lizmccleary/Desktop/Garibaldi R/2024_point_frame_raw.xlsx", sheetName = "2024_pf_cleaned_FGs", append = TRUE)

#Now that the data is clean, we can add functional groups 
# Load 2022 point frame data and functional groups data
point_frame <- read_excel("/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/2022_point_frame_raw.xlsx")
fgs <- read_excel("/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/Species_Func_Groups.xlsx")

# Add a column to track the original order
point_frame <- point_frame %>%
  mutate(original_order = row_number())

# Alter the Canopy Height column so that it is the same measurement of the 2024 data 
point_frame <- point_frame %>%
  mutate(CanopyHeight.mm. = CanopyHeight.mm. * 10) 

# Perform the join between point_frame and fgs
# Join on the correct species and subspecies columns
merged_data2022 <- point_frame %>%
  left_join(fgs, by = c("SPP" = "SPP_code"), relationship = "many-to-many")

# Arrange the data back in its original order
merged_data2022 <- merged_data2022 %>%
  arrange(original_order) %>%
  select(-original_order)

# Write merged data to new sheet in 2022 point frame data 
# write.xlsx(merged_data, "/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/merged_data_2022.xlsx", sheetName = "merged_data_2022")


# Load 2024 point frame data and functional groups data
point_frame <- read_excel("/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/2024_point_frame_raw.xlsx")
fgs <- read_excel("/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/Species_Func_Groups.xlsx")

# Add a column to track the original order
point_frame <- point_frame %>%
  mutate(original_order = row_number())

# Perform the join between point_frame and fgs
# Join on the correct species and subspecies columns
merged_data2024 <- point_frame %>%
  left_join(fgs, by = c("SPP" = "SPP_code"), relationship = "many-to-many")

# Arrange the data back in its original order
merged_data2024 <- merged_data2024 %>%
  arrange(original_order) %>%
  select(-original_order)

#Append merged_data2024 to bottom of merged_data2022 and write to one dataset called ITEX_2025
write.xlsx(merged_data, "/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/ITEX_2025.xlsx", sheetName = "ITEX_2025")

