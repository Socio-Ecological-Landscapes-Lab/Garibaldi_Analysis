# Load necessary library
library(readxl)
library(dplyr)
library(openxlsx)

# Load 2024 point frame data and functional groups data
point_frame <- read_excel("/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/2024_point_frame.xlsx")
fgs <- read_excel("/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/Species_Func_Groups.xlsx")

# Add a column to track the original order
point_frame <- point_frame %>%
  mutate(original_order = row_number())

# Perform the join between point_frame and fgs
# Join on the correct species and subspecies columns
merged_data <- point_frame %>%
  left_join(fgs, by = c("SPP" = "SPP_code"), relationship = "many-to-many")

# Arrange the data back in its original order
merged_data <- merged_data %>%
  arrange(original_order) %>%
  select(-original_order)

# Write merged data to new sheet in 2024 point frame data 
write.xlsx(merged_data, "/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/2024_point_frame.xlsx", sheetName = "merged_data", append = TRUE)
