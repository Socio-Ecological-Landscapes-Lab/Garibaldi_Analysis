# Load necessary library
library(readxl)
library(dplyr)
library(openxlsx)

# Load data into frame 
data <- read_excel("/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/2024_point_frame_raw.xlsx")

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
data$SPP <- plot_naming(data$SPP, "moss.", "moss")

data$SPP_sub <- plot_naming(data$SPP_sub, "carmini.*", "poaspp")
data$SPP_sub <- plot_naming(data$SPP_sub, "soli", "eriper")
data$SPP_sub <- plot_naming(data$SPP_sub, "carmed", "carnig")
data$SPP_sub <- plot_naming(data$SPP_sub, "moss.", "moss")


data$SPP <- data$SPP_sub

# Write changes to the dataset
write.xlsx(data, file = "/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/2024_point_frame_raw.xlsx", sheetName = "2024_pf_cleaned_FGs", append = TRUE)

#Now that the data is clean, we can add functional groups 
# Load 2022 point frame data and functional groups data
point_frame2022 <- read_excel("/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/2022_point_frame_raw.xlsx")
fgs <- read_excel("/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/Species_Func_Groups.xlsx")

# Add a column to track the original order
point_frame2022 <- point_frame2022 %>%
  dplyr::mutate(original_order = row_number())

# Alter the Canopy Height column so that it is the same measurement of the 2024 data 
point_frame2022 <- point_frame2022 %>%
  dplyr::mutate(CanopyHeight.mm. = CanopyHeight.mm. * 10) 

# Make sure DATE is character
point_frame2022$DATE <- as.character(point_frame2022$DATE)

# Replace "6/22" and similar with "2022"
point_frame2022$DATE[point_frame2022$DATE == "6/22"] <- "2022"

# Replace "9/26/22" with "2022"
point_frame2022$DATE[grepl("^\\d{1,2}/\\d{1,2}/\\d{2}$", point_frame2022$DATE)] <- "2022"

# Replace any "2019" with "2022"
point_frame2022$DATE <- plot_naming(point_frame2022$DATE, "2019", "2022")

# Keep last four digits of any unchanged DATE fields 
point_frame2022$DATE <- gsub(".*(\\d{4})$", "\\1", point_frame2022$DATE)

# Perform the join between point_frame and fgs
merged_data2022 <- point_frame2022 %>%
  left_join(fgs, by = c("SPP" = "SPP_code"), relationship = "many-to-many")

# Arrange the data back in its original order
merged_data2022 <- merged_data2022 %>%
  arrange(original_order) %>%
  select(-original_order)

# Load 2024 point frame data and functional groups data
point_frame2024 <- read_excel("/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/2024_point_frame_raw.xlsx")
fgs <- read_excel("/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/Species_Func_Groups.xlsx")

# Add a column to track the original order
point_frame2024 <- point_frame2024 %>%
  dplyr::mutate(original_order = row_number())
point_frame2024 <- point_frame2024 %>%
  dplyr::mutate(DATE = substr(DATE, nchar(DATE) - 3, nchar(DATE)))


# Perform the join between point_frame and fgs
merged_data2024 <- point_frame2024 %>%
  left_join(fgs, by = c("SPP" = "SPP_code"), relationship = "many-to-many")

# Arrange the data back in its original order
merged_data2024 <- merged_data2024 %>%
  arrange(original_order) %>%
  select(-original_order)

# Append merged_data2024 to merged_data2022 and remove unnecessary columns 
combined_data <- bind_rows(merged_data2022, merged_data2024)%>%
  dplyr::select(-filename, -...1, -SPP_sub, -TOMST, -Flag, -Flag.1, -Species, -Validated_23_08, -Description)

write.xlsx(combined_data, "/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/ITEX_2025.xlsx", sheetName = "ITEX_2025")

