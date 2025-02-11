# Load necessary library
library(readxl)
library(dplyr)
library(openxlsx)

# Load data into frame 
data <- read_excel("/Users/lizmccleary/Desktop/Garibaldi_R/2024_point_frame.xlsx")

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
write.xlsx(data, file = "/Users/lizmccleary/Desktop/Garibaldi R/2024_point_frame.xlsx", sheetName = "2024_pf_cleaned_FGs", append = TRUE)

