# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(viridis)

# Load 2024 data
file_path_2024 <- "/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/merged_data_2024.xlsx"
merged_2024 <- read_excel(file_path_2024, sheet = 1)

# Load 2022 data
file_path_2022 <- "/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/merged_data_2022.xlsx"
merged_2022 <- read_excel(file_path_2022, sheet = 1)

# Clean and process 2024 data by removing rows with NA values
data_clean_2024 <- merged_2024 %>%
  filter(!is.na(SPP), !is.na(CanopyHeight.mm.), !is.na(HitOrder), !is.na(SITE)) %>%
  mutate(Year = 2024)  # Add Year column to distinguish datasets

# Clean and process 2022 data by removing rows with NA values
data_clean_2022 <- merged_2022 %>%
  filter(!is.na(SPP), !is.na(CanopyHeight.mm.), !is.na(HitOrder), !is.na(SITE)) %>%
  mutate(Year = 2022)  # Add Year column to distinguish datasets

# Combine both datasets
combined_data <- bind_rows(data_clean_2024, data_clean_2022)

# Group by species, height, hit order, site, and year, and calculate percent cover
percent_cover_data <- combined_data %>%
  group_by(SPP, SITE, CanopyHeight.mm., HitOrder, Year) %>%
  summarise(CoverCount = n(), .groups = "drop") %>%
  group_by(SITE, Year) %>%  # Group by site and year to calculate relative percent cover
  mutate(TotalCover = sum(CoverCount)) %>%
  mutate(RelativePercentCover = (CoverCount / TotalCover) * 100) %>%
  ungroup()  # Ungroup after calculation

# Create the plot with relative percent cover and fixed y-axis scale
ggplot(percent_cover_data, aes(x = factor(CanopyHeight.mm.), y = RelativePercentCover, fill = factor(SPP))) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ SITE + Year, scales = "fixed") +  # Ensure all sites use the same y-axis scale
  labs(title = "Relative Percent Cover of Species by Height, Site Type, and Year", 
       x = "Canopy Height (mm)", 
       y = "Relative Percent Cover (%)") +
  scale_fill_viridis(discrete = TRUE) +  # A color palette suitable for many categories
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text labels
    axis.ticks.x = element_blank(),  # Remove x-axis tick marks
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

