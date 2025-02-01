# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(viridis)

# Load and clean data
file_path <- "/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/merged_data_2024.xlsx"
merged_2024 <- read_excel(file_path, sheet = 1)

# Clean data by removing rows with NA values
data_clean <- merged_2024 %>%
  filter(!is.na(SPP), !is.na(CanopyHeight.mm.), !is.na(HitOrder), !is.na(SITE))

# Group by species, height, hit order, and site, and calculate percent cover
percent_cover_data <- data_clean %>%
  group_by(SPP, SITE, CanopyHeight.mm., HitOrder) %>%
  summarise(CoverCount = n(), .groups = "drop") %>%
  group_by(SITE) %>%  # Group by site to calculate relative percent cover
  mutate(TotalCover = sum(CoverCount)) %>%
  mutate(RelativePercentCover = (CoverCount / TotalCover) * 100) %>%
  ungroup()  # Ungroup after calculation

# Create the plot with relative percent cover
ggplot(percent_cover_data, aes(x = factor(CanopyHeight.mm.), y = RelativePercentCover, fill = factor(SPP))) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ SITE, scales = "free_y") +
  labs(title = "Relative Percent Cover of Species by Height and Site Type", 
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
