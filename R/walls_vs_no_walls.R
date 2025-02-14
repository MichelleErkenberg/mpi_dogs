library(data.table)
library(ggplot2)
library(gridExtra)
library(cowplot)

## Set working directory
setwd("~/github/mpi_dogs/")

## Read in the txt file with dog office and category information
dt.dog_office <- fread('data/dog_samples/R_prep/dog_env_samples_24_v1.txt', na.strings = c('-','NA',''))

## Read in the tsv file with the quicksand data for all samples
dt.tax <- fread('data/env_samples/quicksand.v2/final_report.tsv', na.strings = c('-','NA',''))

## Merge dt.tax with dt.dog_office to include category information
dt.tax_filtered <- merge(dt.tax, dt.dog_office[, .(sample_id, category, category2, office)], by="sample_id", all.x=TRUE)

## Filter for only Hominidae and Canidae families
relevant_families <- c('Hominidae', 'Canidae')
dt.tax_filtered <- dt.tax_filtered[Family %in% relevant_families]

## Add wall information, treating empty categories as "No Wall"
dt.tax_filtered[, `:=`(
  is_wall = ifelse(category == "wall", "Wall", "No Wall")
)]
dt.tax_filtered[is.na(is_wall) | is_wall == "", is_wall := "No Wall"]

## Define offices to exclude
offices_to_exclude <- c("NC", "Mimi/Linda Hallway", "Tracy/Silke Hallway", "Main Entrance")

## Filter out the excluded offices
dt.tax_filtered <- dt.tax_filtered[!office %in% offices_to_exclude]

## Define custom color scheme
custom_colors <- c("Canidae" = "#35978f", "Hominidae" = "#fed976")

## Function to create a single boxplot with custom colors
create_single_boxplot <- function(data, title, x_label) {
  ggplot(data, aes(x = is_wall, y = ReadsDeduped + 1, fill = Family)) +
    geom_boxplot(width = 0.7, alpha = 0.7) +
    stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "black") +
    scale_y_log10(labels = scales::comma) +
    scale_fill_manual(values = custom_colors) +
    theme_bw() +
    labs(title = title, x = x_label, y = "ReadsDeduped (log10 scale)")
}

## Creating the four separate plots
plot_dog_office_canidae <- create_single_boxplot(
  dt.tax_filtered[category2 == "dog_office" & Family == "Canidae"], 
  "Dog Office: Canidae", 
  "Sample Location"
) + theme(legend.position = "none")

plot_dog_office_hominidae <- create_single_boxplot(
  dt.tax_filtered[category2 == "dog_office" & Family == "Hominidae"], 
  "Dog Office: Hominidae", 
  "Sample Location"
) + theme(legend.position = "none")

plot_non_dog_office_canidae <- create_single_boxplot(
  dt.tax_filtered[category2 != "dog_office" & Family == "Canidae"], 
  "Non-Dog Location: Canidae", 
  "Sample Location"
) + theme(legend.position = "none")

plot_non_dog_office_hominidae <- create_single_boxplot(
  dt.tax_filtered[category2 != "dog_office" & Family == "Hominidae"], 
  "Non-Dog Location: Hominidae", 
  "Sample Location"
) + theme(legend.position = "none")

## Combine plots
combined_plots <- plot_grid(
  plot_dog_office_canidae,
  plot_non_dog_office_canidae,
  plot_dog_office_hominidae,
  plot_non_dog_office_hominidae,
  ncol = 2
)

print(combined_plots)

## Saving the combined plot
ggsave("figures/walls_vs_no_walls.png", combined_plots, width = 16, height = 16)




## T-Test
# Data walls vs not wall
wall_data <- dt.tax_filtered[is_wall == "Wall", ReadsDeduped]
no_wall_data <- dt.tax_filtered[is_wall == "No Wall", ReadsDeduped]

t_test_result <- t.test(wall_data, no_wall_data)

print(t_test_result)

