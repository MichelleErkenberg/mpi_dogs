library(data.table)
library(ggplot2)
library(gridExtra)

## Set working directory
setwd("~/github/mpi_dogs/")

## Read in the txt file with dog office and category information
dt.dog_office <- fread('data/dog_samples/R_prep/dog_env_samples_24_v1.txt', na.strings = c('-','NA',''))

## Read in the tsv file with the quicksand data for all samples
dt.tax <- fread('data/env_samples/quicksand.v2/final_report.tsv', na.strings = c('-','NA',''))

## Apply the cap to ReadsDeduped and add a small constant
dt.tax[, ReadsDeduped.cap := pmin(ReadsDeduped, 20000) + 1]

## Merge dt.tax with dt.dog_office to include category information
dt.tax_filtered <- merge(dt.tax, dt.dog_office[, .(sample_id, category, category2)], by="sample_id", all.x=TRUE)

## Filter for only Human (Hominidae) and Dog (Canidae) families
relevant_families <- c('Hominidae', 'Canidae')
dt.tax_filtered <- dt.tax_filtered[Family %in% relevant_families]

## Add wall information and rename families, treating empty categories as "No Wall"
dt.tax_filtered[, `:=`(
  is_wall = ifelse(category == "wall", "Wall", "No Wall"),
  Species = ifelse(Family == "Hominidae", "Human", "Dog")
)]
dt.tax_filtered[is.na(is_wall) | is_wall == "", is_wall := "No Wall"]

## Function to create boxplot without jitter (log scale)
create_boxplot <- function(data, title) {
  ggplot(data, aes(x = is_wall, y = ReadsDeduped.cap, fill = Species)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.7) +
    stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.8)) +
    scale_y_log10(labels = scales::comma) +
    theme_minimal() +
    labs(title = title, x = "Sample Location", y = "ReadsDeduped (capped, log10 scale)") +
    theme(legend.position = "bottom")
}

## Function to create boxplot with jitter (log scale)
create_boxplot_with_jitter <- function(data, title) {
  ggplot(data, aes(x = is_wall, y = ReadsDeduped.cap, fill = Species)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = Species), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha = 0.5) +
    stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.8)) +
    scale_y_log10(labels = scales::comma) +
    theme_minimal() +
    labs(title = title, x = "Sample Location", y = "ReadsDeduped (capped, log10 scale)") +
    theme(legend.position = "bottom")
}

## Create plots without jitter
plot_dog_office <- create_boxplot(dt.tax_filtered[category2 == "dog_office"], "Dog Office: Wall vs No Wall")
plot_non_dog_office <- create_boxplot(dt.tax_filtered[category2 != "dog_office"], "Non-Dog Office: Wall vs No Wall")

## Combine plots without jitter
combined_plot_no_jitter <- grid.arrange(plot_dog_office, plot_non_dog_office, ncol = 2)

## Save combined plot without jitter
ggsave("figures/walls_vs_no_walls.png", combined_plot_no_jitter, width = 16, height = 8)

## Create plots with jitter
plot_dog_office_jitter <- create_boxplot_with_jitter(dt.tax_filtered[category2 == "dog_office"], "Dog Office: Wall vs No Wall")
plot_non_dog_office_jitter <- create_boxplot_with_jitter(dt.tax_filtered[category2 != "dog_office"], "Non-Dog Office: Wall vs No Wall")

## Combine plots with jitter
combined_plot_with_jitter <- grid.arrange(plot_dog_office_jitter, plot_non_dog_office_jitter, ncol = 2)

## Save combined plot with jitter
ggsave("figures/walls_vs_no_walls_with_jitter.png", combined_plot_with_jitter, width = 16, height = 8)

