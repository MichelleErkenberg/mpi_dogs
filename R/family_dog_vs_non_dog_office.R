library(data.table)
library(ggplot2)
library(gridExtra)

## Set working directory
setwd("~/github/mpi_dogs/")

## Read in the txt file to define whether it is a dog office or not
dt.dog_office <- fread('data/dog_samples/R_prep/dog_env_samples_24_v1.txt', na.strings = c('-','NA',''))

## Read in the tsv file with the quicksand data for all samples
dt.tax <- fread('data/env_samples/quicksand.v2/final_report.tsv', na.strings = c('-','NA',''))

## Apply the cap to ReadsDeduped
dt.tax[, ReadsDeduped.cap := pmin(ReadsDeduped, 20000)]

## Split dog_office data
dog_office_samples <- dt.dog_office[category2 == "dog_office", sample_id]
non_dog_office_samples <- dt.dog_office[category2 != "dog_office", sample_id]

## Filter tax data for dog_office and non_dog_office samples
dt.tax_dog_office <- dt.tax[sample_id %in% dog_office_samples]
dt.tax_non_dog_office <- dt.tax[sample_id %in% non_dog_office_samples]

## Function to process data for bar plots and pie charts
process_data <- function(dt, relevant_families, threshold) {
  dt <- dt[Family %in% relevant_families & ReadsDeduped.cap > threshold]
  list(
    bar_data = dt[, .(count = uniqueN(sample_id)), by = Family],
    pie_data = dt[, .(ReadsDeduped = sum(ReadsDeduped.cap)), by = Family]
  )
}

## Set parameters
relevant_families <- c('Hominidae', 'Canidae', 'Felidae', 'Suidae')

## Process data
results_dog_office <- process_data(dt.tax_dog_office, relevant_families, threshold)
results_non_dog_office <- process_data(dt.tax_non_dog_office, relevant_families, threshold)

## Create bar plot function
create_plot <- function(data, title) {
  ggplot(data, aes(x = Family, y = count, fill = Family)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = count), vjust = -0.3) +
    theme_minimal() +
    labs(title = title, x = "Family", y = "Number of Samples") +
    theme(legend.position = "none")
}

## Create pie chart function
create_pie_chart <- function(data, title) {
  ggplot(data, aes(x = "", y = ReadsDeduped, fill = Family)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    theme_void() +
    labs(title = title) +
    geom_text(aes(label = paste0(round(ReadsDeduped/sum(ReadsDeduped)*100, 1), "%")), 
              position = position_stack(vjust = 0.5))
}

## Create plots
plot_dog_office <- create_plot(results_dog_office$bar_data, "Dog Office Samples")
plot_non_dog_office <- create_plot(results_non_dog_office$bar_data, "Non-Dog Office Samples")
pie_dog_office <- create_pie_chart(results_dog_office$pie_data, "Dog Office")
pie_non_dog_office <- create_pie_chart(results_non_dog_office$pie_data, "Non-Dog Office")

## Combine all plots
combined_plot <- grid.arrange(
  plot_dog_office, plot_non_dog_office,
  pie_dog_office, pie_non_dog_office,
  ncol = 2
)

## Save combined plot
#ggsave("figures/combined_samples_count_and_reads.png", combined_plot, width = 16, height = 16)

## Save individual plots
bar_plots <- grid.arrange(plot_dog_office, plot_non_dog_office, ncol = 2)
pie_plots <- grid.arrange(pie_dog_office, pie_non_dog_office, ncol = 2)
#ggsave("figures/combined_samples_count.png", bar_plots, width = 16, height = 8)
#ggsave("figures/combined_samples_reads.png", pie_plots, width = 16, height = 8)
