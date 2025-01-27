library(data.table)
library(ggplot2)
library(dplyr)

## Set working directory
setwd("~/github/mpi_dogs/")

## Read in the txt file with dog office and category information
dt.dog_office <- fread('data/dog_samples/R_prep/dog_env_samples_24_v1.txt', na.strings = c('-','NA',''))

## Read in the tsv file with the quicksand data for all samples
dt.tax <- fread('data/env_samples/quicksand.v2/final_report.tsv', na.strings = c('-','NA',''))

## Filter and prepare the data
dt.tax_filtered <- merge(dt.tax, dt.dog_office[, .(sample_id, category, category2, office)], by="sample_id", all.x=TRUE)
dt.tax_filtered <- dt.tax_filtered[Family %in% c('Hominidae', 'Canidae')]
dt.tax_filtered[, ReadsDeduped.cap := pmin(ReadsDeduped, 20000) + 1]

## Extract unique categories
unique_categories <- unique(dt.dog_office$category2)

## Create initial mapping
category_mapping <- setNames(
  gsub("_", " ", tools::toTitleCase(unique_categories)),
  unique_categories
)

## Print initial mapping
print("Initial category mapping:")
print(category_mapping)

## Manual adjustments
category_mapping["nc"] <- "Negative Control"
category_mapping["nc_office"] <- "No Dog Office"

## Apply the mapping to create the Category column
dt.tax_filtered[, Category := factor(category2, levels = names(category_mapping), labels = category_mapping)]

## Create a new column for combined category labels, only for "lab" category
dt.tax_filtered[, CategoryLabel := ifelse(category2 == "lab", 
                                          paste0(Category, " (", office, ")"), 
                                          as.character(Category))]

## Create a custom order for the categories
custom_order <- c(
  setdiff(unique(dt.tax_filtered$CategoryLabel), dt.tax_filtered[category2 == "lab", unique(CategoryLabel)]),
  dt.tax_filtered[category2 == "lab", sort(unique(CategoryLabel))]
)

## Ensure the order of categories is preserved with the custom order
dt.tax_filtered[, CategoryLabel := factor(CategoryLabel, levels = custom_order)]

## Create the plot
ggplot(dt.tax_filtered, aes(x = ReadsDeduped.cap, y = CategoryLabel, fill = Family)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.7) +
  scale_x_log10(labels = scales::comma) +
  theme_minimal() +
  labs(title = "Distribution of Reads by Category and Family",
       x = "mtDNA (log10 scale)",
       y = "Category") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(angle = 0, hjust = 1),
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor.y = element_blank()) 

## Save the plot
ggsave("figures/category_family_distribution.png", width = 12, height = 10)

