library(data.table)
library(ggplot2)

# Set working directory
setwd("~/github/mpi_dogs/")

# Read the main CSV file with dog data
dt.main <- fread('data/dog_samples/R_prep/all_dogs_ACwoL/R_prep_sample_vs_dog_ACwoL_renamed.csv', na.strings = c('-','NA',''))

# Read the TXT file with location categories
dt.categories <- fread('data/dog_samples/R_prep/sample_location.txt', na.strings = c('-','NA',''))

# Define custom order and custom colors
custom_colors <- c()  # Add your custom colors here if needed
custom_order <- c()  # Leave this empty if you want to sort by sample ID number

# Merge main data with categories
dt.combined <- merge(dt.main, dt.categories, by = "sample_id", all.x = TRUE)

# Identify dog columns (all columns except sample_id and location)
dog_columns <- setdiff(names(dt.combined), c("sample_id", "location"))

# Reshape data to long format
dt.long <- melt(dt.combined, id.vars = c("sample_id", "location"), 
                variable.name = "dog", value.name = "value")

# Replace NA values with 0
dt.long[is.na(value), value := 0]

# Extract number from sample_id
dt.long[, sample_number := as.numeric(sub(".*_", "", sample_id))]

# Sort data by location and sample_number
setorder(dt.long, location, sample_number)

# Apply custom order to location factor if provided, otherwise sort by sample number
if(length(custom_order) > 0) {
  dt.long$location <- factor(dt.long$location, levels = custom_order)
} else {
  dt.long$location <- factor(dt.long$location, levels = unique(dt.long$location[order(dt.long$sample_number)]))
}

# Create the bar plot
p <- ggplot(dt.long, aes(x = factor(sample_number), y = value, fill = dog)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(dog ~ location, scales = "free_x", space = "free_x") +
  theme_bw() +
  labs(title = "Dog Values by Sample and Location",
       x = "Sample Number", y = "Value", fill = "Dog") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "none")

# Apply custom colors if defined
if(length(custom_colors) > 0) {
  p <- p + scale_fill_manual(values = custom_colors)
}

print(p)

# Save the plot
ggsave("figures/dogs_categorized.png", p, width = 24, height = length(dog_columns) * 2, limitsize = FALSE)
