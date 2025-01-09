install.packages(c("reshape2", "gridExtra"))
library(ggplot2)
library(reshape2)
library(gridExtra)

setwd('/github/mpi_dogs/R/')

data <- read.csv("Mimi_Linda.csv")

# Heatmap for each dog
heatmap_lily <- ggplot(data, aes(x = x, y = y, fill = Lily)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
  scale_x_continuous(breaks = unique(data$X)) +
  scale_y_continuous(breaks = unique(data$Y)) +
  coord_fixed() +
  labs(title = "Heatmap for Lily") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),  
    panel.grid = element_blank(),
    axis.text = element_blank()    
  )


heatmap_thorA <- ggplot(data, aes(x = x, y = y, fill = ThorA)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
  scale_x_continuous(breaks = unique(data$X)) +
  scale_y_continuous(breaks = unique(data$Y)) +
  coord_fixed() +
  labs(title = "Heatmap for ThorA") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),  
    panel.grid = element_blank(),
    axis.text = element_blank()    
  )



print(heatmap_lily)
print(heatmap_thorA)
grid.arrange(heatmap_lily, heatmap_thorA, ncol = 1)

