install.packages(c("reshape2", "gridExtra"))
library(ggplot2)
library(reshape2)
library(gridExtra)

setwd('/github/mpi_dogs/R/')

data <- read.csv("Container.csv")

# Heatmap for each Container dog
heatmap_heidi <- ggplot(data, aes(x = x, y = y, fill = Heidi)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
  scale_x_continuous(breaks = unique(data$X)) +
  scale_y_continuous(breaks = unique(data$Y)) +
  coord_fixed() +
  labs(title = "Heatmap for Heidi") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),  
    panel.grid = element_blank(),
    axis.text = element_blank()    
  )


heatmap_vito <- ggplot(data, aes(x = x, y = y, fill = Vito)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
  scale_x_continuous(breaks = unique(data$X)) +
  scale_y_continuous(breaks = unique(data$Y)) +
  coord_fixed() +
  labs(title = "Heatmap for Vito") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),  
    panel.grid = element_blank(),
    axis.text = element_blank()    
  )

heatmap_fritzy <- ggplot(data, aes(x = x, y = y, fill = Fritzy)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
  scale_x_continuous(breaks = unique(data$X)) +
  scale_y_continuous(breaks = unique(data$Y)) +
  coord_fixed() +
  labs(title = "Heatmap for Fritzy") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),  
    panel.grid = element_blank(),
    axis.text = element_blank()    
  )

heatmap_urza <- ggplot(data, aes(x = x, y = y, fill = Urza)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
  scale_x_continuous(breaks = unique(data$X)) +
  scale_y_continuous(breaks = unique(data$Y)) +
  coord_fixed() +
  labs(title = "Heatmap for Urza") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),  
    panel.grid = element_blank(),
    axis.text = element_blank()    
  )



print(heatmap_heidi)
print(heatmap_vito)
print(heatmap_fritzy)
print(heatmap_urza)
grid.arrange(heatmap_fritzy, heatmap_heidi, heatmap_urza, heatmap_vito, ncol = 2)

