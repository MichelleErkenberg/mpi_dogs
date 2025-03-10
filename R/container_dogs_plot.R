#install.packages(c("reshape2", "gridExtra"))
library(ggplot2)
library(reshape2)
library(gridExtra)

setwd("~/github/mpi_dogs/")

data <- read.csv("data/dog_samples/R_prep/all_dogs_AC/Container.csv")

# Heatmap for each Container dog
heatmap_heidi <- ggplot(data, aes(x = x, y = y, fill = AC.Heidi)) +
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


heatmap_vito <- ggplot(data, aes(x = x, y = y, fill = AC.Vito)) +
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

heatmap_fritzy <- ggplot(data, aes(x = x, y = y, fill = AC.Fritzy)) +
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

heatmap_urza <- ggplot(data, aes(x = x, y = y, fill = AC.Urza)) +
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


# creating plots for Cami as our ref dog (in Container and Mimi/Linda office)
heatmap_cami_c <- ggplot(data, aes(x = x, y = y, fill = AC.Cami)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
  scale_x_continuous(breaks = unique(data$X)) +
  scale_y_continuous(breaks = unique(data$Y)) +
  coord_fixed() +
  labs(title = "Heatmap for Cami (Container)") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),  
    panel.grid = element_blank(),
    axis.text = element_blank()    
  )

print(heatmap_cami_c)

data_lily <- read.csv("data/dog_samples/R_prep/all_dogs_AC/Mimi_Linda.csv")

heatmap_lily_l <- ggplot(data_lily, aes(x = x, y = y, fill = AC.Lily)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
  scale_x_continuous(breaks = unique(data$X)) +
  scale_y_continuous(breaks = unique(data$Y)) +
  coord_fixed() +
  labs(title = "Heatmap for Cami (Mimi/Linda)") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),  
    panel.grid = element_blank(),
    axis.text = element_blank()    
  )

print(heatmap_lily_l)

grid.arrange(heatmap_cami_c, heatmap_cami_l, ncol = 2)
