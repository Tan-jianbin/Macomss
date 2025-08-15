# Package
library(parallel)
library(tidyverse)

# Set address
rm(list = ls())
# set working directory to be where the current script is located
mydir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(mydir)

# Load Help_function and package
source("Main_function.R")
source("Help_function.R")

# Basic parameter setting
m1_values <- c(100, 150, 200)
m2 <- 70 # Column number
p1 <- 45 # Missing block size p1*p1
num_cores <- 40 # Core for parallel computation
method_name <- factor(
  c(
    "Complete-data benchmark",
    "MACOMSS",
    "VAE",
    "VAA",
    "PMM",
    "RS",
    "CART",
    "BLR",
    "K-NN"
  ),
  levels = c(
    "Complete-data benchmark",
    "MACOMSS",
    "VAE",
    "VAA",
    "PMM",
    "RS",
    "CART",
    "BLR",
    "K-NN"
  )
)

#####################
# Plot
################################################
# AUC
# Process and analyze results
results <- list()
for (i in 1:length(m1_values)) {
  load(paste0("Result/results_", m1_values[i], ".rda"))
  results[[i]] <- result_t
}

################################################
# AUC
result_auc <- sapply(1:length(results), function(k) {
  sapply(1:(length(results[[k]]) - 1), function(l) {
    mean(results[[k]][[l + 1]][, 2])
  })
})

result_auc <- t(result_auc)
rownames(result_auc) <- m1_values
result_auc <- data.frame(row = m1_values, result_auc)
result_auc <- data.frame(result_auc)
colnames(result_auc)[-1] <- c("Complete-data benchmark",
                              "MACOMSS",
                              "VAE",
                              "VAA",
                              "PMM",
                              "RS",
                              "CART",
                              "BLR",
                              "K-NN")

auc_long <- result_auc %>%
  pivot_longer(cols = -row,
               names_to = "Method",
               values_to = "value")
auc_long$Method <- factor(
  auc_long$Method,
  levels = c(
    "Complete-data benchmark",
    "MACOMSS",
    "VAE",
    "VAA",
    "PMM",
    "RS",
    "CART",
    "BLR",
    "K-NN"
  )
)


method_colors_auc <- c(
  "Complete-data benchmark" = "purple",
  "MACOMSS" = "#e31a1c",
  "VAE" = "cyan3",
  "VAA" = "#838B8B",
  "PMM" = "#fdbf6f",
  "RS" = '#fb9a99',
  "CART" = '#33a02c',
  "BLR" = '#a6761d',
  "K-NN" = "#1f78b4"
)

p1 <- ggplot(auc_long, aes(x = row, y = value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "(B)",
       x = "Number of columns",
       y = "AUC",
       fill = "") +
  scale_x_continuous(breaks = m1_values) +
  theme_bw(base_family = "Times") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0)
  ) +
  coord_cartesian(ylim = c(0.92, 0.98)) +
  scale_fill_manual(values = method_colors_auc)
p1

##############################################
result_mse <- sapply(1:length(results), function(k) {
  sapply(1:(length(results[[k]]) - 1), function(l) {
    mean(results[[k]][[l + 1]][, 4])
  })
})

result_mse <- t(result_mse)
rownames(result_mse) <- m1_values
result_mse <- data.frame(row = m1_values, result_mse)
result_mse <- data.frame(result_mse)
result_mse[,9] <- 3
colnames(result_mse)[-1] <- c("Complete-data benchmark",
                              "MACOMSS",
                              "VAE",
                              "VAA",
                              "PMM",
                              "RS",
                              "CART",
                              "BLR",
                              "K-NN")

mse_long <- result_mse %>%
  pivot_longer(cols = -row,
               names_to = "Method",
               values_to = "value")
mse_long$Method <- factor(
  mse_long$Method,
  levels = c(
    "Complete-data benchmark",
    "MACOMSS",
    "VAE",
    "VAA",
    "PMM",
    "RS",
    "CART",
    "BLR",
    "K-NN"
  )
)

method_colors_nmse <- c(
  "MACOMSS" = "#e31a1c",
  "VAE" = "cyan3",
  "VAA" = "#838B8B",
  "PMM" = "#fdbf6f",
  "RS" = '#fb9a99',
  "CART" = '#33a02c',
  "BLR" = '#a6761d',
  "K-NN" = "#1f78b4"
)

p2 <- ggplot(mse_long, aes(x = row, y = value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(breaks = m1_values) +
  labs(title = "(A)",
       x = "Number of columns",
       y = "NMSE",
       fill = "") +
  theme_bw(base_family = "Times") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0)
  ) +
  coord_cartesian(ylim = c(0, 3)) +
  scale_fill_manual(values = method_colors_nmse)
p2

# Group plot
library(patchwork)

combined <- (p2 + theme(legend.position = "none")) + p1

combined + plot_layout(guides = "collect") 

ggsave(
  "Combine_figure_col.pdf",
  width = 12,
  height = 4,
  dpi = 300
)
