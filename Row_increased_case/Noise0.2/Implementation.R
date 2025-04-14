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
library(foreach)

# Basic parameter setting
m1_values <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
m2 <- 70 # Column number
p1 <- 45 # Missing block size p1*p1
lmax <- 100 # Replicate number to avoid label probalem
num_cores <- 40 # Core for parallel computation

# Run code
results <- parallel_simulation(m1_values, m2, p1, lmax, num_cores = num_cores)
save(results, file = "Result/results.rda", version = 2)

# Process and analyze results
load("Result/results.rda")
for (r in results) {
  cat("Results for m1 =", r$m1, ":\n")
  for (method in c("raw",
                   "MACOMSS",
                   "mice",
                   "mice_s",
                   "mice_cart",
                   "mice_norm",
                   "fill_KNN")) {
    cat(method, "method:\n")
    print(colMeans(r[[paste0("error_", method)]], na.rm = TRUE))
  }
  cat("\n")
}

####################################################################################
# Plot

####################################################################################
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
colnames(result_auc) <- method_name

# Convert into long data
auc_long <- result_auc %>%
  pivot_longer(cols = -row,
               names_to = "method",
               values_to = "value")

p1 <- ggplot(auc_long, aes(x = row, y = value, color = method)) +
  geom_line(size = 1) +
  labs(
    title = "(C)",
    x = "Number of rows",
    y = "AUC",
    color = ""
  ) +
  scale_x_continuous(breaks = m1_values) +
  theme_bw(base_family = "Times") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0)
  ) +
  scale_color_manual(values = c(
    '#1f78b4',
    '#e31a1c',
    '#fdbf6f',
    'purple',
    '#33a02c',
    '#a6761d',
    '#fb9a99'
  )) +
  scale_y_log10()
p1

ggsave("Figure/AUC_plot.pdf",
       p1,
       width = 10,
       height = 6)

############################################################################################
# MSE
result_mse <- sapply(1:length(results), function(k) {
  sapply(1:(length(results[[k]]) - 1), function(l) {
    mean(results[[k]][[l + 1]][, 4])
  })
})

result_mse <- t(result_mse)
rownames(result_mse) <- m1_values
result_mse <- data.frame(row = m1_values, result_mse)
result_mse <- data.frame(result_mse)
colnames(result_mse) <- method_name


mse_long <- result_mse %>%
  pivot_longer(cols = -row,
               names_to = "method",
               values_to = "value")

mse_long <- mse_long %>%
  mutate(value = ifelse(value > 10, NA, value))

p2 <- ggplot(mse_long, aes(x = row, y = value, color = method)) +
  geom_line(size = 1) +
  scale_x_continuous(breaks = m1_values) +
  labs(
    title = "(A)",
    x = "Number of rows",
    y = "NMSE",
    color = ""
  ) +
  theme_bw(base_family = "Times") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0)
  ) +
  scale_color_manual(values = c(
    '#1f78b4',
    '#e31a1c',
    '#fdbf6f',
    'purple',
    '#33a02c',
    '#a6761d',
    '#fb9a99'
  ))
p2

ggsave("Figure/mse_plot_log.pdf",
       p2,
       width = 10,
       height = 6)

# Group plot
library(patchwork)
combined <- p2 + p1 + p1 & theme(legend.position = "bottom")

combined + plot_layout(guides = "collect")

ggsave("Figure/Combine_figure.pdf",
       width = 10,
       height = 4)
