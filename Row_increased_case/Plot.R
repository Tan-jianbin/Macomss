# Set address
rm(list = ls())
# set working directory to be where the current script is located
mydir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(mydir)

# Package
library(reshape2)
library(tidyverse)

# Basic parameter setting
m1_values <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
m2 <- 70 # Column number
p1 <- 45 # Missing block size p1*p1
lmax <- 100 # Replicate number to avoid label problem
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

method_name_adjust <- factor(
  c("MACOMSS", "VAE", "VAA", "PMM", "RS", "CART", "BLR", "K-NN"),
  levels = c("MACOMSS", "VAE", "VAA", "PMM", "RS", "CART", "BLR", "K-NN")
)


#####################
# Plot
################################################
# AUC
result_auc <- data_frame()
for (i in 2:4) {
  load(paste0("Noise0.", i, "/results.rda"))
  results <- t(sapply(1:length(results), function(k) {
    sapply(1:(length(results[[k]]) - 1), function(l) {
      mean(results[[k]][[l + 1]][, 2])
    })
  }))
  colnames(results) <- method_name
  rownames(results) <- m1_values
  results <- melt(results)
  results <- data.frame(SNR = rep(paste0("1 / SNR = 0.", i), nrow(results)), results)
  result_auc <- rbind(result_auc, results)
}
colnames(result_auc) <- c("SNR", "Row", "Method", "Value")

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

# result_auc <- result_auc[result_auc[, 3] != "baseline", ]


p1 <- ggplot(result_auc, aes(x = Row, y = Value, color = Method)) +
  geom_point(aes(shape = Method), size = 2.5) +
  # geom_line() +
  geom_smooth(size = 0.6, se = F) +
  facet_wrap(~ SNR) +
  scale_x_continuous(breaks = seq(100, 1000, 200)) +
  theme_bw(base_family = "Times") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0)
  ) +
  scale_color_manual(
    values = method_colors_auc,
    labels = method_name,
    #method_name_adjust,
    name = "Method",
    guide = FALSE
  ) +
  scale_shape(guide = FALSE) +
  scale_y_sqrt(limits = c(
    min(result_auc$Value, na.rm = T),
    max(result_auc$Value, na.rm = T)
  )) +
  labs(
    title = "(B)",
    x = "Number of rows",
    y = "AUC",
    color = ""
  )
p1

##############################################
result_mse <- data_frame()
for (i in 2:4) {
  load(paste0("Noise0.", i, "/results.rda"))
  results <- t(sapply(1:length(results), function(k) {
    sapply(1:(length(results[[k]]) - 1), function(l) {
      mean(results[[k]][[l + 1]][, 4])
    })
  }))
  colnames(results) <- method_name
  rownames(results) <- m1_values
  results <- melt(results)
  results <- data.frame(SNR = rep(paste0("1 / SNR = 0.", i), nrow(results)), results)
  result_mse <- rbind(result_mse, results)
}
colnames(result_mse) <- c("SNR", "Row", "Method", "Value")

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

p2 <- ggplot(result_mse, aes(x = Row, y = Value, color = Method)) +
  geom_point(aes(shape = Method), size = 2.5) +
  # geom_line() +
  geom_smooth(size = 0.6, se = F) +
  facet_wrap(~ SNR) +
  labs(
    title = "(A)",
    x = "Number of rows",
    y = "NMSE",
    color = ""
  ) +
  scale_x_continuous(breaks = seq(100, 1000, 200)) +
  theme_bw(base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0)
  ) +
  scale_color_manual(
    values = method_colors_nmse,
    labels = method_name_adjust,
    name = "Method",
    guide = FALSE
  ) +
  scale_shape(guide = FALSE) +
  scale_y_log10(limits = c(min(result_mse$Value, na.rm = T), 5))
p2

# Group plot
library(patchwork)
(p2 + p1) +
  scale_color_manual(
    values = method_colors_auc,
    labels = method_name,
    name   = "Method"
  ) +
  guides(
    color = guide_legend(
      nrow  = 1,
      byrow = TRUE
    )
  ) +
  plot_layout(
    guides = "collect",
    ncol = 1
  ) &
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal"
  )

ggsave("Combine_figure.pdf", width = 9, height = 6)
