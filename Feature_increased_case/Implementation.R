library(parallel)
library(tidyverse)

rm(list = ls())
mydir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(mydir)

source("Main_function.R")
source("Help_function.R")
library(foreach)

m1_values <- c(100, 150, 200)
m2 <- 70
p1 <- 45
lmax <- 100
num_cores <- 40

# Run code
results <- parallel_simulation(m1_values, m2, p1, lmax, num_cores = num_cores)
save(results, file = "Result/results.rda", version = 2)

# Process and analyze results
results <- list()
for (i in 1:length(m1_values)) {
  load(paste0("Result/results_", m1_values[i], ".rda"))
  results[[i]] <- result_t
}

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

##############################  Plot AUC ##############################
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

# Convert to long format
auc_long <- result_auc %>%
  pivot_longer(cols = -row,
               names_to = "method",
               values_to = "value")

p2 <- ggplot(auc_long, aes(x = row, y = value, color = method)) +
  geom_line(size = 1) +
  labs(
    title = "(C)",
    x = "Number of columns",
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
p2

ggsave("Figure/AUC_plot.pdf",
       p2,
       width = 10,
       height = 6)

)



##############################  Plot NMSE ##############################
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
               values_to = "value") %>%
  mutate(value = if_else(value > 10, NA_real_, value))

p4 <- ggplot(mse_long, aes(x = row, y = value, color = method)) +
  geom_line(size = 1) +
  scale_y_log10() +
  scale_x_continuous(breaks = m1_values) +
  labs(
    title = "(A)",
    x = "Number of columns",
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
p4

ggsave("Figure/mse_plot_log.pdf",
       p4,
       width = 10,
       height = 6)

# Group plot
library(patchwork)
combined <- p4  + p2 & theme(legend.position = "bottom")

combined + plot_layout(guides = "collect")

ggsave("Figure/Combine_figure.pdf",
       width = 10,
       height = 4)
