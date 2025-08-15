library(parallel)
library(tidyverse)

rm(list = ls())
mydir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(mydir)

source("Main_function.R")
source("Help_function.R")
source("VAE_function.R")
library(foreach)

m1_values <- c(100, 150, 200)
m2 <- 70
p1 <- 45
lmax <- 100
num_cores <- 40

# Run code
results <- parallel_simulation(m1_values, m2, p1, lmax, num_cores = num_cores)
