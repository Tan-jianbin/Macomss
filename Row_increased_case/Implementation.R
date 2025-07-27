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
source("VAE_function.R")
library(foreach)

# Basic parameter setting
m1_values <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
m2 <- 70 # Column number
p1 <- 45 # Missing block size
lmax <- 100 # Replicate number 
num_cores <- 40 # Core for parallel computation

# Run code
for(noise in c(0.2, 0.3, 0.4)){
  results <- parallel_simulation(m1_values, m2, p1, lmax, num_cores = num_cores, noise = noise)
  save(results, file = paste0("Noise", noise, "/results.rda"), version = 2)
}

