# Load required libraries
# library(parallel)
library(matrixcalc)
library(MASS)
library(pROC)
library(glmnet)
library(mice)
library(filling)
library(VIM)

fit_and_calculate_errors <- function(data,
                                     method_name,
                                     train_label,
                                     test_label,
                                     beta,
                                     test_data,
                                     train_data,
                                     train_data_bar,
                                     XB) {
  train_label_noise <- train_label[[1]]
  mask_missing_indicator <- train_label[[2]]
  
  cv_fit <- cv.glmnet(
    data,
    train_label_noise,
    family = "binomial",
    alpha = 0.25,
    intercept = TRUE,
    standardize = TRUE
  )
  best_lambda <- cv_fit$lambda.min
  lrm <- glmnet(
    data,
    train_label_noise,
    family = "binomial",
    alpha = 0.25,
    lambda = best_lambda,
    intercept = TRUE,
    standardize = TRUE
  )
  
  
  beta_method <- as.matrix(coef(lrm, s = "lambda.min"))
  
  beta_method[is.na(beta_method)] <- 0
  beta_method <- matrix(beta_method, ncol = 1)
  
  design_mat <- train_data_bar[[2]]
  XB_bar <- c(design_mat %*% beta_method)
  
  beta_diff <- sum((beta - beta_method) ^ 2, na.rm = TRUE) / sum(beta ^
                                                                   2, na.rm = TRUE)
  zero_mark <- beta == 0
  zero_mark_estimated <- beta_method == 0
  overlap <- sum(zero_mark & zero_mark_estimated)
  
  predictions <- as.vector(predict(lrm, newx = test_data, type = "response"))
  roc_data <- data.frame(test_label = test_label, predictions = predictions)
  g <- roc(test_label ~ predictions, data = roc_data)
  
  
  mae <- ifelse(
    method_name != "raw" && method_name != "omit",
    mean(abs(data - train_data), na.rm = TRUE) / mean(abs(data), na.rm = TRUE),
    NA
  )
  
  data[mask_missing_indicator] <- NA
  true_data <- train_data_bar[[3]]
  true_data[mask_missing_indicator] <- NA
  
  mse <- ifelse(
    method_name != "raw" && method_name != "omit",
    mean((data - true_data) ^ 2, na.rm = TRUE) / mean(true_data ^
                                                        2, na.rm = TRUE),
    NA
  )
  
  c(beta_diff, g$auc, mae, mse)
}



implement_func <- function(l, m1, m2, p1) {
  set.seed(l * 100 * m1 + m2 * p1 + 300) # Set a seed for simulation, then your result is reproduciable.
  
  repeat {
    x <- sample_simulate_1(m1 + 200, m2)
    xx <- x[[1]]
    xxx <- x[[4]]  # Y
    fixed_label <- x[[5]]
    
    
    x_nolabel <- xx[, -ncol(xx)]
    
    label_list <- x
    labeled_x <- label_list[[1]]
    beta <- as.matrix(label_list[[2]])
    label_raw <- as.matrix(label_list[[3]])
    
    total_rows <- nrow(labeled_x)
    indices <- sample(1:total_rows)
    
    # Define split point
    split_point <- floor(total_rows - 200)
    
    # Train-test split
    train_indices <- indices[1:split_point]
    test_indices <- indices[(split_point + 1):total_rows]
    
    train_all <- labeled_x[train_indices, ]
    train_data <- as.matrix(train_all[, -ncol(train_all)])
    train_label <- as.matrix(train_all[, ncol(train_all)])
    
    designm <- x[[7]]
    XB <- x[[6]]
    true_data <- x[[8]]
    
    train_true_data <- true_data[train_indices, ]
    
    ones_column <- rep(1, nrow(train_data))
    train_data_bar <- cbind(ones_column, train_data)
    train_data_bar <- as.matrix(train_data_bar)
    train_data_bar <- list(train_data_bar, designm, train_true_data) #### list here
    
    test_all <- labeled_x[test_indices, ]
    test_label <- as.matrix(fixed_label[test_indices, ])
    
    # Add noise
    train_data_noise <- xxx[train_indices, ]
    test_data <- xxx[test_indices, ]
    test_data_noise <- xxx[test_indices, ]
    
    # Check if train_label or test_label contains only 0 or only 1
    if (length(unique(train_label)) > 1 &&
        length(unique(test_label)) > 1) {
      break
    }
  }
  
  ####################data process
  p2 = m1 * 0.5  # set missing scale
  # MACOMSS method
  train_missing <- missing_s(train_data_noise)
  n_r <- nrow(train_missing)
  n_c <- ncol(train_missing)
  Y11 = as.matrix(train_missing[1:(n_r - p2), 1:(n_c - p1)])
  Y12 = as.matrix(train_missing[1:(n_r - p2), (n_c - p1 + 1):n_c])
  Y21 = as.matrix(train_missing[(n_r - p2 + 1):n_r, 1:(n_c - p1)])
  train_missing[(n_r - p2 + 1):n_r, (n_c - p1 + 1):n_c] <- NA
  train_missing_mat <- as.matrix(train_missing)
  mask_missing_indicator <- !is.na(train_missing)
  train_label <- list(train_label, mask_missing_indicator)
  
  # Raw method
  error_raw <- fit_and_calculate_errors(
    train_data,
    "raw",
    train_label,
    test_label,
    beta,
    test_data,
    train_data,
    train_data_bar,
    XB
  )
  x_MACOMSS <- MACOMSS.parsvd(Y11, Y12, Y21)[[1]]
  error_MACOMSS <- fit_and_calculate_errors(
    x_MACOMSS,
    "MACOMSS",
    train_label,
    test_label,
    beta,
    test_data_noise ,
    train_data_noise,
    train_data_bar,
    XB
  )
  print("MACOMSS end")
  
  # MICE method pmm
  train_missing[(n_r - p2 + 1):n_r, (n_c - p1 + 1):n_c] <- NA
  result_mice <- mice(
    train_missing,
    method = "pmm",
    m = 1,
    maxit = 1,
    printFlag = F
  )
  complete_mice <- as.matrix(mice::complete(result_mice))
  complete_mice <- apply(complete_mice, 2, function(x) {
    ifelse(is.na(x), mean(x, na.rm = TRUE), x)
  })
  error_mice <- fit_and_calculate_errors(
    complete_mice,
    "mice",
    train_label,
    test_label,
    beta,
    test_data_noise ,
    train_data_noise,
    train_data_bar,
    XB
  )
  print("Mice pmm end")
  
  # MICE method sample
  result_mice_s <- mice(
    train_missing,
    method = "sample",
    m = 1,
    maxit = 1,
    printFlag = F
  )
  complete_mice_s <- as.matrix(mice::complete(result_mice_s))
  complete_mice_s <- apply(complete_mice_s, 2, function(x) {
    ifelse(is.na(x), mean(x, na.rm = TRUE), x)
  })
  error_mice_s <- fit_and_calculate_errors(
    complete_mice_s,
    "mice_sample",
    train_label,
    test_label,
    beta,
    test_data_noise,
    train_data_noise,
    train_data_bar,
    XB
  )
  print("MICE sample end")
  # MICE method cart
  result_mice_cart <- mice(
    train_missing,
    method = "cart",
    m = 1,
    maxit = 1,
    printFlag = F
  )
  complete_mice_cart <- as.matrix(mice::complete(result_mice_cart))
  complete_mice_cart <- apply(complete_mice_cart, 2, function(x) {
    ifelse(is.na(x), mean(x, na.rm = TRUE), x)
  })
  error_mice_cart <- fit_and_calculate_errors(
    complete_mice_cart,
    "mice_cart",
    train_label,
    test_label,
    beta,
    test_data_noise ,
    train_data_noise,
    train_data_bar,
    XB
  )
  print("MICE cart end")
  
  # MICE method norm
  result_mice_norm <- mice(
    train_missing,
    method = "norm",
    m = 1,
    maxit = 1,
    printFlag = F
  )
  complete_mice_norm <- as.matrix(mice::complete(result_mice_norm))
  complete_mice_norm <- apply(complete_mice_norm, 2, function(x) {
    ifelse(is.na(x), mean(x, na.rm = TRUE), x)
  })
  error_mice_norm <- fit_and_calculate_errors(
    complete_mice_norm,
    "mice_norm",
    train_label,
    test_label,
    beta,
    test_data_noise,
    train_data_noise,
    train_data_bar,
    XB
  )
  print("MICE norm end")
  
  # KNN
  k_col <- ncol(train_missing)
  # dat_fill_1 <- fill.KNNimpute(train_missing,5)[[1]]
  dat_fill_1 <- fill.KNNimpute(train_missing, 50)[[1]]
  
  dat_fill_2 <- apply(dat_fill_1, 2, function(x) {
    ifelse(is.na(x), mean(x, na.rm = TRUE), x)
  }) #column mean to impute
  
  dat_fill <- as.matrix(dat_fill_2)
  error_fill_KNN <- fit_and_calculate_errors(
    dat_fill,
    "fill_KNN",
    train_label,
    test_label,
    beta,
    test_data_noise ,
    train_data_noise,
    train_data_bar,
    XB
  )
  print("fill KNN end")
  
  
  return(
    list(
      error_raw = error_raw,
      error_MACOMSS = error_MACOMSS,
      error_mice = error_mice,
      error_mice_s = error_mice_s,
      error_mice_cart = error_mice_cart,
      error_mice_norm = error_mice_norm,
      error_fill_KNN = error_fill_KNN
    )
  )
}

parallel_simulation <- function(m1_values, m2, p1, lmax, num_cores = detectCores() - 1) {
  # Create cluster
  cl <- parallel::makeCluster(num_cores)
  invisible(gc())
  doParallel::registerDoParallel(cl)
  
  clusterExport(
    cl,
    c(
      "implement_func",
      "sample_simulate_1",
      "missing_s",
      "fit_and_calculate_errors",
      "MACOMSS.parsvd",
      "MACOMSS",
      "sample_simulate",
      "missing_simulate",
      "label_binary",
      "sigmoid_binary_y"
    )
  )
  
  # Main simulation function
  results <- list()
  
  for (k in 1:length(m1_values)) {
    m1 <- m1_values[k]
    m1_results <- foreach(
      l = 1:lmax,
      .packages = c('matrixcalc', "MASS", "pROC", "glmnet", "mice", "filling")
    ) %dopar% {
      implement_func(
        l = l,
        m1 = m1,
        m2 = m2,
        p1 = p1
      )
    }
    
    # Combine results for this m1 value
    combined_results <- list()
    for (method in c("raw",
                     "MACOMSS",
                     "mice",
                     "mice_s",
                     "mice_cart",
                     "mice_norm",
                     "fill_KNN")) {
      combined_results[[paste0("error_", method)]] <- do.call(rbind, lapply(m1_results, function(x)
        x[[paste0("error_", method)]]))
    }
    
    results[[k]] <- c(list(m1 = m1), combined_results)
  }
  
  stopCluster(cl)
  return(results)
}
