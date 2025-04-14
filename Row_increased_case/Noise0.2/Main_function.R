# Packages
library(matrixcalc)
library(MASS)
library(pROC)
library(glmnet)
library(matrixcalc)
library(MASS)
library(pROC)
library(glmnet)
library(mice)
library(filling)

# Functions
MACOMSS.parsvd <- function(Y11, Y12, Y21){
  m1 <- nrow(Y11) 
  m2 <- ncol(Y11) 
  p1 <- m1 + nrow(Y21) 
  p2 <- m2 + ncol(Y12) 
  M11 <- !is.na(Y11) 
  M12 <- !is.na(Y12) 
  M21 <- !is.na(Y21) 
  M_1dot <- cbind(M11, M12) 
  M_dot1 <- rbind(M11, M21) 
  hat_theta_1 <- pmax(rowSums(M_dot1), 1) %o% pmax(colSums(M_dot1), 1) / pmax(sum(M_dot1), 1) 
  hat_theta_2 <- pmax(rowSums(M_1dot), 1) %o% pmax(colSums(M_1dot), 1) / pmax(sum(M_1dot), 1) 
  hat_theta11 <- 2/(1/hat_theta_1[1:m1,1:m2] + 1/hat_theta_2[1:m1,1:m2]) 
  hat_theta12 <- hat_theta_2[,(m2+1):p2] 
  hat_theta21 <- hat_theta_1[(m1+1):p1,] 
  tilde_Y11 <- Y11
  tilde_Y12 <- Y12
  tilde_Y21 <- Y21
  tilde_Y11[is.na(tilde_Y11)] <- 0 
  tilde_Y12[is.na(tilde_Y12)] <- 0 
  tilde_Y21[is.na(tilde_Y21)] <- 0 
  tilde_Y11 <- tilde_Y11 / hat_theta11 
  tilde_Y21 <- tilde_Y21 / hat_theta21 
  tilde_Y12 <- tilde_Y12 / hat_theta12 
  
  
  list1 <- svd(rbind(tilde_Y11, tilde_Y21))
  list2 <- svd(cbind(tilde_Y11, tilde_Y12)) 
  B11 <- t(list2$u) %*% tilde_Y11 %*% list1$v 
  B12 <- t(list2$u) %*% tilde_Y12
  B21 <- tilde_Y21 %*% list1$v 
  m <- min(m1, m2)
  hatr = 0
  for (s in m:1){
    if(!is.singular.matrix(as.matrix(B11[1:s,1:s]))){
      if (norm(B21[, 1:s] %*% ginv(B11[1:s, 1:s]), "2") <= 2*sqrt(p1/m1) | norm(ginv(B11[1:s, 1:s]) %*% B12[1:s, ], "2") <= 2*sqrt(p2/m2)){
        hatr <- s; break
      }
    }
  }
  if (hatr >0 ){
    B11r <- matrix(B11[1:hatr, 1:hatr], hatr, hatr)
    B12r <- matrix(B12[1:hatr, ], hatr, p2-m2)
    B21r <- matrix(B21[,1:hatr], p1-m1, hatr)
    hatA11 <- list2$u[,1:hatr] %*% B11r %*% t(list1$v[,1:hatr])
    hatA12 <- list2$u[,1:hatr] %*% B12r
    hatA21 <- B21r %*% t(list1$v[,1:hatr])
    hatA22 <- B21r %*% ginv(B11r) %*% B12r
    hatA <- rbind(cbind(hatA11, hatA12), cbind(hatA21, hatA22))
  } else {
    r_new <- nrow(Y21)
    c_new <- ncol(Y12)
    Y22 <- matrix(0, r_new,c_new)
    Y11[is.na(Y11)] <- 0
    Y12[is.na(Y12)] <- 0
    Y21[is.na(Y21)] <- 0
    Y_1 <- rbind(Y11,Y21)
    Y_2 <- rbind(Y12,Y22)
    hatA <- cbind(Y_1,Y_2)
    # hatA22 <- matrix(0, p1 - m1, p2 - m2)
    # hatA <- matrix(0, p1, p2)
  }
  return(list(hatA, hatr))
}

MACOMSS <- function(Y11, Y12, Y21){
  m1 <- nrow(Y11)
  m2 <- ncol(Y11)
  p1 <- m1 + nrow(Y21)
  p2 <- m2 + ncol(Y12)
  M11 <- !is.na(Y11)
  M12 <- !is.na(Y12)
  M21 <- !is.na(Y21)
  M_1dot <- cbind(M11, M12)
  M_dot1 <- rbind(M11, M21)
  hat_theta_1 <- rowSums(M_dot1) %o% colSums(M_dot1) / sum(M_dot1)
  hat_theta_2 <- rowSums(M_1dot) %o% colSums(M_1dot) / sum(M_1dot)
  hat_theta11 <- (hat_theta_1[1:m1,1:m2] + hat_theta_2[1:m1,1:m2])/2
  hat_theta12 <- hat_theta_2[,(m2+1):p2]
  hat_theta21 <- hat_theta_1[(m1+1):p1,]
  tilde_Y11 <- Y11
  tilde_Y12 <- Y12
  tilde_Y21 <- Y21
  tilde_Y11[is.na(tilde_Y11)] <- 0
  tilde_Y12[is.na(tilde_Y12)] <- 0
  tilde_Y21[is.na(tilde_Y21)] <- 0 
  tilde_Y11 <- tilde_Y11 / hat_theta11
  tilde_Y21 <- tilde_Y21 / hat_theta21
  tilde_Y12 <- tilde_Y12 / hat_theta12
  
  # Apply recursive algorithm to recover A. Hard thresholding
  list1 <- svd(rbind(tilde_Y11, tilde_Y21))
  list2 <- svd(cbind(tilde_Y11, tilde_Y12))
  B11 <- t(list2$u) %*% tilde_Y11 %*% list1$v
  B12 <- t(list2$u) %*% tilde_Y12
  B21 <- tilde_Y21 %*% list1$v
  m <- min(m1, m2)
  hatr = 0
  for (s in m:1){
    if(!is.singular.matrix(as.matrix(B11[1:s,1:s]))){
      if (norm(B21[, 1:s] %*% ginv(B11[1:s, 1:s]), "2") <= 2*sqrt(p1/m1) | norm(ginv(B11[1:s, 1:s]) %*% B21[1:s, ], "2") <= 2*sqrt(p2/m2)){
        hatr <- s; break
      }
    }
  }
  if (hatr >0 ){
    B11r <- matrix(B11[1:hatr, 1:hatr], hatr, hatr)
    B12r <- matrix(B12[1:hatr, ], hatr, p2-m2)
    B21r <- matrix(B21[,1:hatr], p1-m1, hatr)
    hatA11 <- list2$u[,1:hatr] %*% B11r %*% t(list1$v[,1:hatr])
    hatA12 <- list2$u[,1:hatr] %*% B12r
    hatA21 <- B21r %*% t(list1$v[,1:hatr])
    hatA22 <- B21r %*% ginv(B11r) %*% B12r
    hatA <- rbind(cbind(hatA11, hatA12), cbind(hatA21, hatA22))
  } else {
    hatA22 <- matrix(0, p1 - m1, p2 - m2)
    hatA <- matrix(0, p1, p2)
  }
  return(list(hatA, hatr))
}

sample_simulate <- function(nr, nc) { 
  r=3
  U = svd(matrix(rnorm(nr*r), nrow=nr, ncol=r))$u
  V = svd(matrix(rnorm(nc*r), nrow=nr, ncol=r))$u
  D = diag(r)
  A = U %*% D %*% t(V) 
  return(A)
}


sample_simulate_1 <- function(nr, nc) {
  # Generate predictors
  r =5
  U = svd(matrix(rnorm(nr*r,sd=5), nrow=nr, ncol=r))$u
  V = svd(matrix(rnorm((nc)*r), nrow=nc, ncol=r))$u  # Changed nc to nc-1
  D = diag(rgamma(5,2,2))
  A = U %*% D %*% t(V)
  A = scale(A)
  sig = norm(A, 'F')/sqrt(nrow(A)*ncol(A)) * 0.2 # 0.1,0.2,0.3
  Y = A + matrix(sig*rnorm(nrow(A)*ncol(A)), nrow=nrow(A), ncol=ncol(A)) 
  X <- as.matrix(A)  # X is now nr * (nc-1)
  
  ones_column <- rep(1, nrow(A))
  bipro <- cbind(ones_column, A)
  design_mat <- as.matrix(bipro)
  
  beta_norandom <- c(rnorm(15), rnorm(nc-15,0,0))
  beta <- sample(beta_norandom)
  random_intercp <- rnorm(1, mean = 0, sd = 1)
  beta_all <- c(random_intercp,beta)

  # Generate outcome
  logit <- X %*% beta + random_intercp  # X and beta dimensions now match
  prob <- 1 / (1 + exp(-logit))
  y <- rbinom(nr, 1, prob)
  
  label_nonbernouli <- ifelse(prob < 0.5, 0, 1)
  # Combine predictors and outcome
  result <- cbind(X, y)
  
  return(list(result, beta_all, y, Y, label_nonbernouli,logit,design_mat,X))
}


missing_s <- function(x){
  theta1 <- 1-runif(nrow(x))*0.1
  theta2 <- 1-runif(ncol(x))*0.1
  theta <- theta1%o%theta2
  M <- matrix(rbinom(nrow(x)*ncol(x), 1, vec(theta)), nrow(x), ncol(x))
  x[M==0] <- NA
  return(x)
}

missing_simulate <- function(x,missing_rate){
  indic <- runif((nrow(x)*ncol(x)), min=0, max=1)
  result_missing <- ifelse(indic < missing_rate, NA, x)
  result_missing <- matrix(result_missing, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  return(result_missing)
}


label_binary <- function(x) {
  nn <-ncol(x)
  random_intercp <- rnorm(1, mean = 0, sd = 1 )
  random_number <- rnorm(nn, mean = 0.5, sd = 1)
  # random_intercp <- runif(1, min = -10, max = 10)
  # random_number <- runif(nn, min = -10, max = 10)
  random_number <- as.matrix(random_number)
  prob <- x%*%random_number + random_intercp
  label_c <- 1/(1+exp(-prob))
  
  class_labels <- matrix(0,length(label_c))
  for (i in 1:length(label_c)){
    class_labels[i] <- rbinom(1, 1, prob = label_c[i])
  }
  
  class_x <- data.frame(x, class_labels)
  class_x <- as.matrix(class_x)
  random_beta <- c(random_intercp, random_number)
  return(list(class_x, random_beta,label_c))
}

sigmoid_binary_y <- function(x) {
  y <- 1 / (1 + exp(-x))
  A <- ifelse(y >= 0.5, 1, 0)
  return(A)
}
