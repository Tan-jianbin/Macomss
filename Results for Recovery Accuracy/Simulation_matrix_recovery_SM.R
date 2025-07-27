# Noisy SMC Simulation

## Structured matrix completion with noisy and missing entries

library(matrixcalc)
library(MASS)

################################### Setting 1 ##################################

# Parameters that are allowed to varies, p_1, p_2, m_1, m_2, theta. 
# In the first setting, we consider the effect of m_1, m_2

p1 = 300; p2 = 300; m = (1:10)*10

lmax = 1000;
loss_F = matrix(0, nrow = length(m), ncol = length(m)); 
loss_spe = matrix(0, nrow = length(m), ncol = length(m)); 

ptm <- proc.time() # Timing start!
for (l in 1:lmax){
  r=3
  U = svd(matrix(rnorm(p1*r), nrow=p1, ncol=r))$u
  V = svd(matrix(rnorm(p2*r), nrow=p1, ncol=r))$u
  D = diag(r)
  A = U %*% D %*% t(V) # Original matrix
  sig = norm(A, 'F')/sqrt(p1*p2) * 0.3
  
  
  for (i in 1:length(m)){
    for (j in 1:length(m)){
      m1 = m[i]; m2 = m[j]; r = 3;
      
      
      Y = A + matrix(sig*rnorm(p1*p2), nrow=p1, ncol=p2) # Observations without missing values
      theta1 <- 1-runif(p1)*0.05
      theta2 <- 1-runif(p2)*0.05
      theta <- theta1%o%theta2
      M <- matrix(rbinom(p1*p2, 1, vec(theta)), p1, p2)
      Y[M==0] <- NA
      Y11 = Y[1:m1, 1:m2]; Y12 = Y[1:m1, (m2+1):p2]; Y21 = Y[(m1+1):p1, 1:m2]
      
      result = smc.parsvd(Y11, Y12, Y21)
      hatA = result[[1]]
      loss_F[i,j] <- loss_F[i,j] + norm(hatA - A, 'F')/lmax
      loss_spe[i,j] <- loss_spe[i,j] + norm(hatA - A, '2')/lmax
    }
  }
  print(l)
}
proc.time() - ptm

loss_dat <- data.frame(z = as.vector(loss_F_setting1))
loss_dat$x <- rep(seq(10, 100, 10), 10)
loss_dat$y <- rep(seq(10, 100, 10), 10)
wireframe(z~x*y, loss_dat,
          xlab = expression(m[1]), ylab = expression(m[2]), zlab = "",
          main = "Average Frobenius norm loss",
          scales = list(arrows=FALSE,cex=.5,tick.number = 5),
          drape = TRUE,
          colorkey = TRUE,
          screen = list(z = -110, x = -60)
)

loss_dat <- data.frame(z = as.vector(loss_spe_setting1))
loss_dat$x <- rep(seq(10, 100, 10), 10)
loss_dat$y <- rep(seq(10, 100, 10), 10)
wireframe(z~x*y, loss_dat,
          xlab = expression(m[1]), ylab = expression(m[2]), zlab = "",
          main = "Average spectral norm loss",
          scales = list(arrows=FALSE,cex=.5,tick.number = 5),
          drape = TRUE,
          colorkey = TRUE,
          screen = list(z = -110, x = -60)
)

################################### Setting 2 ##################################

# Parameters that are allowed to varies, p_1, p_2, m_1, m_2, theta. 
# In the first setting, we consider the effect of theta and sigma

p1 = 300; p2 =300; m1 = 50; m2 = 50;
sig_choice = 0.2 + (0:8)*0.2
theta_choice = (0:8)*0.03

lmax = 1000;
loss_F = matrix(0, nrow = length(sig_choice), ncol = length(theta_choice)); 
loss_spe = matrix(0, nrow = length(sig_choice), ncol = length(theta_choice)); 

ptm <- proc.time() # Timing start
for (l in 1:lmax){
  U = svd(matrix(rnorm(p1*r), nrow=p1, ncol=r))$u
  V = svd(matrix(rnorm(p2*r), nrow=p2, ncol=r))$u
  D = diag(r)
  A = U %*% D %*% t(V) # Original matrix
  
  for (i in 1:length(sig_choice)){
    for (j in 1:length(theta_choice)){
      r = 3;
      
      sig = norm(A, 'F')/sqrt(p1*p2) * sig_choice[i]
      Y = A + matrix(sig*rnorm(p1*p2), nrow=p1, ncol=p2) # Observations without missing values
      theta1 <- 1-runif(p1)*theta_choice[j]
      theta2 <- 1-runif(p2)*theta_choice[j]
      theta <- theta1%o%theta2
      M <- matrix(rbinom(p1*p2, 1, vec(theta)), p1, p2)
      Y[M==0] <- NA
      Y11 = Y[1:m1, 1:m2]; Y12 = Y[1:m1, (m2+1):p2]; Y21 = Y[(m1+1):p1, 1:m2]
      
      #result = smc(Y11, Y12, Y21)
      #hatA = result[[1]]
      #hatr1[l] = result[[2]]
      #loss_F[l] <- norm(hatA - A, 'F')
      
      result = smc.parsvd(Y11, Y12, Y21)
      hatA = result[[1]]
      loss_F[i,j] <- loss_F[i,j] + norm(hatA - A, 'F')/lmax
      loss_spe[i,j] <- loss_spe[i,j] + norm(hatA - A, '2')/lmax
    }
  }
  print(l)
}
proc.time() - ptm
#norm(hatA - A, '2')/norm(A, '2')

# Plot algorithms
loss_dat <- data.frame(z = as.vector(loss_F_setting2))
loss_dat$x <- rep(seq(0.2, 1.8, 0.2), 9)
loss_dat$y <- rep(seq(0, 0.24, 0.03), 9)
wireframe(z~x*y, loss_dat,
          xlab = expression(sigma), ylab = expression(eta), zlab = "",
          main = "Average Frobenius norm loss",
          scales = list(arrows=FALSE,cex=.5,tick.number = 5),
          drape = TRUE,
          colorkey = TRUE,
          screen = list(z = 45, x = -60)
)

loss_dat <- data.frame(z = as.vector(loss_spe_setting2))
loss_dat$x <- rep(seq(0.2, 1.8, 0.2), 9)
loss_dat$y <- rep(seq(0, 0.24, 0.03), 9)
wireframe(z~x*y, loss_dat,
          xlab = expression(sigma), ylab = expression(eta), zlab = "",
          main = "Average spectral norm loss",
          scales = list(arrows=FALSE,cex=.5,tick.number = 5),
          drape = TRUE,
          colorkey = TRUE,
          screen = list(z = 45, x = -60)
)


################################### Setting 3 ##################################

# In this setting, we consider the more general setting when A is approximately low-rank instead of exactly low-rank
# To be specific, let 
library(pracma)
p1 = 300; p2 =300; m = c(10, 20, 50, 100); alpha = seq(.6, 2.4, 0.2); r = 3;

lmax = 1000;
loss_F = matrix(0, nrow = length(m), ncol = length(alpha)); 
loss_spe = matrix(0, nrow = length(m), ncol = length(alpha)); 

ptm <- proc.time() # Timing start!
for (l in 1:lmax){
  U = randortho(p1)
  V = randortho(p2)
  for (i in 1:length(m)){
    for (j in 1:length(alpha)){
      m1 = m[i]; m2=m[i]
      D = diag(c(rep(1, r), (1:(min(p1, p2) - r))^(-alpha[j])))
      A = U %*% D %*% t(V) # Original matrix
      sig = norm(A, 'F')/sqrt(p1*p2) * 0.3
      Y = A + matrix(sig*rnorm(p1*p2), nrow=p1, ncol=p2) # Observations without missing values
      theta1 <- 1-runif(p1)*0.05
      theta2 <- 1-runif(p2)*0.05
      theta <- theta1%o%theta2
      M <- matrix(rbinom(p1*p2, 1, vec(theta)), p1, p2)
      Y[M==0] <- NA
      Y11 = Y[1:m1, 1:m2]; Y12 = Y[1:m1, (m2+1):p2]; Y21 = Y[(m1+1):p1, 1:m2]
      
      #result = smc(Y11, Y12, Y21)
      #hatA = result[[1]]
      #hatr1[l] = result[[2]]
      #loss_F[l] <- norm(hatA - A, 'F')
      
      result = smc.parsvd(Y11, Y12, Y21)
      hatA = result[[1]]
      loss_F[i,j] <- loss_F[i,j] + norm(hatA - A, 'F')/lmax
      loss_spe[i,j] <- loss_spe[i,j] + norm(hatA - A, '2')/lmax
    }
  }
  print(l)
}
proc.time() - ptm

plot(seq(0.6, 2.4, 0.2), loss_F_setting3[1,], type="b", col="red", pch=0, lty=1, ylim=c(0.2, 3), xlab=expression(alpha), ylab="Average Frobenius norm loss")
lines(seq(0.6, 2.4, 0.2), loss_F_setting3[2,], type = "b", pch=1, lty=2, col="green")
lines(seq(0.6, 2.4, 0.2), loss_F_setting3[3,], type = "b", pch=2, lty=3, col="blue")
lines(seq(0.6, 2.4, 0.2), loss_F_setting3[4,], type = "b", pch=3, lty=4, col="purple")
legend("topright",legend=c("m=10", "m=20", "m=50", "m=100"), col=c("red", "green", "blue", "purple"), lty=c(1,2,3,4), pch=c(1,2,3,4))

plot(seq(0.6, 2.4, 0.2), loss_spe_setting3[1,], type="b", col="red", pch=0, lty=1, ylim=c(0.1, 2), xlab=expression(alpha), ylab="Average spectral norm loss")
lines(seq(0.6, 2.4, 0.2), loss_spe_setting3[2,], type = "b", pch=1, lty=2, col="green")
lines(seq(0.6, 2.4, 0.2), loss_spe_setting3[3,], type = "b", pch=2, lty=3, col="blue")
lines(seq(0.6, 2.4, 0.2), loss_spe_setting3[4,], type = "b", pch=3, lty=4, col="purple")
legend("topright",legend=c("m=10", "m=20", "m=50", "m=100"), col=c("red", "green", "blue", "purple"), lty=c(1,2,3,4), pch=c(1,2,3,4))

################################### Setting 4 ##################################

# In this Setting 4, we consider the special setting of Poisson noise
# Randomly construct A as non-negative low-rank matrix by A = U * t(V).

p1 = 300; p2 =300; m = c(10, 20, 50, 100); lambda0 = (1:10)

lmax = 1000;
loss_F = matrix(0, nrow = length(m), ncol = length(lambda0)); 
loss_spe = matrix(0, nrow = length(m), ncol = length(lambda0)); 

ptm <- proc.time() # Timing start!
for (l in 1:lmax){
  U = abs(matrix(rnorm(p1*r), nrow=p1, ncol=r))/sqrt(p1)
  V = abs(matrix(rnorm(p2*r), nrow=p2, ncol=r))/sqrt(p2)
  D = diag(r)
  A = U %*% D %*% t(V) # Original matrix
  
  for (i in 1:length(m)){
    for (j in 1:length(lambda0)){
      m1 = m[i]; m2 = m[i]; r = 3;
      lambda = p1*p2/sum(A) * lambda0[j]
      Y = matrix(rpois(p1*p2, lambda*A), nrow=p1, ncol=p2) # Observations without missing values
      theta1 <- 1-runif(p1)*0.05
      theta2 <- 1-runif(p2)*0.05
      theta <- theta1%o%theta2
      M <- matrix(rbinom(p1*p2, 1, vec(theta)), p1, p2)
      Y[M==0] <- NA
      Y11 = Y[1:m1, 1:m2]; Y12 = Y[1:m1, (m2+1):p2]; Y21 = Y[(m1+1):p1, 1:m2]
      
      #result = smc(Y11, Y12, Y21)
      #hatA = result[[1]]
      #hatr1[l] = result[[2]]
      #loss_F[l] <- norm(hatA - A, 'F')
      
      result = smc.parsvd(Y11, Y12, Y21)
      hatA = result[[1]]/lambda
      loss_F[i,j] <- loss_F[i,j] + norm(hatA - A, 'F')/lmax
      loss_spe[i,j] <- loss_spe[i,j] + norm(hatA - A, '2')/lmax
    }
  }
  print(l)
}
proc.time() - ptm
#norm(hatA - A, '2')/norm(A, '2')

# Plot code
#persp(seq(10, 100, 10), seq(10, 100, 10), loss_F_setting1, phi = 45, theta = 135,
#      xlab = "m_1", ylab = "m_2", zlab = "Average Frobenius norm loss",
#      main = ""
#)
plot(1:10, loss_F_setting4[1,], type="b", col="red", pch=0, lty=1, ylim=c(0.2, 1.5), xlab=expression(lambda[0]), ylab="Average Frobenius norm loss")
lines(1:10, loss_F_setting4[2,], type = "b", pch=1, lty=2, col="green")
lines(1:10, loss_F_setting4[3,], type = "b", pch=2, lty=3, col="blue")
lines(1:10, loss_F_setting4[4,], type = "b", pch=3, lty=4, col="purple")
legend("topright",legend=c("m=10", "m=20", "m=50", "m=100"), col=c("red", "green", "blue", "purple"), lty=c(1,2,3,4), pch=c(1,2,3,4))

plot(1:10, loss_spe_setting4[1,], type="b", col="red", pch=0, lty=1, ylim=c(0.1, 1.2), xlab=expression(lambda[0]), ylab="Average spectral norm loss")
lines(1:10, loss_spe_setting4[2,], type = "b", pch=1, lty=2, col="green")
lines(1:10, loss_spe_setting4[3,], type = "b", pch=2, lty=3, col="blue")
lines(1:10, loss_spe_setting4[4,], type = "b", pch=3, lty=4, col="purple")
legend("topright",legend=c("m=10", "m=20", "m=50", "m=100"), col=c("red", "green", "blue", "purple"), lty=c(1,2,3,4), pch=c(1,2,3,4))

############################ Ploting code ###############

loss_dat <- data.frame(z = as.vector(loss_F_setting1))
loss_dat$x <- rep(seq(10, 100, 10), 10)
loss_dat$y <- rep(seq(10, 100, 10), 10)
wireframe(z~x*y, loss_dat,
          xlab = expression(m[1]), ylab = expression(m[2]), zlab = "",
          main = "Average Frobenius norm loss",
          scales = list(arrows=FALSE,cex=.5,tick.number = 5),
          drape = TRUE,
          colorkey = TRUE,
          screen = list(z = -110, x = -60)
)

loss_dat <- data.frame(z = as.vector(loss_spe_setting1))
loss_dat$x <- rep(seq(10, 100, 10), 10)
loss_dat$y <- rep(seq(10, 100, 10), 10)
wireframe(z~x*y, loss_dat,
          xlab = expression(m[1]), ylab = expression(m[2]), zlab = "",
          main = "Average spectral norm loss",
          scales = list(arrows=FALSE,cex=.5,tick.number = 5),
          drape = TRUE,
          colorkey = TRUE,
          screen = list(z = -110, x = -60)
)


################################ Recovery procedure ###########################
smc.parsvd <- function(Y11, Y12, Y21){
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
  
  # Apply recursive algorithm to recover A. Hard thresholding
  #list1 <- svd(rbind(tilde_Y11, tilde_Y21))
  list1 <- svd(rbind(tilde_Y21))
  #list2 <- svd(cbind(tilde_Y11, tilde_Y12))
  list2 <- svd(cbind(tilde_Y12))
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
    hatA22 <- matrix(0, p1 - m1, p2 - m2)
    hatA <- matrix(0, p1, p2)
  }
  return(list(hatA, hatr))
}

smc <- function(Y11, Y12, Y21){
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

################################### Setting 3 ##################################

# In this setting, we consider the more general setting when A is approximately low-rank instead of exactly low-rank
# To be specific, let 
library(pracma)
p1 = 300; p2 =300; m = c(10, 20, 50, 100); alpha = seq(.6, 2.4, 0.2); r = 3;

lmax = 1000;
loss_F = matrix(0, nrow = length(m), ncol = length(alpha)); 
loss_spe = matrix(0, nrow = length(m), ncol = length(alpha)); 

ptm <- proc.time() # Timing start!
for (l in 1:lmax){
  U = randortho(p1)
  V = randortho(p2)
  for (i in 1:length(m)){
    for (j in 1:length(alpha)){
      m1 = m[i]; m2=m[i]
      D = diag(c(rep(1, r), (1:(min(p1, p2) - r))^(-alpha[j])))
      A = U %*% D %*% t(V) # Original matrix
      sig = norm(A, 'F')/sqrt(p1*p2) * 0.3
      Y = A + matrix(sig*rnorm(p1*p2), nrow=p1, ncol=p2) # Observations without missing values
      theta1 <- 1-runif(p1)*0.05
      theta2 <- 1-runif(p2)*0.05
      theta <- theta1%o%theta2
      M <- matrix(rbinom(p1*p2, 1, vec(theta)), p1, p2)
      Y[M==0] <- NA
      Y11 = Y[1:m1, 1:m2]; Y12 = Y[1:m1, (m2+1):p2]; Y21 = Y[(m1+1):p1, 1:m2]
      
      #result = smc(Y11, Y12, Y21)
      #hatA = result[[1]]
      #hatr1[l] = result[[2]]
      #loss_F[l] <- norm(hatA - A, 'F')
      
      result = smc.parsvd(Y11, Y12, Y21)
      hatA = result[[1]]
      loss_F[i,j] <- loss_F[i,j] + norm(hatA - A, 'F')/lmax
      loss_spe[i,j] <- loss_spe[i,j] + norm(hatA - A, '2')/lmax
    }
  }
  print(l)
}
proc.time() - ptm

plot(seq(0.6, 2.4, 0.2), loss_F_setting3[1,], type="b", col="red", pch=0, lty=1, ylim=c(0.2, 3), xlab=expression(alpha), ylab="Average Frobenius norm loss")
lines(seq(0.6, 2.4, 0.2), loss_F_setting3[2,], type = "b", pch=1, lty=2, col="green")
lines(seq(0.6, 2.4, 0.2), loss_F_setting3[3,], type = "b", pch=2, lty=3, col="blue")
lines(seq(0.6, 2.4, 0.2), loss_F_setting3[4,], type = "b", pch=3, lty=4, col="purple")
legend("topright",legend=c("m=10", "m=20", "m=50", "m=100"), col=c("red", "green", "blue", "purple"), lty=c(1,2,3,4), pch=c(1,2,3,4))

plot(seq(0.6, 2.4, 0.2), loss_spe_setting3[1,], type="b", col="red", pch=0, lty=1, ylim=c(0.1, 2), xlab=expression(alpha), ylab="Average spectral norm loss")
lines(seq(0.6, 2.4, 0.2), loss_spe_setting3[2,], type = "b", pch=1, lty=2, col="green")
lines(seq(0.6, 2.4, 0.2), loss_spe_setting3[3,], type = "b", pch=2, lty=3, col="blue")
lines(seq(0.6, 2.4, 0.2), loss_spe_setting3[4,], type = "b", pch=3, lty=4, col="purple")
legend("topright",legend=c("m=10", "m=20", "m=50", "m=100"), col=c("red", "green", "blue", "purple"), lty=c(1,2,3,4), pch=c(1,2,3,4))

################################### Setting 4 ##################################

# In this Setting 4, we consider the special setting of Poisson noise
# Randomly construct A as non-negative low-rank matrix by A = U * t(V).

p1 = 300; p2 =300; m = c(10, 20, 50, 100); lambda0 = (1:10)

lmax = 1000;
loss_F = matrix(0, nrow = length(m), ncol = length(lambda0)); 
loss_spe = matrix(0, nrow = length(m), ncol = length(lambda0)); 

ptm <- proc.time() # Timing start!
for (l in 1:lmax){
  U = abs(matrix(rnorm(p1*r), nrow=p1, ncol=r))/sqrt(p1)
  V = abs(matrix(rnorm(p2*r), nrow=p2, ncol=r))/sqrt(p2)
  D = diag(r)
  A = U %*% D %*% t(V) # Original matrix
  
  for (i in 1:length(m)){
    for (j in 1:length(lambda0)){
      m1 = m[i]; m2 = m[i]; r = 3;
      lambda = p1*p2/sum(A) * lambda0[j]
      Y = matrix(rpois(p1*p2, lambda*A), nrow=p1, ncol=p2) # Observations without missing values
      theta1 <- 1-runif(p1)*0.05
      theta2 <- 1-runif(p2)*0.05
      theta <- theta1%o%theta2
      M <- matrix(rbinom(p1*p2, 1, vec(theta)), p1, p2)
      Y[M==0] <- NA
      Y11 = Y[1:m1, 1:m2]; Y12 = Y[1:m1, (m2+1):p2]; Y21 = Y[(m1+1):p1, 1:m2]
      
      #result = smc(Y11, Y12, Y21)
      #hatA = result[[1]]
      #hatr1[l] = result[[2]]
      #loss_F[l] <- norm(hatA - A, 'F')
      
      result = smc.parsvd(Y11, Y12, Y21)
      hatA = result[[1]]/lambda
      loss_F[i,j] <- loss_F[i,j] + norm(hatA - A, 'F')/lmax
      loss_spe[i,j] <- loss_spe[i,j] + norm(hatA - A, '2')/lmax
    }
  }
  print(l)
}
proc.time() - ptm
#norm(hatA - A, '2')/norm(A, '2')

# Plot code
#persp(seq(10, 100, 10), seq(10, 100, 10), loss_F_setting1, phi = 45, theta = 135,
#      xlab = "m_1", ylab = "m_2", zlab = "Average Frobenius norm loss",
#      main = ""
#)
plot(1:10, loss_F_setting4[1,], type="b", col="red", pch=0, lty=1, ylim=c(0.2, 1.5), xlab=expression(lambda[0]), ylab="Average Frobenius norm loss")
lines(1:10, loss_F_setting4[2,], type = "b", pch=1, lty=2, col="green")
lines(1:10, loss_F_setting4[3,], type = "b", pch=2, lty=3, col="blue")
lines(1:10, loss_F_setting4[4,], type = "b", pch=3, lty=4, col="purple")
legend("topright",legend=c("m=10", "m=20", "m=50", "m=100"), col=c("red", "green", "blue", "purple"), lty=c(1,2,3,4), pch=c(1,2,3,4))

plot(1:10, loss_spe_setting4[1,], type="b", col="red", pch=0, lty=1, ylim=c(0.1, 1.2), xlab=expression(lambda[0]), ylab="Average spectral norm loss")
lines(1:10, loss_spe_setting4[2,], type = "b", pch=1, lty=2, col="green")
lines(1:10, loss_spe_setting4[3,], type = "b", pch=2, lty=3, col="blue")
lines(1:10, loss_spe_setting4[4,], type = "b", pch=3, lty=4, col="purple")
legend("topright",legend=c("m=10", "m=20", "m=50", "m=100"), col=c("red", "green", "blue", "purple"), lty=c(1,2,3,4), pch=c(1,2,3,4))

############################ Ploting code ###############

loss_dat <- data.frame(z = as.vector(loss_F_setting1))
loss_dat$x <- rep(seq(10, 100, 10), 10)
loss_dat$y <- rep(seq(10, 100, 10), 10)
wireframe(z~x*y, loss_dat,
          xlab = expression(m[1]), ylab = expression(m[2]), zlab = "",
          main = "Average Frobenius norm loss",
          scales = list(arrows=FALSE,cex=.5,tick.number = 5),
          drape = TRUE,
          colorkey = TRUE,
          screen = list(z = -110, x = -60)
)

loss_dat <- data.frame(z = as.vector(loss_spe_setting1))
loss_dat$x <- rep(seq(10, 100, 10), 10)
loss_dat$y <- rep(seq(10, 100, 10), 10)
wireframe(z~x*y, loss_dat,
          xlab = expression(m[1]), ylab = expression(m[2]), zlab = "",
          main = "Average spectral norm loss",
          scales = list(arrows=FALSE,cex=.5,tick.number = 5),
          drape = TRUE,
          colorkey = TRUE,
          screen = list(z = -110, x = -60)
)


################################ Recovery procedure ###########################
smc.parsvd <- function(Y11, Y12, Y21){
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
  
  # Apply recursive algorithm to recover A. Hard thresholding
  #list1 <- svd(rbind(tilde_Y11, tilde_Y21))
  list1 <- svd(rbind(tilde_Y21))
  #list2 <- svd(cbind(tilde_Y11, tilde_Y12))
  list2 <- svd(cbind(tilde_Y12))
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
    hatA22 <- matrix(0, p1 - m1, p2 - m2)
    hatA <- matrix(0, p1, p2)
  }
  return(list(hatA, hatr))
}

smc <- function(Y11, Y12, Y21){
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

print(l)
proc.time() - ptm

plot(seq(0.6, 2.4, 0.2), loss_F_setting3[1,], type="b", col="red", pch=0, lty=1, ylim=c(0.2, 3), xlab=expression(alpha), ylab="Average Frobenius norm loss")
lines(seq(0.6, 2.4, 0.2), loss_F_setting3[2,], type = "b", pch=1, lty=2, col="green")
lines(seq(0.6, 2.4, 0.2), loss_F_setting3[3,], type = "b", pch=2, lty=3, col="blue")
lines(seq(0.6, 2.4, 0.2), loss_F_setting3[4,], type = "b", pch=3, lty=4, col="purple")
legend("topright",legend=c("m=10", "m=20", "m=50", "m=100"), col=c("red", "green", "blue", "purple"), lty=c(1,2,3,4), pch=c(1,2,3,4))

plot(seq(0.6, 2.4, 0.2), loss_spe_setting3[1,], type="b", col="red", pch=0, lty=1, ylim=c(0.1, 2), xlab=expression(alpha), ylab="Average spectral norm loss")
lines(seq(0.6, 2.4, 0.2), loss_spe_setting3[2,], type = "b", pch=1, lty=2, col="green")
lines(seq(0.6, 2.4, 0.2), loss_spe_setting3[3,], type = "b", pch=2, lty=3, col="blue")
lines(seq(0.6, 2.4, 0.2), loss_spe_setting3[4,], type = "b", pch=3, lty=4, col="purple")
legend("topright",legend=c("m=10", "m=20", "m=50", "m=100"), col=c("red", "green", "blue", "purple"), lty=c(1,2,3,4), pch=c(1,2,3,4))

################################### Setting 4 ##################################

# In this Setting 4, we consider the special setting of Poisson noise
# Randomly construct A as non-negative low-rank matrix by A = U * t(V).

p1 = 300; p2 =300; m = c(10, 20, 50, 100); lambda0 = (1:10)

lmax = 1000;
loss_F = matrix(0, nrow = length(m), ncol = length(lambda0)); 
loss_spe = matrix(0, nrow = length(m), ncol = length(lambda0)); 

ptm <- proc.time() # Timing start!
for (l in 1:lmax){
  U = abs(matrix(rnorm(p1*r), nrow=p1, ncol=r))/sqrt(p1)
  V = abs(matrix(rnorm(p2*r), nrow=p2, ncol=r))/sqrt(p2)
  D = diag(r)
  A = U %*% D %*% t(V) # Original matrix
  
  for (i in 1:length(m)){
    for (j in 1:length(lambda0)){
      m1 = m[i]; m2 = m[i]; r = 3;
      lambda = p1*p2/sum(A) * lambda0[j]
      Y = matrix(rpois(p1*p2, lambda*A), nrow=p1, ncol=p2) # Observations without missing values
      theta1 <- 1-runif(p1)*0.05
      theta2 <- 1-runif(p2)*0.05
      theta <- theta1%o%theta2
      M <- matrix(rbinom(p1*p2, 1, vec(theta)), p1, p2)
      Y[M==0] <- NA
      Y11 = Y[1:m1, 1:m2]; Y12 = Y[1:m1, (m2+1):p2]; Y21 = Y[(m1+1):p1, 1:m2]
      
      #result = smc(Y11, Y12, Y21)
      #hatA = result[[1]]
      #hatr1[l] = result[[2]]
      #loss_F[l] <- norm(hatA - A, 'F')
      
      result = smc.parsvd(Y11, Y12, Y21)
      hatA = result[[1]]/lambda
      loss_F[i,j] <- loss_F[i,j] + norm(hatA - A, 'F')/lmax
      loss_spe[i,j] <- loss_spe[i,j] + norm(hatA - A, '2')/lmax
    }
  }
  print(l)
}
proc.time() - ptm
#norm(hatA - A, '2')/norm(A, '2')

# Plot code
#persp(seq(10, 100, 10), seq(10, 100, 10), loss_F_setting1, phi = 45, theta = 135,
#      xlab = "m_1", ylab = "m_2", zlab = "Average Frobenius norm loss",
#      main = ""
#)
plot(1:10, loss_F_setting4[1,], type="b", col="red", pch=0, lty=1, ylim=c(0.2, 1.5), xlab=expression(lambda[0]), ylab="Average Frobenius norm loss")
lines(1:10, loss_F_setting4[2,], type = "b", pch=1, lty=2, col="green")
lines(1:10, loss_F_setting4[3,], type = "b", pch=2, lty=3, col="blue")
lines(1:10, loss_F_setting4[4,], type = "b", pch=3, lty=4, col="purple")
legend("topright",legend=c("m=10", "m=20", "m=50", "m=100"), col=c("red", "green", "blue", "purple"), lty=c(1,2,3,4), pch=c(1,2,3,4))

plot(1:10, loss_spe_setting4[1,], type="b", col="red", pch=0, lty=1, ylim=c(0.1, 1.2), xlab=expression(lambda[0]), ylab="Average spectral norm loss")
lines(1:10, loss_spe_setting4[2,], type = "b", pch=1, lty=2, col="green")
lines(1:10, loss_spe_setting4[3,], type = "b", pch=2, lty=3, col="blue")
lines(1:10, loss_spe_setting4[4,], type = "b", pch=3, lty=4, col="purple")
legend("topright",legend=c("m=10", "m=20", "m=50", "m=100"), col=c("red", "green", "blue", "purple"), lty=c(1,2,3,4), pch=c(1,2,3,4))

############################ Ploting code ###############

loss_dat <- data.frame(z = as.vector(loss_F_setting1))
loss_dat$x <- rep(seq(10, 100, 10), 10)
loss_dat$y <- rep(seq(10, 100, 10), 10)
wireframe(z~x*y, loss_dat,
          xlab = expression(m[1]), ylab = expression(m[2]), zlab = "",
          main = "Average Frobenius norm loss",
          scales = list(arrows=FALSE,cex=.5,tick.number = 5),
          drape = TRUE,
          colorkey = TRUE,
          screen = list(z = -110, x = -60)
)

loss_dat <- data.frame(z = as.vector(loss_spe_setting1))
loss_dat$x <- rep(seq(10, 100, 10), 10)
loss_dat$y <- rep(seq(10, 100, 10), 10)
wireframe(z~x*y, loss_dat,
          xlab = expression(m[1]), ylab = expression(m[2]), zlab = "",
          main = "Average spectral norm loss",
          scales = list(arrows=FALSE,cex=.5,tick.number = 5),
          drape = TRUE,
          colorkey = TRUE,
          screen = list(z = -110, x = -60)
)


################################ Recovery procedure ###########################
smc.parsvd <- function(Y11, Y12, Y21){
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
  
  # Apply recursive algorithm to recover A. Hard thresholding
  #list1 <- svd(rbind(tilde_Y11, tilde_Y21))
  list1 <- svd(rbind(tilde_Y21))
  #list2 <- svd(cbind(tilde_Y11, tilde_Y12))
  list2 <- svd(cbind(tilde_Y12))
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
    hatA22 <- matrix(0, p1 - m1, p2 - m2)
    hatA <- matrix(0, p1, p2)
  }
  return(list(hatA, hatr))
}

smc <- function(Y11, Y12, Y21){
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