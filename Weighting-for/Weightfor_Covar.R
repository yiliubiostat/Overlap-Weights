### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

####### Data generating process (DGP): covariates

### by Yi Liu
### Sept 11, 2021

# Sample size
N <- 1000

# --- No. of replicates: 2000
M0 <- 2000
M <- round(1.5*M0)

# Generate covariates
library(MASS)
set.seed(13456)

covar <- data.frame(
  Ite = rep(1:M, rep(N, M)),
  X1 = rep(0, M*N),
  X2 = rep(0, M*N),
  X3 = rep(0, M*N),
  X4 = rep(0, M*N),
  X5 = rep(0, M*N),
  X6 = rep(0, M*N),
  X7 = rep(0, M*N)  )

for(i in 1:M) {
  
  # X3, X4: binary
  covar$X4[(1+N*(i-1)):(N*i)] = rbinom(n = N, size = 1, prob = 0.5)
  X4 <- covar$X4[(1+N*(i-1)):(N*i)]
  covar$X3[(1+N*(i-1)):(N*i)] = rbinom(n = N, size = 1, prob = 0.6*X4+0.4*(1-X4))
  X3 <- covar$X3[(1+N*(i-1)):(N*i)]
  
  mu <- cbind(-X3 + X4 + 0.5*X3*X4, X3 - X4 + X3*X4)
  
  A1 <- matrix(c(1, .5, .5, 1), ncol = 2)
  A2 <- matrix(c(2, .25, .25, 2), ncol = 2)
  
  BN <- data.frame(X1 = rep(NA, N), X2 = rep(NA, N))
  
  for(k in 1:N) {
    BN[k,] = mvrnorm(n = 1, mu = mu[k,], Sigma = X3[k]*A1 + (1-X3[k])*A2)
  }
  
  # (X1, X2) joint Gaussian
  covar$X1[(1+N*(i-1)):(N*i)] <- BN$X1
  covar$X2[(1+N*(i-1)):(N*i)] <- BN$X2
  
  # Interactions X5~X7
  covar$X5[(1+N*(i-1)):(N*i)] <- (BN$X1)^2
  covar$X6[(1+N*(i-1)):(N*i)] <- (BN$X1)*(BN$X2)
  covar$X7[(1+N*(i-1)):(N*i)] <- (BN$X2)^2
  
  rm(BN, A1, A2)
}

# save(file = "covar_sims.RData", covar)
