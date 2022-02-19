### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Variance estimations ATE ATT ATC ~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### Data generating process (DGP): PS and OR models

### by Yi Liu
### Oct 24, 2021

load("covar_sims.RData")

M = max(covar$Ite)
M0 = round(M/1.5)
N = nrow(covar[covar$Ite==1,])

# A function for generating PS model (Z)
Z.reps <- function(alpha0, alpha1) {
  
  # prop: proportion of people who receive the treatment
  # alpha0: intercept of the linear model
  # alpha1: slope vector of X
  # Covar: covariates data, should be a matrix or data frame
  
  # Simulate correct Z
  # --- M and N are from "Weightfor_Covar_Sims.R" program
  Z <- rep(NA, M*N)
  
  for(i in 1:M) {
    
    logit <- alpha1 %*% t(as.matrix(covar[(1+N*(i-1)):(N*i),-1])) + rep(alpha0, N)
    e <- as.vector(exp(logit)/(1+exp(logit)))
    
    Z[(1+N*(i-1)):(N*i)] <- rbinom(n = N, size = 1, prob = e)
  }
  
  ZZ <- ifelse(Z==1, 0, 1)
  PS.data <- cbind(covar, Z, ZZ)
  
  return(PS.data)
  rm(logit, e)
}

# Generate outcome models
PS.model.reps <- function(alpha0, alpha1) {
  
  # First, use the PS.model.reps function we write to generate the data containing Z
  # --- Note that this function generates both correct and misspecified models
  PS.data <- Z.reps(alpha0 = alpha0, alpha1 = alpha1)
  
  # Generate Y, similarly to true value program
  Y.c <- rep(NA, M*N)
  Y.h <- rep(NA, M*N)
 
  for(i in 1:M) {
    
    # --- Y is correct and constant
    epsilon <- rnorm(n = N, mean = 0, sd = 2)
    Y.c[(1+N*(i-1)):(N*i)] <- 4*PS.data[(1+N*(i-1)):(N*i), "Z"] + 0.5 + PS.data[(1+N*(i-1)):(N*i), "X1"] + 
      0.6*PS.data[(1+N*(i-1)):(N*i), "X2"] + 2.2*PS.data[(1+N*(i-1)):(N*i), "X3"] - 
      1.2*PS.data[(1+N*(i-1)):(N*i), "X4"] + 
      3*(PS.data[(1+N*(i-1)):(N*i), "X1"] + PS.data[(1+N*(i-1)):(N*i), "X2"])^2 + 
      epsilon
    
    # --- Y is correct and heterogeneous
    epsilon <- rnorm(n = N, mean = 0, sd = 2)
    Y.h[(1+N*(i-1)):(N*i)] <- 0.5 + PS.data[(1+N*(i-1)):(N*i), "X1"] + 0.6*PS.data[(1+N*(i-1)):(N*i), "X2"] + 
      2.2*PS.data[(1+N*(i-1)):(N*i), "X3"] - 1.2*PS.data[(1+N*(i-1)):(N*i), "X4"] + 
      (PS.data[(1+N*(i-1)):(N*i), "X1"] + PS.data[(1+N*(i-1)):(N*i), "X2"])^2 + 
      (4 + 3*(PS.data[(1+N*(i-1)):(N*i), "X1"] + PS.data[(1+N*(i-1)):(N*i), "X2"])^2 + PS.data[(1+N*(i-1)):(N*i), "X1"]*PS.data[(1+N*(i-1)):(N*i), "X3"])*PS.data[(1+N*(i-1)):(N*i), "Z"] + 
      epsilon
  }
  PSweight.data.reps <- cbind(PS.data, Y.c, Y.h)
  return(PSweight.data.reps)
}
