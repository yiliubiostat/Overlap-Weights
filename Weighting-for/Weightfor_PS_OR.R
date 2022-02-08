### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

####### Data generating process (DGP): treatment and outcome

### by Yi Liu
### Sept 11, 2021

library(PSweight)
library(ggplot2)
library(dplyr)

load("covar_sims.RData")

M = max(covar$Ite)
M0 = round(M/1.5)
N = nrow(covar[covar$Ite==1,])

# Generating propensity score model
Z.reps <- function(alpha0, alpha1) {
  
  # alpha0: intercept
  # alpha1: slope vector of covariates
  # covar: covariates by above programs
  
  # The purpose of using this function is to customize different alpha's
  
  # Simulate Z using all X1~X7, so this will be the correctly specified model later
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

# Outcome regression models
PS.model.reps <- function(alpha0, alpha1) {
  
  PS.data <- Z.reps(alpha0 = alpha0, alpha1 = alpha1)
  
  # Y is the outcome vector
  Y.c <- rep(NA, M*N)
  Y.h <- rep(NA, M*N)
  
  for(i in 1:M) {
    
    # --- Constant treatment effect, WATE = 4
    epsilon <- rnorm(n = N, mean = 0, sd = 2)
    Y.c[(1+N*(i-1)):(N*i)] <- 4*PS.data[(1+N*(i-1)):(N*i), "Z"] + 0.5 + PS.data[(1+N*(i-1)):(N*i), "X1"] + 
      0.6*PS.data[(1+N*(i-1)):(N*i), "X2"] + 2.2*PS.data[(1+N*(i-1)):(N*i), "X3"] - 
      1.2*PS.data[(1+N*(i-1)):(N*i), "X4"] + 
      3*(PS.data[(1+N*(i-1)):(N*i), "X1"] + PS.data[(1+N*(i-1)):(N*i), "X2"])^2 + 
      epsilon
    
    # --- Heterogeneous treatment effect
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

