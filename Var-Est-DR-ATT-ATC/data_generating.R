### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Variance estimations ATE ATT ATC ~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### Data generating process (DGP)

### by Yi Liu
### Oct 24, 2021

# Model data for simulation
source("ps_out_md_func.R")
alpha1 <- c(0.3, 0.4, 0.4, 0.4, -0.1, -0.1, 0.1)

### Data for model 1: alpha0 = -2.17
md1.simdat <- PS.model.reps(alpha0 = -2.17, alpha1)

### Data for model 2: alpha0 = -0.78
md2.simdat <- PS.model.reps(alpha0 = -0.78, alpha1)

### Data for model 3: alpha0 = 0.98
md3.simdat <- PS.model.reps(alpha0 = 0.98, alpha1)

save(file = "model_data_sim.Rdata", md1.simdat, md2.simdat, md3.simdat)

# Truth of estimands
### Constant effect: tau = 4 for all tilting functions
### Heterogeneous effect:
# --- Choose a large N to simulate true population
N_true <- 1000000

# Generate covariates
set.seed(13244)

# --- X3 and X4: Bern(p)
X4 <- rbinom(n = N_true, size = 1, prob = 0.5)
X3 <- rbinom(n = N_true, size = 1, prob = 0.6*X4+0.4*(1-X4))

# --- X1 and X2: Bivariate Normal
library(MASS)
mu <- cbind(-X3 + X4 + 0.5*X3*X4, X3 - X4 + X3*X4)

A1 <- matrix(c(1, .5, .5, 1), ncol = 2)
A2 <- matrix(c(2, .25, .25, 2), ncol = 2)

BN <- data.frame(X1 = rep(NA, N_true), X2 = rep(NA, N_true))

for(i in 1:N_true) {
  BN[i,] = mvrnorm(n = 1, mu = mu[i,], Sigma = X3[i]*A1 + (1-X3[i])*A2)
}

# --- X5, X6 and X7
X5 <- (BN$X1)^2
X6 <- (BN$X1)*(BN$X2)
X7 <- (BN$X2)^2

covar_true <- data.frame(X1 = BN$X1, X2 = BN$X2,
                         X3 = X3, X4 = X4, X5 = X5, X6 = X6, X7 = X7)
rm(BN, A1, A2)

# load("covar_true.RData")

# A function for finding alpha0
# If we already have alpha0 we can use, just skip this function
treatment <- function(prop, alpha, Covar = covar_true){
  
  # prop: proportion of people who receive the treatment
  # alpha: coefficient vector of X, not including alpha0
  # Covar: covariates data, should be a matrix or data frame
  
  alpha1 <- alpha
  
  # Find an alpha0
  eta = 1
  # The initial value is very important
  alpha0 = -5
  
  while(eta > 0.01) {
    
    # Simulate Z
    logit <- alpha1 %*% t(as.matrix(Covar)) + rep(alpha0, nrow(Covar))
    e <- as.vector(exp(logit)/(1+exp(logit)))
    Z <- rbinom(n = nrow(Covar), size = 1, prob = e)
    
    n1 = length(Z[Z==1])
    eta = abs(n1/length(Z) - prop)
    
    # Update the value of alpha0
    alpha0 = alpha0 + 0.01
  }
  PS.data <- cbind(Covar, Z)
  return(list(PS.data, alpha0))
  rm(logit, e)
}

# Generate PS and outcome models
# --- For correctly specified models
PS.model <- function(prop = NA, alpha, alpha0 = NA, Covar = covar_true) {
  
  # First, generating Z
  
  if(is.na(alpha0)) {
    
    if(is.na(prop)) { error("Either prop or alpha0 must be specified!") }
    
    PS <- treatment(prop = prop, alpha = alpha, Covar = Covar)
    PS.data <- PS[[1]]
    alpha0 <- PS[[2]]
    
  } else {
    logit <- alpha1 %*% t(as.matrix(Covar)) + rep(alpha0, nrow(Covar))
    e <- as.vector(exp(logit)/(1+exp(logit)))
    Z <- rbinom(n = nrow(Covar), size = 1, prob = e)
    PS.data <- cbind(Covar, Z)
  }
  return(PS.data)
}

md1.dat <- PS.model(alpha = alpha1, alpha0 = -2.17)
md2.dat <- PS.model(alpha = alpha1, alpha0 = -0.78)
md3.dat <- PS.model(alpha = alpha1, alpha0 = 0.98)

save(file = "md_true.RData", md1.dat, md2.dat, md3.dat)
source("trueEst_func.R")

### Note that the outcome Y is generated in the Heter_eff function
md1.true_eff <- Heter_eff(md1.dat)
md2.true_eff <- Heter_eff(md2.dat)
md3.true_eff <- Heter_eff(md3.dat)

save(file = "truth.Rdata", md1.true_eff, md2.true_eff, md3.true_eff)
