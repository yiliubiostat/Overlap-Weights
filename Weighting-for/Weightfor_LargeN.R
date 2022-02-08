### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### ~~~~ Simulate a large population

### by Yi Liu
### Sept 5, 2021

# --- The purpose of this program is to calculate (heterogeneous) treatment effect truth
# --- using the very large N as follows
N_true <- 1000000

# Generate covariates
set.seed(13244)

# --- X3 and X4: binary
X4 <- rbinom(n = N_true, size = 1, prob = 0.5)
X3 <- rbinom(n = N_true, size = 1, prob = 0.6*X4+0.4*(1-X4))

# --- X1 and X2: joint Gaussian, continuous
library(MASS)
mu <- cbind(-X3 + X4 + 0.5*X3*X4, X3 - X4 + X3*X4)

A1 <- matrix(c(1, .5, .5, 1), ncol = 2)
A2 <- matrix(c(2, .25, .25, 2), ncol = 2)

BN <- data.frame(X1 = rep(NA, N_true), X2 = rep(NA, N_true))

for(i in 1:N_true) {
  BN[i,] = mvrnorm(n = 1, mu = mu[i,], Sigma = X3[i]*A1 + (1-X3[i])*A2)
}

# --- Interactions X5~X7
X5 <- (BN$X1)^2
X6 <- (BN$X1)*(BN$X2)
X7 <- (BN$X2)^2

covar_true <- data.frame(X1 = BN$X1, X2 = BN$X2,
                         X3 = X3, X4 = X4, X5 = X5, X6 = X6, X7 = X7)
rm(BN, A1, A2)

# save(file = "covar_true.RData", covar_true)
# load("covar_true.RData")

# A function for finding alpha0 for varying proportions of treated subjects
# If we already have alpha0 we can use, there is no need to run this function
treatment <- function(prop, alpha, Covar = covar_true){
  
  # prop: proportion of people who receive the treatment
  # alpha: coefficient vector of X, not including alpha0
  # Covar: covariates data, should be a matrix or data frame
  
  alpha1 <- alpha
  
  # Find an alpha0
  eta = 1
  # The initial value is important, this is a monotone selection process
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

PS.model <- function(prop = NA, alpha, alpha0 = NA, Covar = covar_true) {
  
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

alpha <- c(0.3, 0.4, 0.4, 0.4, -0.1, -0.1, 0.1)

# These alpha0's are chosen by the iterative loop above
md1.dat <- PS.model(alpha = alpha, alpha0 = -3.07)
md2.dat <- PS.model(alpha = alpha, alpha0 = -2.17)
md3.dat <- PS.model(alpha = alpha, alpha0 = -0.78)
md4.dat <- PS.model(alpha = alpha, alpha0 = 0.98)
md5.dat <- PS.model(alpha = alpha, alpha0 = 1.86)

# Y will be generated in "TrueEst.R" program
save(file = "md_true.RData", md1.dat, md2.dat, md3.dat, md4.dat, md5.dat)
