### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

####### Coding point estimations for all WATEs

### by Yi Liu
### Oct 02, 2021

library(dplyr)

# ATE via IPW weighting
IPW <- function(data, delta = 0) {
  
  ps.mult.c <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
  ps.mult.m <- Z ~ X1 + X2 + X3 + X4
  
  ##################### Hajek like (Weighted) estimators ##################### 
  
  # --- In non-augmented part we only consider correctly specified models
  # --- Propensity score using logistic regression
  fit = glm(formula = ps.mult.c, data = data, family = binomial(link = "logit"))
  ps = as.numeric(fit$fitted.values)

  # --- Trim the dataset by delta which is the threshold of the extreme propensity either too close to 0 or 1
  # --- e.g., delta = 0.05 means we only keep observations with PS ranged in [0.05, 0.95]
  
  data.trim.c <- data %>% mutate(PS = ps) %>% filter(PS>=delta & PS<=(1-delta))
  
  fit = glm(formula = ps.mult.c, data = data.trim.c, family = binomial(link = "logit"))
  data.trim.c$PS = as.numeric(fit$fitted.values)
  
  Y.c = data.trim.c$Y.c
  Y.h = data.trim.c$Y.h
  Z = data.trim.c$Z
  PS = data.trim.c$PS
  
  # Weighted estimators
  mu1 <- sum(Z*Y.c*(1/PS)) / sum(Z/PS)
  mu0 <- sum( (1-Z)*Y.c*(1/(1-PS)) ) / sum((1-Z)/(1-PS))
  IPW.wt.c <- mu1 - mu0
  
  mu1 <- sum(Z*Y.h*(1/PS)) / sum(Z/PS)
  mu0 <- sum( (1-Z)*Y.h*(1/(1-PS)) ) / sum((1-Z)/(1-PS))
  IPW.wt.h <- mu1 - mu0

  ########### Augmented (in this case doubly-robust) estimators #############
  
  # Outcome regression (OR) models: used only for augmented part
  # --- Constant effect and correctly specified
  outform.c.c <- Y.c ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2)
  # --- Heterogeneous effect and correctly specified
  outform.h.c <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2) + I(X1*X3)
  
  # --- Constant effect and misspecified
  outform.c.m <- Y.c ~ X1 + X2 + X3 + X4
  # --- Heterogeneous effect and misspecified
  outform.h.m <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X3)
  
  # Function for getting augmented IPW estimator
  IPW.aug <- function(outform.c, outform.h, data) {
    
    data0 = data[data$Z==0, ]
    data1 = data[data$Z==1, ]
    
    PS = data$PS
    
    Y.c = data$Y.c
    Y.h = data$Y.h
    Z = data$Z
    PS = data$PS
    
    fit = lm(formula = outform.c, data = data0)
    out.c.0 = as.numeric(predict(fit, data))
    fit = lm(formula = outform.c, data = data1)
    out.c.1 = as.numeric(predict(fit, data))
    
    fit = lm(formula = outform.h, data = data0)
    out.h.0 = as.numeric(predict(fit, data))
    fit = lm(formula = outform.h, data = data1)
    out.h.1 = as.numeric(predict(fit, data))
    
    # --- Estimators
    mu1 <- sum( Z*(Y.c-out.c.1)*(1/PS) ) / sum(Z/PS) + sum(out.c.1)/nrow(data)
    mu0 <- sum( (1-Z)*(Y.c-out.c.0)*(1/(1-PS)) ) / sum((1-Z)/(1-PS)) + sum(out.c.0)/nrow(data)
    IPW.aug.c <- mu1 - mu0
    
    mu1 <- sum( Z*(Y.h-out.h.1)*(1/PS) ) / sum(Z/PS) + sum(out.h.1)/nrow(data)
    mu0 <- sum( (1-Z)*(Y.h-out.h.0)*(1/(1-PS)) ) / sum((1-Z)/(1-PS)) + sum(out.h.0)/nrow(data)
    IPW.aug.h <- mu1 - mu0
    
    return(list(aug.c = IPW.aug.c, aug.h = IPW.aug.h))
    
  }
  
  ###### Case 1: Both PS and OR are correctly specified
  aug.est <- IPW.aug(outform.c = outform.c.c, outform.h = outform.h.c, data = data.trim.c)
  IPW.aug.c.cc = aug.est$aug.c
  IPW.aug.h.cc = aug.est$aug.h
  
  ###### Case 2: PS correctly specified, OR misspecified
  aug.est <- IPW.aug(outform.c = outform.c.m, outform.h = outform.h.m, data = data.trim.c)
  IPW.aug.c.cm = aug.est$aug.c
  IPW.aug.h.cm = aug.est$aug.h
  
  ###### Case 3: PS misspecified, OR correctly specified
  fit = glm(formula = ps.mult.m, data = data, family = binomial(link = "logit"))
  ps = as.numeric(fit$fitted.values)
  
  data.trim.m <- data %>% mutate(PS = ps) %>% filter(PS>=delta & PS<=(1-delta))
  
  fit = glm(formula = ps.mult.m, data = data.trim.m, family = binomial(link = "logit"))
  data.trim.m$PS = as.numeric(fit$fitted.values)
  
  aug.est <- IPW.aug(outform.c = outform.c.c, outform.h = outform.h.c, data = data.trim.m)
  IPW.aug.c.mc = aug.est$aug.c
  IPW.aug.h.mc = aug.est$aug.h
  
  ###### Case 4: Both PS and OR misspecified
  aug.est <- IPW.aug(outform.c = outform.c.m, outform.h = outform.h.m, data = data.trim.m)
  IPW.aug.c.mm = aug.est$aug.c
  IPW.aug.h.mm = aug.est$aug.h

  return(data.frame(IPW.wt.c = IPW.wt.c, IPW.wt.h = IPW.wt.h,
                    IPW.aug.c.cc = IPW.aug.c.cc, IPW.aug.h.cc = IPW.aug.h.cc,
                    IPW.aug.c.cm = IPW.aug.c.cm, IPW.aug.h.cm = IPW.aug.h.cm,
                    IPW.aug.c.mc = IPW.aug.c.mc, IPW.aug.h.mc = IPW.aug.h.mc,
                    IPW.aug.c.mm = IPW.aug.c.mm, IPW.aug.h.mm = IPW.aug.h.mm  ))
}


# ATE via overlap weighting
overlap <- function(data) {
  
  ps.mult.c <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
  ps.mult.m <- Z ~ X1 + X2 + X3 + X4
  
  ##################### Hajek like (Weighted) estimators ##################### 
  
  # --- In non-augmented part we only consider correctly specified models
  # --- Propensity score using logistic regression
  fit = glm(formula = ps.mult.c, data = data, family = binomial(link = "logit"))
  data$PS = as.numeric(fit$fitted.values)
  
  Y.c = data$Y.c
  Y.h = data$Y.h
  Z = data$Z
  PS = data$PS
  
  # Weighted estimators
  mu1 <- sum(Z*Y.c*(1-PS)) / sum(Z*(1-PS))
  mu0 <- sum( (1-Z)*Y.c*PS ) / sum((1-Z)*PS)
  ATO.wt.c <- mu1 - mu0
  
  mu1 <- sum(Z*Y.h*(1-PS)) / sum(Z*(1-PS))
  mu0 <- sum( (1-Z)*Y.h*PS ) / sum((1-Z)*PS)
  ATO.wt.h <- mu1 - mu0
  
  ################## Augmented (NOT doubly-robust) estimators ##################
  
  # Outcome regression (OR) models: used only for augmented part
  # --- Constant effect and correctly specified
  outform.c.c <- Y.c ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2) 
  # --- Heterogeneous effect and correctly specified
  outform.h.c <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2) + I(X1*X3)
  
  # --- Constant effect and misspecified
  outform.c.m <- Y.c ~ X1 + X2 + X3 + X4
  # --- Heterogeneous effect and misspecified
  outform.h.m <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X3)
  
  # Function for getting augmented overlap estimator
  ATO.aug <- function(outform.c, outform.h, data) {
    
    data0 = data[data$Z==0, ]
    data1 = data[data$Z==1, ]
    
    Y.c = data$Y.c
    Y.h = data$Y.h
    Z = data$Z
    PS = data$PS
    
    fit = lm(formula = outform.c, data = data0)
    out.c.0 = as.numeric(predict(fit, data))
    fit = lm(formula = outform.c, data = data1)
    out.c.1 = as.numeric(predict(fit, data))
    
    fit = lm(formula = outform.h, data = data0)
    out.h.0 = as.numeric(predict(fit, data))
    fit = lm(formula = outform.h, data = data1)
    out.h.1 = as.numeric(predict(fit, data))
    
    # --- Estimators
    mu1 <- sum( Z*(Y.c-out.c.1)*(1-PS) ) / sum(Z*(1-PS)) + sum(out.c.1*PS*(1-PS))/sum(PS*(1-PS))
    mu0 <- sum( (1-Z)*(Y.c-out.c.0)*PS ) / sum((1-Z)*PS) + sum(out.c.0*PS*(1-PS))/sum(PS*(1-PS))
    ATO.aug.c <- mu1 - mu0
    
    mu1 <- sum( Z*(Y.h-out.h.1)*(1-PS) ) / sum(Z*(1-PS)) + sum(out.h.1*PS*(1-PS))/sum(PS*(1-PS))
    mu0 <- sum( (1-Z)*(Y.h-out.h.0)*PS ) / sum((1-Z)*PS) + sum(out.h.0*PS*(1-PS))/sum(PS*(1-PS))
    ATO.aug.h <- mu1 - mu0
    
    return(list(aug.c = ATO.aug.c, aug.h = ATO.aug.h))
    
  }
  
  ###### Case 1: Both PS and OR correctly specified
  aug.est <- ATO.aug(outform.c = outform.c.c, outform.h = outform.h.c, data = data)
  ATO.aug.c.cc = aug.est$aug.c
  ATO.aug.h.cc = aug.est$aug.h
  
  ###### Case 2: PS correctly specified, OR misspecified
  aug.est <- ATO.aug(outform.c = outform.c.m, outform.h = outform.h.m, data = data)
  ATO.aug.c.cm = aug.est$aug.c
  ATO.aug.h.cm = aug.est$aug.h
  
  ###### Case 3: PS misspecified, OR correctly specified
  fit = glm(formula = ps.mult.m, data = data, family = binomial(link = "logit"))
  data$PS = as.numeric(fit$fitted.values)
  
  aug.est <- ATO.aug(outform.c = outform.c.c, outform.h = outform.h.c, data = data)
  ATO.aug.c.mc = aug.est$aug.c
  ATO.aug.h.mc = aug.est$aug.h
  
  ###### Case 4: Both PS and OR misspecified
  aug.est <- ATO.aug(outform.c = outform.c.m, outform.h = outform.h.m, data = data)
  ATO.aug.c.mm = aug.est$aug.c
  ATO.aug.h.mm = aug.est$aug.h

  return(data.frame(ATO.wt.c = ATO.wt.c, ATO.wt.h = ATO.wt.h,
                    ATO.aug.c.cc = ATO.aug.c.cc, ATO.aug.h.cc = ATO.aug.h.cc,
                    ATO.aug.c.cm = ATO.aug.c.cm, ATO.aug.h.cm = ATO.aug.h.cm,
                    ATO.aug.c.mc = ATO.aug.c.mc, ATO.aug.h.mc = ATO.aug.h.mc,
                    ATO.aug.c.mm = ATO.aug.c.mm, ATO.aug.h.mm = ATO.aug.h.mm ))
}

# ATE via matching weight
matching <- function(data) {
  
  ps.mult.c <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
  ps.mult.m <- Z ~ X1 + X2 + X3 + X4
  
  ##################### Hajek like (Weighted) estimators ##################### 
  
  # --- In non-augmented part we only consider correctly specified models
  # --- Propensity score using logistic regression
  fit = glm(formula = ps.mult.c, data = data, family = binomial(link = "logit"))
  data$PS = as.numeric(fit$fitted.values)
  
  Y.c = data$Y.c
  Y.h = data$Y.h
  Z = data$Z
  PS = data$PS
  
  U = apply(cbind(PS, 1-PS), 1, min)
  
  # Weighted estimators
  mu1 <- sum(Z*Y.c*U*(1/PS)) / sum(U*Z/PS)
  mu0 <- sum( (1-Z)*Y.c*U*(1/(1-PS)) ) / sum(U*(1-Z)/(1-PS))
  ATM.wt.c <- mu1 - mu0
  
  mu1 <- sum(Z*Y.h*U*(1/PS)) / sum(U*Z/PS)
  mu0 <- sum( (1-Z)*Y.h*U*(1/(1-PS)) ) / sum(U*(1-Z)/(1-PS))
  ATM.wt.h <- mu1 - mu0
  
  ################## Augmented (NOT doubly-robust) estimators ###################
  
  # Outcome regression (OR) models: Used only for augmented part
  # --- Constant effect and correctly specified
  outform.c.c <- Y.c ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2)
  # --- Heterogeneous effect and correctly specified
  outform.h.c <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2) + I(X1*X3)
  
  # --- Constant effect and misspecified
  outform.c.m <- Y.c ~ X1 + X2 + X3 + X4
  # --- Heterogeneous effect and misspecified
  outform.h.m <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X3)
  
  # Function for getting augmented matching estimator
  ATM.aug <- function(outform.c, outform.h, data) {
    
    data0 = data[data$Z==0, ]
    data1 = data[data$Z==1, ]
    
    Y.c = data$Y.c
    Y.h = data$Y.h
    Z = data$Z
    PS = data$PS
    
    U = apply(cbind(PS, 1-PS), 1, min)
    
    fit = lm(formula = outform.c, data = data0)
    out.c.0 = as.numeric(predict(fit, data))
    fit = lm(formula = outform.c, data = data1)
    out.c.1 = as.numeric(predict(fit, data))
    
    fit = lm(formula = outform.h, data = data0)
    out.h.0 = as.numeric(predict(fit, data))
    fit = lm(formula = outform.h, data = data1)
    out.h.1 = as.numeric(predict(fit, data))
    
    # --- Estimators
    mu1 <- sum( Z*(Y.c-out.c.1)*U*(1/PS) ) / sum(Z*U/PS) + sum(U*out.c.1)/sum(U)
    mu0 <- sum( (1-Z)*(Y.c-out.c.0)*U*(1/(1-PS)) ) / sum((1-Z)*U/(1-PS)) + sum(U*out.c.0)/sum(U)
    ATM.aug.c <- mu1 - mu0
    
    mu1 <- sum( Z*(Y.h-out.h.1)*U*(1/PS) ) / sum(Z*U/PS) + sum(U*out.h.1)/sum(U)
    mu0 <- sum( (1-Z)*(Y.h-out.h.0)*U*(1/(1-PS)) ) / sum((1-Z)*U/(1-PS)) + sum(U*out.h.0)/sum(U)
    ATM.aug.h <- mu1 - mu0
    
    return(list(aug.c = ATM.aug.c, aug.h = ATM.aug.h))
    
  }
  
  ###### Case 1: Both PS and OR correctly specified
  aug.est <- ATM.aug(outform.c = outform.c.c, outform.h = outform.h.c, data = data)
  ATM.aug.c.cc = aug.est$aug.c
  ATM.aug.h.cc = aug.est$aug.h
  
  ###### Case 2: PS correctly specified, OR misspecified
  aug.est <- ATM.aug(outform.c = outform.c.m, outform.h = outform.h.m, data = data)
  ATM.aug.c.cm = aug.est$aug.c
  ATM.aug.h.cm = aug.est$aug.h
  
  ###### Case 3: PS misspecified, OR correctly specified
  fit = glm(formula = ps.mult.m, data = data, family = binomial(link = "logit"))
  data$PS = as.numeric(fit$fitted.values)
  
  aug.est <- ATM.aug(outform.c = outform.c.c, outform.h = outform.h.c, data = data)
  ATM.aug.c.mc = aug.est$aug.c
  ATM.aug.h.mc = aug.est$aug.h
  
  ###### Case 4: Both PS and OR misspecified
  aug.est <- ATM.aug(outform.c = outform.c.m, outform.h = outform.h.m, data = data)
  ATM.aug.c.mm = aug.est$aug.c
  ATM.aug.h.mm = aug.est$aug.h
  
  return(data.frame(ATM.wt.c = ATM.wt.c, ATM.wt.h = ATM.wt.h,
                    ATM.aug.c.cc = ATM.aug.c.cc, ATM.aug.h.cc = ATM.aug.h.cc,
                    ATM.aug.c.cm = ATM.aug.c.cm, ATM.aug.h.cm = ATM.aug.h.cm,
                    ATM.aug.c.mc = ATM.aug.c.mc, ATM.aug.h.mc = ATM.aug.h.mc,
                    ATM.aug.c.mm = ATM.aug.c.mm, ATM.aug.h.mm = ATM.aug.h.mm ))
}

# ATE via entropy weight
entropy <- function(data) {
  
  ps.mult.c <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
  ps.mult.m <- Z ~ X1 + X2 + X3 + X4
  
  ##################### Hajek like (Weighted) estimators ##################### 
  
  # --- In non-augmented part we only consider correctly specified models
  # --- Propensity score using logistic regression
  fit = glm(formula = ps.mult.c, data = data, family = binomial(link = "logit"))
  data$PS = as.numeric(fit$fitted.values)
  
  Y.c = data$Y.c
  Y.h = data$Y.h
  Z = data$Z
  PS = data$PS
  
  U = -(PS*log(PS) + (1-PS)*log(1-PS))
  
  # Weighted estimators
  mu1 <- sum(Z*Y.c*U*(1/PS)) / sum(U*Z/PS)
  mu0 <- sum( (1-Z)*Y.c*U*(1/(1-PS)) ) / sum(U*(1-Z)/(1-PS))
  ATEN.wt.c <- mu1 - mu0
  
  mu1 <- sum(Z*Y.h*U*(1/PS)) / sum(U*Z/PS)
  mu0 <- sum( (1-Z)*Y.h*U*(1/(1-PS)) ) / sum(U*(1-Z)/(1-PS))
  ATEN.wt.h <- mu1 - mu0
  
  ################## Augmented (NOT doubly-robust) estimators ###################
  
  # Outcome regression (OR) models: Used only for augmented part
  # --- Constant effect and correctly specified
  outform.c.c <- Y.c ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2) 
  # --- Heterogeneous effect and correctly specified
  outform.h.c <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2) + I(X1*X3)
  
  # --- Constant effect and misspecified
  outform.c.m <- Y.c ~ X1 + X2 + X3 + X4
  # --- Heterogeneous effect and misspecified
  outform.h.m <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X3)
  
  # Function for getting augmented entropy estimator
  ATEN.aug <- function(outform.c, outform.h, data) {
    
    data0 = data[data$Z==0, ]
    data1 = data[data$Z==1, ]
    
    Y.c = data$Y.c
    Y.h = data$Y.h
    Z = data$Z
    PS = data$PS
    
    U = -(PS*log(PS) + (1-PS)*log(1-PS))
    
    fit = lm(formula = outform.c, data = data0)
    out.c.0 = as.numeric(predict(fit, data))
    fit = lm(formula = outform.c, data = data1)
    out.c.1 = as.numeric(predict(fit, data))
    
    fit = lm(formula = outform.h, data = data0)
    out.h.0 = as.numeric(predict(fit, data))
    fit = lm(formula = outform.h, data = data1)
    out.h.1 = as.numeric(predict(fit, data))
    
    # --- Estimators
    mu1 <- sum( Z*(Y.c-out.c.1)*U*(1/PS) ) / sum(Z*U/PS) + sum(U*out.c.1)/sum(U)
    mu0 <- sum( (1-Z)*(Y.c-out.c.0)*U*(1/(1-PS)) ) / sum((1-Z)*U/(1-PS)) + sum(U*out.c.0)/sum(U)
    ATEN.aug.c <- mu1 - mu0
    
    mu1 <- sum( Z*(Y.h-out.h.1)*U*(1/PS) ) / sum(Z*U/PS) + sum(U*out.h.1)/sum(U)
    mu0 <- sum( (1-Z)*(Y.h-out.h.0)*U*(1/(1-PS)) ) / sum((1-Z)*U/(1-PS)) + sum(U*out.h.0)/sum(U)
    ATEN.aug.h <- mu1 - mu0
    
    return(list(aug.c = ATEN.aug.c, aug.h = ATEN.aug.h))
    
  }
  
  ###### Case 1: Both PS and OR correctly specified
  aug.est <- ATEN.aug(outform.c = outform.c.c, outform.h = outform.h.c, data = data)
  ATEN.aug.c.cc = aug.est$aug.c
  ATEN.aug.h.cc = aug.est$aug.h
  
  ###### Case 2: PS correctly specified, OR misspecified
  aug.est <- ATEN.aug(outform.c = outform.c.m, outform.h = outform.h.m, data = data)
  ATEN.aug.c.cm = aug.est$aug.c
  ATEN.aug.h.cm = aug.est$aug.h
  
  ###### Case 3: PS misspecified, OR correctly specified
  fit = glm(formula = ps.mult.m, data = data, family = binomial(link = "logit"))
  data$PS = as.numeric(fit$fitted.values)
  
  aug.est <- ATEN.aug(outform.c = outform.c.c, outform.h = outform.h.c, data = data)
  ATEN.aug.c.mc = aug.est$aug.c
  ATEN.aug.h.mc = aug.est$aug.h
  
  ###### Case 4: Both PS and OR misspecified
  aug.est <- ATEN.aug(outform.c = outform.c.m, outform.h = outform.h.m, data = data)
  ATEN.aug.c.mm = aug.est$aug.c
  ATEN.aug.h.mm = aug.est$aug.h
  
  ##### Return values
  return(data.frame(ATEN.wt.c = ATEN.wt.c, ATEN.wt.h = ATEN.wt.h,
                    ATEN.aug.c.cc = ATEN.aug.c.cc, ATEN.aug.h.cc = ATEN.aug.h.cc,
                    ATEN.aug.c.cm = ATEN.aug.c.cm, ATEN.aug.h.cm = ATEN.aug.h.cm,
                    ATEN.aug.c.mc = ATEN.aug.c.mc, ATEN.aug.h.mc = ATEN.aug.h.mc,
                    ATEN.aug.c.mm = ATEN.aug.c.mm, ATEN.aug.h.mm = ATEN.aug.h.mm ))
}

# ATE on controls
control <- function(data) {
  
  ps.mult.c <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
  ps.mult.m <- Z ~ X1 + X2 + X3 + X4
  
  ##################### Hajek like (Weighted) estimators ##################### 
  
  # --- In non-augmented part we only consider correctly specified models
  # --- Propensity score using logistic regression
  fit = glm(formula = ps.mult.c, data = data, family = binomial(link = "logit"))
  data$PS = as.numeric(fit$fitted.values)

  Y.c = data$Y.c
  Y.h = data$Y.h
  Z = data$Z
  PS = data$PS
  
  # Weighted estimators
  mu1 <- sum(Z*Y.c*(1-PS)/PS) / sum((1-PS)*Z/PS)
  mu0 <- sum( (1-Z)*Y.c ) / sum(1-Z)
  ATC.wt.c <- mu1 - mu0
  
  mu1 <- sum(Z*Y.h*(1-PS)/PS) / sum((1-PS)*Z/PS)
  mu0 <- sum( (1-Z)*Y.h ) / sum(1-Z)
  ATC.wt.h <- mu1 - mu0
  
  
  ################ Augmented (NOT doubly-robust) estimators #################
  
  # Outcome regression (OR) models: Used only for augmented part
  # --- Constant effect and correctly specified
  outform.c.c <- Y.c ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2)
  # --- Heterogeneous effect and correctly specified
  outform.h.c <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2) + I(X1*X3)
  
  # --- Constant effect and misspecified
  outform.c.m <- Y.c ~ X1 + X2 + X3 + X4
  # --- Heterogeneous effect and misspecified
  outform.h.m <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X3)
  
  # Function for getting augmented ATC estimator
  ATC.aug <- function(outform.c, outform.h, data) {
    
    data0 = data[data$Z==0, ]
    data1 = data[data$Z==1, ]
    
    Y.c = data$Y.c
    Y.h = data$Y.h
    Z = data$Z
    PS = data$PS
    
    fit = lm(formula = outform.c, data = data0)
    out.c.0 = as.numeric(predict(fit, data))
    fit = lm(formula = outform.c, data = data1)
    out.c.1 = as.numeric(predict(fit, data))
    
    fit = lm(formula = outform.h, data = data0)
    out.h.0 = as.numeric(predict(fit, data))
    fit = lm(formula = outform.h, data = data1)
    out.h.1 = as.numeric(predict(fit, data))
    
    # --- Estimators
    mu1 <- sum( Z*(Y.c-out.c.1)*(1-PS)/PS ) / sum(Z*(1-PS)/PS) + sum((1-PS)*out.c.1)/sum(1-PS)
    mu0 <- sum( (1-Z)*(Y.c-out.c.0)) / sum(1-Z) + sum((1-PS)*out.c.0)/sum(1-PS)
    ATC.aug.c <- mu1 - mu0
    
    mu1 <- sum( Z*(Y.h-out.h.1)*(1-PS)/PS ) / sum(Z*(1-PS)/PS) + sum((1-PS)*out.h.1)/sum(1-PS)
    mu0 <- sum( (1-Z)*(Y.h-out.h.0) ) / sum(1-Z) + sum((1-PS)*out.h.0)/sum(1-PS)
    ATC.aug.h <- mu1 - mu0
    
    return(list(aug.c = ATC.aug.c, aug.h = ATC.aug.h))
    
  }
  
  ###### Case 1: Both PS and OR correctly specified
  aug.est <- ATC.aug(outform.c = outform.c.c, outform.h = outform.h.c, data = data)
  ATC.aug.c.cc = aug.est$aug.c
  ATC.aug.h.cc = aug.est$aug.h
  
  ###### Case 2: Both PS correctly specified, OR misspecified
  aug.est <- ATC.aug(outform.c = outform.c.m, outform.h = outform.h.m, data = data)
  ATC.aug.c.cm = aug.est$aug.c
  ATC.aug.h.cm = aug.est$aug.h
  
  ###### Case 3: PS misspecified, OR correctly specified
  fit = glm(formula = ps.mult.m, data = data, family = binomial(link = "logit"))
  data$PS  = as.numeric(fit$fitted.values)
  
  aug.est <- ATC.aug(outform.c = outform.c.c, outform.h = outform.h.c, data = data)
  ATC.aug.c.mc = aug.est$aug.c
  ATC.aug.h.mc = aug.est$aug.h
  
  ###### Case 4: Both PS and OR misspecified
  aug.est <- ATC.aug(outform.c = outform.c.m, outform.h = outform.h.m, data = data)
  ATC.aug.c.mm = aug.est$aug.c
  ATC.aug.h.mm = aug.est$aug.h

  return(data.frame(ATC.wt.c = ATC.wt.c, ATC.wt.h = ATC.wt.h,
                    ATC.aug.c.cc = ATC.aug.c.cc, ATC.aug.h.cc = ATC.aug.h.cc,
                    ATC.aug.c.cm = ATC.aug.c.cm, ATC.aug.h.cm = ATC.aug.h.cm,
                    ATC.aug.c.mc = ATC.aug.c.mc, ATC.aug.h.mc = ATC.aug.h.mc,
                    ATC.aug.c.mm = ATC.aug.c.mm, ATC.aug.h.mm = ATC.aug.h.mm ))
}

# ATE on treated
treated <- function(data) {
  
  ps.mult.c <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
  ps.mult.m <- Z ~ X1 + X2 + X3 + X4
  
  ##################### Hajek like (Weighted) estimators ##################### 
  
  # --- In non-augmented part we only consider correctly specified models
  # --- Propensity score using logistic regression
  fit = glm(formula = ps.mult.c, data = data, family = binomial(link = "logit"))
  data$PS = as.numeric(fit$fitted.values)
  
  Y.c = data$Y.c
  Y.h = data$Y.h
  Z = data$Z
  PS = data$PS
  
  # Weighted estimators
  mu1 <- sum(Z*Y.c) / sum(Z)
  mu0 <- sum( (1-Z)*Y.c*PS/(1-PS) ) / sum((1-Z)*PS/(1-PS))
  ATT.wt.c <- mu1 - mu0
  
  mu1 <- sum(Z*Y.h) / sum(Z)
  mu0 <- sum( (1-Z)*Y.h*PS/(1-PS) ) / sum((1-Z)*PS/(1-PS))
  ATT.wt.h <- mu1 - mu0
  
  
  ##################### Augmented (NOT doubly-robust) estimators ######################
  
  # Outcome regression (OR) models: Used only for augmented part
  # --- Constant effect and correctly specified
  outform.c.c <- Y.c ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2) 
  # --- Heterogeneous effect and correctly specified
  outform.h.c <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2) + I(X1*X3)
  
  # --- Constant effect and misspecified
  outform.c.m <- Y.c ~ X1 + X2 + X3 + X4
  # --- Heterogeneous effect and misspecified
  outform.h.m <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X3)
  
  # Function for getting augmented IPW estimator
  ATT.aug <- function(outform.c, outform.h, data) {
    
    data0 = data[data$Z==0, ]
    data1 = data[data$Z==1, ]
    
    Y.c = data$Y.c
    Y.h = data$Y.h
    Z = data$Z
    PS = data$PS
    
    fit = lm(formula = outform.c, data = data0)
    out.c.0 = as.numeric(predict(fit, data))
    fit = lm(formula = outform.c, data = data1)
    out.c.1 = as.numeric(predict(fit, data))
    
    fit = lm(formula = outform.h, data = data0)
    out.h.0 = as.numeric(predict(fit, data))
    fit = lm(formula = outform.h, data = data1)
    out.h.1 = as.numeric(predict(fit, data))
    
    # --- Estimators
    mu1 <- sum( Z*(Y.c-out.c.1) ) / sum(Z) + sum(PS*out.c.1)/sum(PS)
    mu0 <- sum( (1-Z)*(Y.c-out.c.0)*PS/(1-PS)) / sum((1-Z)*PS/(1-PS)) + sum(PS*out.c.0)/sum(PS)
    ATT.aug.c <- mu1 - mu0
    
    mu1 <- sum( Z*(Y.h-out.h.1) ) / sum(Z) + sum(PS*out.h.1)/sum(PS)
    mu0 <- sum( (1-Z)*(Y.h-out.h.0)*PS/(1-PS) ) / sum((1-Z)*PS/(1-PS)) + sum(PS*out.h.0)/sum(PS)
    ATT.aug.h <- mu1 - mu0
    
    return(list(aug.c = ATT.aug.c, aug.h = ATT.aug.h))
    
  }
  
  ###### Case 1: Both PS and OR correctly specified
  aug.est <- ATT.aug(outform.c = outform.c.c, outform.h = outform.h.c, data = data)
  ATT.aug.c.cc = aug.est$aug.c
  ATT.aug.h.cc = aug.est$aug.h
  
  ###### Case 2: PS correctly specified, OR misspecified
  aug.est <- ATT.aug(outform.c = outform.c.m, outform.h = outform.h.m, data = data)
  ATT.aug.c.cm = aug.est$aug.c
  ATT.aug.h.cm = aug.est$aug.h
  
  ###### Case 3: PS misspecified, OR correctly specified
  fit = glm(formula = ps.mult.m, data = data, family = binomial(link = "logit"))
  data$PS = as.numeric(fit$fitted.values)
  
  aug.est <- ATT.aug(outform.c = outform.c.c, outform.h = outform.h.c, data = data)
  ATT.aug.c.mc = aug.est$aug.c
  ATT.aug.h.mc = aug.est$aug.h
  
  ###### Case 4: Both PS and OR models are misspecified
  aug.est <- ATT.aug(outform.c = outform.c.m, outform.h = outform.h.m, data = data)
  ATT.aug.c.mm = aug.est$aug.c
  ATT.aug.h.mm = aug.est$aug.h
  
  return(data.frame(ATT.wt.c = ATT.wt.c, ATT.wt.h = ATT.wt.h,
                    ATT.aug.c.cc = ATT.aug.c.cc, ATT.aug.h.cc = ATT.aug.h.cc,
                    ATT.aug.c.cm = ATT.aug.c.cm, ATT.aug.h.cm = ATT.aug.h.cm,
                    ATT.aug.c.mc = ATT.aug.c.mc, ATT.aug.h.mc = ATT.aug.h.mc,
                    ATT.aug.c.mm = ATT.aug.c.mm, ATT.aug.h.mm = ATT.aug.h.mm ))
}
