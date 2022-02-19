### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Variance estimations ATE ATT ATC ~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### Truth for ATE, ATT, ATC

### by Yi Liu
### Create date: Oct 09, 2021

Heter_eff <- function(data) {
  
  ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
  fit = glm(formula = ps.mult, data = data, family = binomial(link = "logit"))
  ps = as.numeric(fit$fitted.values)
  
  Y.h <- 4 + 3*(data$X1 + data$X2)^2 + data$X1*data$X3
  
  ATE <- sum(1*Y.h)/nrow(data)
  ATT <- sum(Y.h*ps)/sum(ps)
  ATC <- sum(Y.h*(1-ps))/sum(1-ps)
  
  return(data.frame(True.IPW.h = ATE, True.ATT.h = ATT, True.ATC.h = ATC))
}

# For constant effect, based on our model, we use 4 for all cases
