### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### ~~~~ Truth of (heterogeneous) WATE

### ~~~~ by Yi Liu
### ~~~~ Oct 09, 2021

# --- We get the truth of heterogenous WATE using the following function
# --- The data input has a very large size, and Z is correctly specified

Heter_eff <- function(data) {
  
  ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
  fit = glm(formula = ps.mult, data = data, family = binomial(link = "logit"))
  ps = as.numeric(fit$fitted.values)
  
  Y.h <- 4 + 3*(data$X1 + data$X2)^2 + data$X1*data$X3
  
  IPW <- sum(1*Y.h)/nrow(data)
  IPW.5 <- sum(Y.h[ps >= 0.05 & ps <= 0.95]*1)/sum(ps >= 0.05 & ps <= 0.95)
  IPW.10 <- sum(Y.h[ps >= 0.1 & ps <= 0.9]*1)/sum(ps >= 0.1 & ps <= 0.9)
  IPW.15 <- sum(Y.h[ps >= 0.15 & ps <= 0.85]*1)/sum(ps >= 0.15 & ps <= 0.85)
  
  ATO <- sum(Y.h*ps*(1-ps))/sum(ps*(1-ps))
  
  u = apply(cbind(ps, 1-ps), 1, min)
  ATM <- sum(Y.h*u)/sum(u)
  
  ksi = -(ps*log(ps) + (1-ps)*log(1-ps))
  ATEN <- sum(Y.h*ksi)/sum(ksi)
  
  ATT <- sum(Y.h*ps)/sum(ps)
  ATC <- sum(Y.h*(1-ps))/sum(1-ps)
  
  return(data.frame(True.IPW.h = IPW, True.IPW.5.h = IPW.5, True.IPW.10.h = IPW.10,
                    True.IPW.15.h = IPW.15, True.ATO.h = ATO, True.ATM.h = ATM,
                    True.ATEN.h = ATEN, True.ATC.h = ATC, True.ATT.h = ATT))
}

# Heter_eff(md1.dat)
# Heter_eff(md2.dat)
# Heter_eff(md3.dat)
# Heter_eff(md4.dat)
# Heter_eff(md5.dat)
