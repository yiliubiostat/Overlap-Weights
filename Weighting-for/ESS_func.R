### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

####### Effective sample sizes (ESS)

### by Yi Liu
### Oct 9, 2021

library(PSweight)
library(dplyr)

ESS <- function(md_sim.data) {
  
  ess.IPW.1 <- c()
  ess.IPW.0 <- c()
  
  ess.IPW.5.1 <- c()
  ess.IPW.5.0 <- c()
  
  ess.IPW.10.1 <- c()
  ess.IPW.10.0 <- c()
  
  ess.IPW.15.1 <- c()
  ess.IPW.15.0 <- c()
  
  ess.ATO.1 <- c()
  ess.ATO.0 <- c()
  
  ess.ATM.1 <- c()
  ess.ATM.0 <- c()
  
  ess.ATEN.1 <- c()
  ess.ATEN.0 <- c()
  
  ess.ATC.1 <- c()
  ess.ATC.0 <- c()
  
  ess.ATT.1 <- c()
  ess.ATT.0 <- c()
  
  R <- unique(md_sim.data$Ite)
  for(i in R) {
    
    data <- md_sim.data %>% filter(Ite == i)
    
    ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
    fit = glm(formula = ps.mult, data = data, family = binomial(link = "logit"))
    ps = as.numeric(fit$fitted.values)
    Z = data$Z
    
    ess.IPW.1 <- c(ess.IPW.1, (sum(Z*(1/ps)))^2 / sum((Z*(1/ps))^2) )
    ess.IPW.0 <- c(ess.IPW.0, (sum((1-Z)*(1/(1-ps))))^2 / sum(((1-Z)*(1/(1-ps)))^2) )
    
    ess.IPW.5.1 <- c(ess.IPW.5.1, (sum(Z[ps >= 0.05 & ps <= 0.95]*(1/ps[ps >= 0.05 & ps <= 0.95])))^2 / sum((Z[ps >= 0.05 & ps <= 0.95]*(1/ps[ps >= 0.05 & ps <= 0.95]))^2) )
    ess.IPW.5.0 <- c(ess.IPW.5.0, (sum((1-Z[ps >= 0.05 & ps <= 0.95])*(1/(1-ps[ps >= 0.05 & ps <= 0.95]))))^2 / sum(((1-Z[ps >= 0.05 & ps <= 0.95])*(1/(1-ps[ps >= 0.05 & ps <= 0.95])))^2) )
    
    ess.IPW.10.1 <- c(ess.IPW.10.1, (sum(Z[ps >= 0.1 & ps <= 0.9]*(1/ps[ps >= 0.1 & ps <= 0.9])))^2 / sum((Z[ps >= 0.1 & ps <= 0.9]*(1/ps[ps >= 0.1 & ps <= 0.9]))^2) )
    ess.IPW.10.0 <- c(ess.IPW.10.0, (sum((1-Z[ps >= 0.1 & ps <= 0.9])*(1/(1-ps[ps >= 0.1 & ps <= 0.9]))))^2 / sum(((1-Z[ps >= 0.1 & ps <= 0.9])*(1/(1-ps[ps >= 0.1 & ps <= 0.9])))^2) )
    
    ess.IPW.15.1 <- c(ess.IPW.15.1, (sum(Z[ps >= 0.15 & ps <= 0.85]*(1/ps[ps >= 0.15 & ps <= 0.85])))^2 / sum((Z[ps >= 0.15 & ps <= 0.85]*(1/ps[ps >= 0.15 & ps <= 0.85]))^2) )
    ess.IPW.15.0 <- c(ess.IPW.15.0, (sum((1-Z[ps >= 0.15 & ps <= 0.85])*(1/(1-ps[ps >= 0.15 & ps <= 0.85]))))^2 / sum(((1-Z[ps >= 0.15 & ps <= 0.85])*(1/(1-ps[ps >= 0.15 & ps <= 0.85])))^2) )
    
    ess.ATO.1 <- c(ess.ATO.1, (sum(Z*(1-ps)))^2 / sum((Z*(1-ps))^2) )
    ess.ATO.0 <- c(ess.ATO.0, (sum((1-Z)*ps))^2 / sum(((1-Z)*ps)^2) )
    
    u = apply(cbind(ps, 1-ps), 1, min)
    ess.ATM.1 <- c(ess.ATM.1, (sum(Z*u*(1/ps)))^2 / sum((Z*u*(1/ps))^2) )
    ess.ATM.0 <- c(ess.ATM.0, (sum((1-Z)*u*(1/(1-ps))))^2 / sum(((1-Z)*u*(1/(1-ps)))^2) )
    
    ksi = -(ps*log(ps) + (1-ps)*log(1-ps))
    ess.ATEN.1 <- c(ess.ATEN.1, (sum(Z*ksi*(1/ps)))^2 / sum((Z*ksi*(1/ps))^2) )
    ess.ATEN.0 <- c(ess.ATEN.0, (sum((1-Z)*ksi*(1/(1-ps))))^2 / sum(((1-Z)*ksi*(1/(1-ps)))^2) )
    
    ess.ATC.1 <- c(ess.ATC.1, (sum(Z*((1-ps)/ps)))^2 / sum((Z*((1-ps)/ps))^2) )
    ess.ATC.0 <- c(ess.ATC.0, (sum((1-Z)*1))^2 / sum(((1-Z)*1)^2) )
    
    ess.ATT.1 <- c(ess.ATT.1, (sum(Z*1))^2 / sum((Z*1)^2) )
    ess.ATT.0 <- c(ess.ATT.0, (sum((1-Z)*(ps/(1-ps))))^2 / sum(((1-Z)*(ps/(1-ps)))^2) )
    
  }
  
  ESS <- data.frame(Z = c(0,1), 
                    IPW = c(mean(ess.IPW.0), mean(ess.IPW.1)), 
                    IPW.5 = c(mean(ess.IPW.5.0), mean(ess.IPW.5.1)), 
                    IPW.10 = c(mean(ess.IPW.10.0), mean(ess.IPW.10.1)),
                    IPW.15 = c(mean(ess.IPW.15.0), mean(ess.IPW.15.1)), 
                    ATO = c(mean(ess.ATO.0), mean(ess.ATO.1)), 
                    ATM = c(mean(ess.ATM.0), mean(ess.ATM.1)), 
                    ATEN = c(mean(ess.ATEN.0), mean(ess.ATEN.1)), 
                    ATC = c(mean(ess.ATC.0), mean(ess.ATC.1)), 
                    ATT = c(mean(ess.ATT.0), mean(ess.ATT.1)) )
  return(ESS)
}

# For example,
# load("md1_sims.rdata")
# ESS(weight.data.reps1)
