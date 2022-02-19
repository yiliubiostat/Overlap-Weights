### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

####### Effective sample sizes (ESS)

### by Yi Liu
### Oct 9, 2021

library(PSweight)
library(dplyr)

# Function for returning effective sample size in this study
ESS <- function(md_sim.data) {
  
  ess.IPW.1 <- c()
  ess.IPW.0 <- c()
  
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
    
    ess.ATC.1 <- c(ess.ATC.1, (sum(Z*((1-ps)/ps)))^2 / sum((Z*((1-ps)/ps))^2) )
    ess.ATC.0 <- c(ess.ATC.0, (sum((1-Z)*1))^2 / sum(((1-Z)*1)^2) )
    
    ess.ATT.1 <- c(ess.ATT.1, (sum(Z*1))^2 / sum((Z*1)^2) )
    ess.ATT.0 <- c(ess.ATT.0, (sum((1-Z)*(ps/(1-ps))))^2 / sum(((1-Z)*(ps/(1-ps)))^2) )
  }
  ESS <- data.frame(Z = c(0,1), 
                    IPW = c(mean(ess.IPW.0), mean(ess.IPW.1)), 
                    ATC = c(mean(ess.ATC.0), mean(ess.ATC.1)), 
                    ATT = c(mean(ess.ATT.0), mean(ess.ATT.1)) )
  
  return(ESS) 
}

# For example,
# load("md1_sims.rdata")
# ESS(weight.data.reps1)
