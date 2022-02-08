### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

####### Model 1 simulation study

### by Yi Liu
### Oct 9, 2021

load("md_true.RData")
source("Weightfor_PS_OR.R")
source("TrueEst.R")
source("newSand_func.R")

alpha1 <- c(0.3, 0.4, 0.4, 0.4, -0.1, -0.1, 0.1)
md1.alpha0 <- -3.07
weight.data.reps <- PS.model.reps(alpha0 = md1.alpha0, alpha1)

# Initialize the data frame for recording results
results.c <- data.frame(IPW.c = rep(NA, M),
                        IPW.c.var = rep(NA, M),
                        IPW.c.lwr = rep(NA, M),
                        IPW.c.upr = rep(NA, M),
                        IPW.c.ifci = rep(NA, M),
                        
                        IPW.5.c = rep(NA, M),
                        IPW.5.c.var = rep(NA, M),
                        IPW.5.c.lwr = rep(NA, M),
                        IPW.5.c.upr = rep(NA, M),
                        IPW.5.c.ifci = rep(NA, M),
                        
                        IPW.10.c = rep(NA, M),
                        IPW.10.c.var = rep(NA, M),
                        IPW.10.c.lwr = rep(NA, M),
                        IPW.10.c.upr = rep(NA, M),
                        IPW.10.c.ifci = rep(NA, M),
                        
                        IPW.15.c = rep(NA, M),
                        IPW.15.c.var = rep(NA, M),
                        IPW.15.c.lwr = rep(NA, M),
                        IPW.15.c.upr = rep(NA, M),
                        IPW.15.c.ifci = rep(NA, M),
                        
                        ATO.c = rep(NA, M),
                        ATO.c.var = rep(NA, M),
                        ATO.c.lwr = rep(NA, M),
                        ATO.c.upr = rep(NA, M),
                        ATO.c.ifci = rep(NA, M),
                        
                        ATM.c = rep(NA, M),
                        ATM.c.var = rep(NA, M),
                        ATM.c.lwr = rep(NA, M),
                        ATM.c.upr = rep(NA, M),
                        ATM.c.ifci = rep(NA, M),
                        
                        ATEN.c = rep(NA, M),
                        ATEN.c.var = rep(NA, M),
                        ATEN.c.lwr = rep(NA, M),
                        ATEN.c.upr = rep(NA, M),
                        ATEN.c.ifci = rep(NA, M),
                        
                        ATC.c = rep(NA, M),
                        ATC.c.var = rep(NA, M),
                        ATC.c.lwr = rep(NA, M),
                        ATC.c.upr = rep(NA, M),
                        ATC.c.ifci = rep(NA, M),
                        
                        ATT.c = rep(NA, M),
                        ATT.c.var = rep(NA, M),
                        ATT.c.lwr = rep(NA, M),
                        ATT.c.upr = rep(NA, M),
                        ATT.c.ifci = rep(NA, M),
                        
                        Row.Num = 1:M)

results.h <- data.frame(IPW.h = rep(NA, M),
                        IPW.h.var = rep(NA, M),
                        IPW.h.lwr = rep(NA, M),
                        IPW.h.upr = rep(NA, M),
                        IPW.h.ifci = rep(NA, M),
                        
                        IPW.5.h = rep(NA, M),
                        IPW.5.h.var = rep(NA, M),
                        IPW.5.h.lwr = rep(NA, M),
                        IPW.5.h.upr = rep(NA, M),
                        IPW.5.h.ifci = rep(NA, M),
                        
                        IPW.10.h = rep(NA, M),
                        IPW.10.h.var = rep(NA, M),
                        IPW.10.h.lwr = rep(NA, M),
                        IPW.10.h.upr = rep(NA, M),
                        IPW.10.h.ifci = rep(NA, M),
                        
                        IPW.15.h = rep(NA, M),
                        IPW.15.h.var = rep(NA, M),
                        IPW.15.h.lwr = rep(NA, M),
                        IPW.15.h.upr = rep(NA, M),
                        IPW.15.h.ifci = rep(NA, M),
                        
                        ATO.h = rep(NA, M),
                        ATO.h.var = rep(NA, M),
                        ATO.h.lwr = rep(NA, M),
                        ATO.h.upr = rep(NA, M),
                        ATO.h.ifci = rep(NA, M),
                        
                        ATM.h = rep(NA, M),
                        ATM.h.var = rep(NA, M),
                        ATM.h.lwr = rep(NA, M),
                        ATM.h.upr = rep(NA, M),
                        ATM.h.ifci = rep(NA, M),
                        
                        ATEN.h = rep(NA, M),
                        ATEN.h.var = rep(NA, M),
                        ATEN.h.lwr = rep(NA, M),
                        ATEN.h.upr = rep(NA, M),
                        ATEN.h.ifci = rep(NA, M),
                        
                        ATC.h = rep(NA, M),
                        ATC.h.var = rep(NA, M),
                        ATC.h.lwr = rep(NA, M),
                        ATC.h.upr = rep(NA, M),
                        ATC.h.ifci = rep(NA, M),
                        
                        ATT.h = rep(NA, M),
                        ATT.h.var = rep(NA, M),
                        ATT.h.lwr = rep(NA, M),
                        ATT.h.upr = rep(NA, M),
                        ATT.h.ifci = rep(NA, M),
                        
                        Row.Num = 1:M)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~ Weighted/Non-augmented Estimators ~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim.wt.c <- results.c
sim.wt.h <- results.h

ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
atc.mult <- ZZ ~ X1 + X2 + X3 + X4 + X5 + X6 + X7

PE.c <- 4
PE.h <- Heter_eff(md1.dat)

for(i in 1:M) {
  data <- weight.data.reps %>% filter(Ite == i)
  
  # Constant treatment effect
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = data), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.wt.c$IPW.c[i] <- PE
  sim.wt.c$IPW.c.var[i] <- Var
  sim.wt.c$IPW.c.lwr[i] <- Lower
  sim.wt.c$IPW.c.upr[i] <- Upper
  sim.wt.c$IPW.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- IPW with 0.05 trimming
  trim1 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.05)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = trim1$data), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.wt.c$IPW.5.c[i] <- PE
  sim.wt.c$IPW.5.c.var[i] <- Var
  sim.wt.c$IPW.5.c.lwr[i] <- Lower
  sim.wt.c$IPW.5.c.upr[i] <- Upper
  sim.wt.c$IPW.5.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- IPW with 0.10 trimming
  trim2 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.10)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = trim2$data), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.wt.c$IPW.10.c[i] <- PE
  sim.wt.c$IPW.10.c.var[i] <- Var
  sim.wt.c$IPW.10.c.lwr[i] <- Lower
  sim.wt.c$IPW.10.c.upr[i] <- Upper
  sim.wt.c$IPW.10.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- IPW with 0.15 trimming
  trim3 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.15)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = trim3$data), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.wt.c$IPW.15.c[i] <- PE
  sim.wt.c$IPW.15.c.var[i] <- Var
  sim.wt.c$IPW.15.c.lwr[i] <- Lower
  sim.wt.c$IPW.15.c.upr[i] <- Upper
  sim.wt.c$IPW.15.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATO
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "overlap", data = data), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.wt.c$ATO.c[i] <- PE
  sim.wt.c$ATO.c.var[i] <- Var
  sim.wt.c$ATO.c.lwr[i] <- Lower
  sim.wt.c$ATO.c.upr[i] <- Upper
  sim.wt.c$ATO.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATM
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "matching", data = data), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.wt.c$ATM.c[i] <- PE
  sim.wt.c$ATM.c.var[i] <- Var
  sim.wt.c$ATM.c.lwr[i] <- Lower
  sim.wt.c$ATM.c.upr[i] <- Upper
  sim.wt.c$ATM.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATEN
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "entropy", data = data), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.wt.c$ATEN.c[i] <- PE
  sim.wt.c$ATEN.c.var[i] <- Var
  sim.wt.c$ATEN.c.lwr[i] <- Lower
  sim.wt.c$ATEN.c.upr[i] <- Upper
  sim.wt.c$ATEN.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.c", weight = "treated", data = data), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  sim.wt.c$ATC.c[i] <- PE
  sim.wt.c$ATC.c.var[i] <- Var
  sim.wt.c$ATC.c.lwr[i] <- Lower
  sim.wt.c$ATC.c.upr[i] <- Upper
  sim.wt.c$ATC.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "treated", data = data), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.wt.c$ATT.c[i] <- PE
  sim.wt.c$ATT.c.var[i] <- Var
  sim.wt.c$ATT.c.lwr[i] <- Lower
  sim.wt.c$ATT.c.upr[i] <- Upper
  sim.wt.c$ATT.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # Heterogeneous treatment effect
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = data), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.wt.h$IPW.h[i] <- PE
  sim.wt.h$IPW.h.var[i] <- Var
  sim.wt.h$IPW.h.lwr[i] <- Lower
  sim.wt.h$IPW.h.upr[i] <- Upper
  sim.wt.h$IPW.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.h < Upper & PE.h$True.IPW.h > Lower), 1, 0))
  
  # --- IPW with 0.05 trimming
  trim1 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.05)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = trim1$data), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.wt.h$IPW.5.h[i] <- PE
  sim.wt.h$IPW.5.h.var[i] <- Var
  sim.wt.h$IPW.5.h.lwr[i] <- Lower
  sim.wt.h$IPW.5.h.upr[i] <- Upper
  sim.wt.h$IPW.5.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.5.h < Upper & PE.h$True.IPW.5.h > Lower), 1, 0))
  
  # --- IPW with 0.10 trimming
  trim2 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.10)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = trim2$data), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.wt.h$IPW.10.h[i] <- PE
  sim.wt.h$IPW.10.h.var[i] <- Var
  sim.wt.h$IPW.10.h.lwr[i] <- Lower
  sim.wt.h$IPW.10.h.upr[i] <- Upper
  sim.wt.h$IPW.10.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.10.h < Upper & PE.h$True.IPW.10.h > Lower), 1, 0))
  
  # --- IPW with 0.15 trimming
  trim3 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.15)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = trim3$data), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.wt.h$IPW.15.h[i] <- PE
  sim.wt.h$IPW.15.h.var[i] <- Var
  sim.wt.h$IPW.15.h.lwr[i] <- Lower
  sim.wt.h$IPW.15.h.upr[i] <- Upper
  sim.wt.h$IPW.15.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.15.h < Upper & PE.h$True.IPW.15.h > Lower), 1, 0))
  
  # --- ATO
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "overlap", data = data), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.wt.h$ATO.h[i] <- PE
  sim.wt.h$ATO.h.var[i] <- Var
  sim.wt.h$ATO.h.lwr[i] <- Lower
  sim.wt.h$ATO.h.upr[i] <- Upper
  sim.wt.h$ATO.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATO.h < Upper & PE.h$True.ATO.h > Lower), 1, 0))
  
  # --- ATM
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "matching", data = data), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.wt.h$ATM.h[i] <- PE
  sim.wt.h$ATM.h.var[i] <- Var
  sim.wt.h$ATM.h.lwr[i] <- Lower
  sim.wt.h$ATM.h.upr[i] <- Upper
  sim.wt.h$ATM.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATM.h < Upper & PE.h$True.ATM.h > Lower), 1, 0))
  
  # --- ATEN
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "entropy", data = data), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.wt.h$ATEN.h[i] <- PE
  sim.wt.h$ATEN.h.var[i] <- Var
  sim.wt.h$ATEN.h.lwr[i] <- Lower
  sim.wt.h$ATEN.h.upr[i] <- Upper
  sim.wt.h$ATEN.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATEN.h < Upper & PE.h$True.ATEN.h > Lower), 1, 0))
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.h", weight = "treated", data = data), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  sim.wt.h$ATC.h[i] <- PE
  sim.wt.h$ATC.h.var[i] <- Var
  sim.wt.h$ATC.h.lwr[i] <- Lower
  sim.wt.h$ATC.h.upr[i] <- Upper
  sim.wt.h$ATC.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "treated", data = data), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.wt.h$ATT.h[i] <- PE
  sim.wt.h$ATT.h.var[i] <- Var
  sim.wt.h$ATT.h.lwr[i] <- Lower
  sim.wt.h$ATT.h.upr[i] <- Upper
  sim.wt.h$ATT.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~ Augmented Estimators ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Case 1: Both PS and OR models are correctly specified
sim.aug.c.cc <- results.c
sim.aug.h.cc <- results.h

ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
atc.mult <- ZZ ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
out.form.c <- Y.c ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2)
out.form.h <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2) + I(X1*X3)

for(i in 1:M) {
  data <- weight.data.reps %>% filter(Ite == i)
  
  # Constant treatment effect
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.cc$IPW.c[i] <- PE
  sim.aug.c.cc$IPW.c.var[i] <- Var
  sim.aug.c.cc$IPW.c.lwr[i] <- Lower
  sim.aug.c.cc$IPW.c.upr[i] <- Upper
  sim.aug.c.cc$IPW.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- IPW with 0.05 trimming
  trim1 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.05)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = trim1$data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.cc$IPW.5.c[i] <- PE
  sim.aug.c.cc$IPW.5.c.var[i] <- Var
  sim.aug.c.cc$IPW.5.c.lwr[i] <- Lower
  sim.aug.c.cc$IPW.5.c.upr[i] <- Upper
  sim.aug.c.cc$IPW.5.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- IPW with 0.10 trimming
  trim2 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.10)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = trim2$data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.cc$IPW.10.c[i] <- PE
  sim.aug.c.cc$IPW.10.c.var[i] <- Var
  sim.aug.c.cc$IPW.10.c.lwr[i] <- Lower
  sim.aug.c.cc$IPW.10.c.upr[i] <- Upper
  sim.aug.c.cc$IPW.10.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- IPW with 0.15 trimming
  trim3 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.15)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = trim3$data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.cc$IPW.15.c[i] <- PE
  sim.aug.c.cc$IPW.15.c.var[i] <- Var
  sim.aug.c.cc$IPW.15.c.lwr[i] <- Lower
  sim.aug.c.cc$IPW.15.c.upr[i] <- Upper
  sim.aug.c.cc$IPW.15.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATO
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "overlap", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.cc$ATO.c[i] <- PE
  sim.aug.c.cc$ATO.c.var[i] <- Var
  sim.aug.c.cc$ATO.c.lwr[i] <- Lower
  sim.aug.c.cc$ATO.c.upr[i] <- Upper
  sim.aug.c.cc$ATO.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATM
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "matching", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.cc$ATM.c[i] <- PE
  sim.aug.c.cc$ATM.c.var[i] <- Var
  sim.aug.c.cc$ATM.c.lwr[i] <- Lower
  sim.aug.c.cc$ATM.c.upr[i] <- Upper
  sim.aug.c.cc$ATM.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATEN
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "entropy", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.cc$ATEN.c[i] <- PE
  sim.aug.c.cc$ATEN.c.var[i] <- Var
  sim.aug.c.cc$ATEN.c.lwr[i] <- Lower
  sim.aug.c.cc$ATEN.c.upr[i] <- Upper
  sim.aug.c.cc$ATEN.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATC
  ps.cov <- as.matrix( data[,c("X1","X2","X3","X4","X5","X6","X7")])
  out.cov <- as.matrix(cbind(data[,c("X1","X2","X3","X4")],
                             X1X2 = data[,"X1"]*data[,"X2"], X1sq = data[,"X1"]^2, X2sq = data[,"X2"]^2))
  
  result <- ATC(y = data$Y.c, z = data$Z, X = ps.cov, DR = TRUE, X.out = out.cov)
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.c.cc$ATC.c[i] <- PE
  sim.aug.c.cc$ATC.c.var[i] <- Var
  sim.aug.c.cc$ATC.c.lwr[i] <- Lower
  sim.aug.c.cc$ATC.c.upr[i] <- Upper
  sim.aug.c.cc$ATC.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- ATT(y = data$Y.c, z = data$Z, X = ps.cov, DR = TRUE, X.out = out.cov)
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.c.cc$ATT.c[i] <- PE
  sim.aug.c.cc$ATT.c.var[i] <- Var
  sim.aug.c.cc$ATT.c.lwr[i] <- Lower
  sim.aug.c.cc$ATT.c.upr[i] <- Upper
  sim.aug.c.cc$ATT.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # Heterogeneous treatment effect
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.cc$IPW.h[i] <- PE
  sim.aug.h.cc$IPW.h.var[i] <- Var
  sim.aug.h.cc$IPW.h.lwr[i] <- Lower
  sim.aug.h.cc$IPW.h.upr[i] <- Upper
  sim.aug.h.cc$IPW.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.h < Upper & PE.h$True.IPW.h > Lower), 1, 0))
  
  # --- IPW with 0.05 trimming
  trim1 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.05)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = trim1$data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.cc$IPW.5.h[i] <- PE
  sim.aug.h.cc$IPW.5.h.var[i] <- Var
  sim.aug.h.cc$IPW.5.h.lwr[i] <- Lower
  sim.aug.h.cc$IPW.5.h.upr[i] <- Upper
  sim.aug.h.cc$IPW.5.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.5.h < Upper & PE.h$True.IPW.5.h > Lower), 1, 0))
  
  # --- IPW with 0.10 trimming
  trim2 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.10)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = trim2$data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.cc$IPW.10.h[i] <- PE
  sim.aug.h.cc$IPW.10.h.var[i] <- Var
  sim.aug.h.cc$IPW.10.h.lwr[i] <- Lower
  sim.aug.h.cc$IPW.10.h.upr[i] <- Upper
  sim.aug.h.cc$IPW.10.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.10.h < Upper & PE.h$True.IPW.10.h > Lower), 1, 0))
  
  # --- IPW with 0.15 trimming
  trim3 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.15)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = trim3$data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.cc$IPW.15.h[i] <- PE
  sim.aug.h.cc$IPW.15.h.var[i] <- Var
  sim.aug.h.cc$IPW.15.h.lwr[i] <- Lower
  sim.aug.h.cc$IPW.15.h.upr[i] <- Upper
  sim.aug.h.cc$IPW.15.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.15.h < Upper & PE.h$True.IPW.15.h > Lower), 1, 0))
  
  # --- ATO
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "overlap", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.cc$ATO.h[i] <- PE
  sim.aug.h.cc$ATO.h.var[i] <- Var
  sim.aug.h.cc$ATO.h.lwr[i] <- Lower
  sim.aug.h.cc$ATO.h.upr[i] <- Upper
  sim.aug.h.cc$ATO.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATO.h < Upper & PE.h$True.ATO.h > Lower), 1, 0))
  
  # --- ATM
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "matching", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.cc$ATM.h[i] <- PE
  sim.aug.h.cc$ATM.h.var[i] <- Var
  sim.aug.h.cc$ATM.h.lwr[i] <- Lower
  sim.aug.h.cc$ATM.h.upr[i] <- Upper
  sim.aug.h.cc$ATM.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATM.h < Upper & PE.h$True.ATM.h > Lower), 1, 0))
  
  # --- ATEN
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "entropy", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.cc$ATEN.h[i] <- PE
  sim.aug.h.cc$ATEN.h.var[i] <- Var
  sim.aug.h.cc$ATEN.h.lwr[i] <- Lower
  sim.aug.h.cc$ATEN.h.upr[i] <- Upper
  sim.aug.h.cc$ATEN.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATEN.h < Upper & PE.h$True.ATEN.h > Lower), 1, 0))
  
  # --- ATC
  ps.cov <- as.matrix( data[,c("X1","X2","X3","X4","X5","X6","X7")])
  out.cov <- as.matrix(cbind(data[,c("X1","X2","X3","X4")],
                             X1X2 = data[,"X1"]*data[,"X2"], X1sq = data[,"X1"]^2,
                             X2sq = data[,"X2"]^2, X1X3 = data[,"X1"]*data[,"X3"]))
  
  result <- ATC(y = data$Y.h, z = data$Z, X = ps.cov, DR = TRUE, X.out = out.cov)
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.h.cc$ATC.h[i] <- PE
  sim.aug.h.cc$ATC.h.var[i] <- Var
  sim.aug.h.cc$ATC.h.lwr[i] <- Lower
  sim.aug.h.cc$ATC.h.upr[i] <- Upper
  sim.aug.h.cc$ATC.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))
  
  # --- ATT
  result <- ATT(y = data$Y.h, z = data$Z, X = ps.cov, DR = TRUE, X.out = out.cov)
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.h.cc$ATT.h[i] <- PE
  sim.aug.h.cc$ATT.h.var[i] <- Var
  sim.aug.h.cc$ATT.h.lwr[i] <- Lower
  sim.aug.h.cc$ATT.h.upr[i] <- Upper
  sim.aug.h.cc$ATT.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
}

# Case 2: PS model is correctly specified but OR model is misspecified
sim.aug.c.cm <- results.c
sim.aug.h.cm <- results.h

ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
atc.mult <- ZZ ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
out.form.c <- Y.c ~ X1 + X2 + X3 + X4
out.form.h <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X3)

for(i in 1:M) {
  
  data <- weight.data.reps %>% filter(Ite == i)
  
  # Constant treatment effect
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.cm$IPW.c[i] <- PE
  sim.aug.c.cm$IPW.c.var[i] <- Var
  sim.aug.c.cm$IPW.c.lwr[i] <- Lower
  sim.aug.c.cm$IPW.c.upr[i] <- Upper
  sim.aug.c.cm$IPW.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- IPW with 0.05 trimming
  trim1 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.05)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = trim1$data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.cm$IPW.5.c[i] <- PE
  sim.aug.c.cm$IPW.5.c.var[i] <- Var
  sim.aug.c.cm$IPW.5.c.lwr[i] <- Lower
  sim.aug.c.cm$IPW.5.c.upr[i] <- Upper
  sim.aug.c.cm$IPW.5.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- IPW with 0.10 trimming
  trim2 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.10)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = trim2$data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.cm$IPW.10.c[i] <- PE
  sim.aug.c.cm$IPW.10.c.var[i] <- Var
  sim.aug.c.cm$IPW.10.c.lwr[i] <- Lower
  sim.aug.c.cm$IPW.10.c.upr[i] <- Upper
  sim.aug.c.cm$IPW.10.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- IPW with 0.15 trimming
  trim3 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.15)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = trim3$data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.cm$IPW.15.c[i] <- PE
  sim.aug.c.cm$IPW.15.c.var[i] <- Var
  sim.aug.c.cm$IPW.15.c.lwr[i] <- Lower
  sim.aug.c.cm$IPW.15.c.upr[i] <- Upper
  sim.aug.c.cm$IPW.15.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATO
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "overlap", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.cm$ATO.c[i] <- PE
  sim.aug.c.cm$ATO.c.var[i] <- Var
  sim.aug.c.cm$ATO.c.lwr[i] <- Lower
  sim.aug.c.cm$ATO.c.upr[i] <- Upper
  sim.aug.c.cm$ATO.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATM
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "matching", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.cm$ATM.c[i] <- PE
  sim.aug.c.cm$ATM.c.var[i] <- Var
  sim.aug.c.cm$ATM.c.lwr[i] <- Lower
  sim.aug.c.cm$ATM.c.upr[i] <- Upper
  sim.aug.c.cm$ATM.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATEN
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "entropy", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.cm$ATEN.c[i] <- PE
  sim.aug.c.cm$ATEN.c.var[i] <- Var
  sim.aug.c.cm$ATEN.c.lwr[i] <- Lower
  sim.aug.c.cm$ATEN.c.upr[i] <- Upper
  sim.aug.c.cm$ATEN.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATC
  ps.cov <- as.matrix( data[,c("X1","X2","X3","X4","X5","X6","X7")])
  out.cov <- as.matrix(cbind(data[,c("X1","X2","X3","X4")]))
  
  result <- ATC(y = data$Y.c, z = data$Z, X = ps.cov, DR = TRUE, X.out = out.cov)
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.c.cm$ATC.c[i] <- PE
  sim.aug.c.cm$ATC.c.var[i] <- Var
  sim.aug.c.cm$ATC.c.lwr[i] <- Lower
  sim.aug.c.cm$ATC.c.upr[i] <- Upper
  sim.aug.c.cm$ATC.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- ATT(y = data$Y.c, z = data$Z, X = ps.cov, DR = TRUE, X.out = out.cov)
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.c.cm$ATT.c[i] <- PE
  sim.aug.c.cm$ATT.c.var[i] <- Var
  sim.aug.c.cm$ATT.c.lwr[i] <- Lower
  sim.aug.c.cm$ATT.c.upr[i] <- Upper
  sim.aug.c.cm$ATT.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # Heterogeneous treatment effect
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.cm$IPW.h[i] <- PE
  sim.aug.h.cm$IPW.h.var[i] <- Var
  sim.aug.h.cm$IPW.h.lwr[i] <- Lower
  sim.aug.h.cm$IPW.h.upr[i] <- Upper
  sim.aug.h.cm$IPW.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.h < Upper & PE.h$True.IPW.h > Lower), 1, 0))
  
  # --- IPW with 0.05 trimming
  trim1 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.05)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = trim1$data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.cm$IPW.5.h[i] <- PE
  sim.aug.h.cm$IPW.5.h.var[i] <- Var
  sim.aug.h.cm$IPW.5.h.lwr[i] <- Lower
  sim.aug.h.cm$IPW.5.h.upr[i] <- Upper
  sim.aug.h.cm$IPW.5.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.5.h < Upper & PE.h$True.IPW.5.h > Lower), 1, 0))
  
  # --- IPW with 0.10 trimming
  trim2 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.10)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = trim2$data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.cm$IPW.10.h[i] <- PE
  sim.aug.h.cm$IPW.10.h.var[i] <- Var
  sim.aug.h.cm$IPW.10.h.lwr[i] <- Lower
  sim.aug.h.cm$IPW.10.h.upr[i] <- Upper
  sim.aug.h.cm$IPW.10.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.10.h < Upper & PE.h$True.IPW.10.h > Lower), 1, 0))
  
  # --- IPW with 0.15 trimming
  trim3 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.15)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = trim3$data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.cm$IPW.15.h[i] <- PE
  sim.aug.h.cm$IPW.15.h.var[i] <- Var
  sim.aug.h.cm$IPW.15.h.lwr[i] <- Lower
  sim.aug.h.cm$IPW.15.h.upr[i] <- Upper
  sim.aug.h.cm$IPW.15.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.15.h < Upper & PE.h$True.IPW.15.h > Lower), 1, 0))
  
  # --- ATO
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "overlap", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.cm$ATO.h[i] <- PE
  sim.aug.h.cm$ATO.h.var[i] <- Var
  sim.aug.h.cm$ATO.h.lwr[i] <- Lower
  sim.aug.h.cm$ATO.h.upr[i] <- Upper
  sim.aug.h.cm$ATO.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATO.h < Upper & PE.h$True.ATO.h > Lower), 1, 0))
  
  # --- ATM
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "matching", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.cm$ATM.h[i] <- PE
  sim.aug.h.cm$ATM.h.var[i] <- Var
  sim.aug.h.cm$ATM.h.lwr[i] <- Lower
  sim.aug.h.cm$ATM.h.upr[i] <- Upper
  sim.aug.h.cm$ATM.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATM.h < Upper & PE.h$True.ATM.h > Lower), 1, 0))
  
  # --- ATEN
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "entropy", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.cm$ATEN.h[i] <- PE
  sim.aug.h.cm$ATEN.h.var[i] <- Var
  sim.aug.h.cm$ATEN.h.lwr[i] <- Lower
  sim.aug.h.cm$ATEN.h.upr[i] <- Upper
  sim.aug.h.cm$ATEN.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATEN.h < Upper & PE.h$True.ATEN.h > Lower), 1, 0))
  
  # --- ATC
  ps.cov <- as.matrix( data[,c("X1","X2","X3","X4","X5","X6","X7")])
  out.cov <- as.matrix(cbind(data[,c("X1","X2","X3","X4")], X1X3 = data[,"X1"]*data[,"X3"]))
  
  result <- ATC(y = data$Y.h, z = data$Z, X = ps.cov, DR = TRUE, X.out = out.cov)
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.h.cm$ATC.h[i] <- PE
  sim.aug.h.cm$ATC.h.var[i] <- Var
  sim.aug.h.cm$ATC.h.lwr[i] <- Lower
  sim.aug.h.cm$ATC.h.upr[i] <- Upper
  sim.aug.h.cm$ATC.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))
  
  # --- ATT
  result <- ATT(y = data$Y.h, z = data$Z, X = ps.cov, DR = TRUE, X.out = out.cov)
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.h.cm$ATT.h[i] <- PE
  sim.aug.h.cm$ATT.h.var[i] <- Var
  sim.aug.h.cm$ATT.h.lwr[i] <- Lower
  sim.aug.h.cm$ATT.h.upr[i] <- Upper
  sim.aug.h.cm$ATT.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
}

# Case 3: PS model is misspecified but OR model is correctly specified
sim.aug.c.mc <- results.c
sim.aug.h.mc <- results.h

ps.mult <- Z ~ X1 + X2 + X3 + X4
atc.mult <- ZZ ~ X1 + X2 + X3 + X4
out.form.c <- Y.c ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2)
out.form.h <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2) + I(X1*X3)

for(i in 1:M) {
  
  data <- weight.data.reps %>% filter(Ite == i)
  
  # Constant treatment effect
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.mc$IPW.c[i] <- PE
  sim.aug.c.mc$IPW.c.var[i] <- Var
  sim.aug.c.mc$IPW.c.lwr[i] <- Lower
  sim.aug.c.mc$IPW.c.upr[i] <- Upper
  sim.aug.c.mc$IPW.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- IPW with 0.05 trimming
  trim1 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.05)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = trim1$data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.mc$IPW.5.c[i] <- PE
  sim.aug.c.mc$IPW.5.c.var[i] <- Var
  sim.aug.c.mc$IPW.5.c.lwr[i] <- Lower
  sim.aug.c.mc$IPW.5.c.upr[i] <- Upper
  sim.aug.c.mc$IPW.5.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- IPW with 0.10 trimming
  trim2 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.10)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = trim2$data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.mc$IPW.10.c[i] <- PE
  sim.aug.c.mc$IPW.10.c.var[i] <- Var
  sim.aug.c.mc$IPW.10.c.lwr[i] <- Lower
  sim.aug.c.mc$IPW.10.c.upr[i] <- Upper
  sim.aug.c.mc$IPW.10.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- IPW with 0.15 trimming
  trim3 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.15)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = trim3$data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.mc$IPW.15.c[i] <- PE
  sim.aug.c.mc$IPW.15.c.var[i] <- Var
  sim.aug.c.mc$IPW.15.c.lwr[i] <- Lower
  sim.aug.c.mc$IPW.15.c.upr[i] <- Upper
  sim.aug.c.mc$IPW.15.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATO
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "overlap", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.mc$ATO.c[i] <- PE
  sim.aug.c.mc$ATO.c.var[i] <- Var
  sim.aug.c.mc$ATO.c.lwr[i] <- Lower
  sim.aug.c.mc$ATO.c.upr[i] <- Upper
  sim.aug.c.mc$ATO.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATM
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "matching", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.mc$ATM.c[i] <- PE
  sim.aug.c.mc$ATM.c.var[i] <- Var
  sim.aug.c.mc$ATM.c.lwr[i] <- Lower
  sim.aug.c.mc$ATM.c.upr[i] <- Upper
  sim.aug.c.mc$ATM.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATEN
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "entropy", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.mc$ATEN.c[i] <- PE
  sim.aug.c.mc$ATEN.c.var[i] <- Var
  sim.aug.c.mc$ATEN.c.lwr[i] <- Lower
  sim.aug.c.mc$ATEN.c.upr[i] <- Upper
  sim.aug.c.mc$ATEN.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATC
  ps.cov <- as.matrix( data[,c("X1","X2","X3","X4")])
  out.cov <- as.matrix(cbind(data[,c("X1","X2","X3","X4")],
                             X1X2 = data[,"X1"]*data[,"X2"], X1sq = data[,"X1"]^2,X2sq = data[,"X2"]^2))
  
  result <- ATC(y = data$Y.c, z = data$Z, X = ps.cov, DR = TRUE, X.out = out.cov)
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.c.mc$ATC.c[i] <- PE
  sim.aug.c.mc$ATC.c.var[i] <- Var
  sim.aug.c.mc$ATC.c.lwr[i] <- Lower
  sim.aug.c.mc$ATC.c.upr[i] <- Upper
  sim.aug.c.mc$ATC.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- ATT(y = data$Y.c, z = data$Z, X = ps.cov, DR = TRUE, X.out = out.cov)
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.c.mc$ATT.c[i] <- PE
  sim.aug.c.mc$ATT.c.var[i] <- Var
  sim.aug.c.mc$ATT.c.lwr[i] <- Lower
  sim.aug.c.mc$ATT.c.upr[i] <- Upper
  sim.aug.c.mc$ATT.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # Heterogeneous treatment effect
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.mc$IPW.h[i] <- PE
  sim.aug.h.mc$IPW.h.var[i] <- Var
  sim.aug.h.mc$IPW.h.lwr[i] <- Lower
  sim.aug.h.mc$IPW.h.upr[i] <- Upper
  sim.aug.h.mc$IPW.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.h < Upper & PE.h$True.IPW.h > Lower), 1, 0))
  
  # --- IPW with 0.05 trimming
  trim1 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.05)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = trim1$data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.mc$IPW.5.h[i] <- PE
  sim.aug.h.mc$IPW.5.h.var[i] <- Var
  sim.aug.h.mc$IPW.5.h.lwr[i] <- Lower
  sim.aug.h.mc$IPW.5.h.upr[i] <- Upper
  sim.aug.h.mc$IPW.5.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.5.h < Upper & PE.h$True.IPW.5.h > Lower), 1, 0))
  
  # --- IPW with 0.10 trimming
  trim2 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.10)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = trim2$data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.mc$IPW.10.h[i] <- PE
  sim.aug.h.mc$IPW.10.h.var[i] <- Var
  sim.aug.h.mc$IPW.10.h.lwr[i] <- Lower
  sim.aug.h.mc$IPW.10.h.upr[i] <- Upper
  sim.aug.h.mc$IPW.10.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.10.h < Upper & PE.h$True.IPW.10.h > Lower), 1, 0))
  
  # --- IPW with 0.15 trimming
  trim3 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.15)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = trim3$data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.mc$IPW.15.h[i] <- PE
  sim.aug.h.mc$IPW.15.h.var[i] <- Var
  sim.aug.h.mc$IPW.15.h.lwr[i] <- Lower
  sim.aug.h.mc$IPW.15.h.upr[i] <- Upper
  sim.aug.h.mc$IPW.15.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.15.h < Upper & PE.h$True.IPW.15.h > Lower), 1, 0))
  
  # --- ATO
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "overlap", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.mc$ATO.h[i] <- PE
  sim.aug.h.mc$ATO.h.var[i] <- Var
  sim.aug.h.mc$ATO.h.lwr[i] <- Lower
  sim.aug.h.mc$ATO.h.upr[i] <- Upper
  sim.aug.h.mc$ATO.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATO.h < Upper & PE.h$True.ATO.h > Lower), 1, 0))
  
  # --- ATM
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "matching", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.mc$ATM.h[i] <- PE
  sim.aug.h.mc$ATM.h.var[i] <- Var
  sim.aug.h.mc$ATM.h.lwr[i] <- Lower
  sim.aug.h.mc$ATM.h.upr[i] <- Upper
  sim.aug.h.mc$ATM.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATM.h < Upper & PE.h$True.ATM.h > Lower), 1, 0))
  
  # --- ATEN
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "entropy", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.mc$ATEN.h[i] <- PE
  sim.aug.h.mc$ATEN.h.var[i] <- Var
  sim.aug.h.mc$ATEN.h.lwr[i] <- Lower
  sim.aug.h.mc$ATEN.h.upr[i] <- Upper
  sim.aug.h.mc$ATEN.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATEN.h < Upper & PE.h$True.ATEN.h > Lower), 1, 0))
  
  # --- ATC
  ps.cov <- as.matrix( data[,c("X1","X2","X3","X4")])
  out.cov <- as.matrix(cbind(data[,c("X1","X2","X3","X4")],
                             X1X2 = data[,"X1"]*data[,"X2"], X1sq = data[,"X1"]^2, 
                             X2sq = data[,"X2"]^2, X1X3 = data[,"X1"]*data[,"X3"]))
  
  result <- ATC(y = data$Y.h, z = data$Z, X = ps.cov, DR = TRUE, X.out = out.cov)
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.h.mc$ATC.h[i] <- PE
  sim.aug.h.mc$ATC.h.var[i] <- Var
  sim.aug.h.mc$ATC.h.lwr[i] <- Lower
  sim.aug.h.mc$ATC.h.upr[i] <- Upper
  sim.aug.h.mc$ATC.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))
  
  # --- ATT
  result <- ATT(y = data$Y.h, z = data$Z, X = ps.cov, DR = TRUE, X.out = out.cov)
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.h.mc$ATT.h[i] <- PE
  sim.aug.h.mc$ATT.h.var[i] <- Var
  sim.aug.h.mc$ATT.h.lwr[i] <- Lower
  sim.aug.h.mc$ATT.h.upr[i] <- Upper
  sim.aug.h.mc$ATT.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
  
}

# Case 4: Both PS and OR models are misspecified

sim.aug.c.mm <- results.c
sim.aug.h.mm <- results.h

ps.mult <- Z ~ X1 + X2 + X3 + X4
atc.mult <- ZZ ~ X1 + X2 + X3 + X4
out.form.c <- Y.c ~ X1 + X2 + X3 + X4 
out.form.h <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X3)

for(i in 1:M) {
  data <- weight.data.reps %>% filter(Ite == i)
  
  # Constant treatment effect
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.mm$IPW.c[i] <- PE
  sim.aug.c.mm$IPW.c.var[i] <- Var
  sim.aug.c.mm$IPW.c.lwr[i] <- Lower
  sim.aug.c.mm$IPW.c.upr[i] <- Upper
  sim.aug.c.mm$IPW.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- IPW with 0.05 trimming
  trim1 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.05)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = trim1$data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.mm$IPW.5.c[i] <- PE
  sim.aug.c.mm$IPW.5.c.var[i] <- Var
  sim.aug.c.mm$IPW.5.c.lwr[i] <- Lower
  sim.aug.c.mm$IPW.5.c.upr[i] <- Upper
  sim.aug.c.mm$IPW.5.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- IPW with 0.10 trimming
  trim2 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.10)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = trim2$data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.mm$IPW.10.c[i] <- PE
  sim.aug.c.mm$IPW.10.c.var[i] <- Var
  sim.aug.c.mm$IPW.10.c.lwr[i] <- Lower
  sim.aug.c.mm$IPW.10.c.upr[i] <- Upper
  sim.aug.c.mm$IPW.10.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- IPW with 0.15 trimming
  trim3 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.15)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = trim3$data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.mm$IPW.15.c[i] <- PE
  sim.aug.c.mm$IPW.15.c.var[i] <- Var
  sim.aug.c.mm$IPW.15.c.lwr[i] <- Lower
  sim.aug.c.mm$IPW.15.c.upr[i] <- Upper
  sim.aug.c.mm$IPW.15.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATO
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "overlap", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.mm$ATO.c[i] <- PE
  sim.aug.c.mm$ATO.c.var[i] <- Var
  sim.aug.c.mm$ATO.c.lwr[i] <- Lower
  sim.aug.c.mm$ATO.c.upr[i] <- Upper
  sim.aug.c.mm$ATO.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATM
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "matching", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.mm$ATM.c[i] <- PE
  sim.aug.c.mm$ATM.c.var[i] <- Var
  sim.aug.c.mm$ATM.c.lwr[i] <- Lower
  sim.aug.c.mm$ATM.c.upr[i] <- Upper
  sim.aug.c.mm$ATM.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATEN
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "entropy", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.mm$ATEN.c[i] <- PE
  sim.aug.c.mm$ATEN.c.var[i] <- Var
  sim.aug.c.mm$ATEN.c.lwr[i] <- Lower
  sim.aug.c.mm$ATEN.c.upr[i] <- Upper
  sim.aug.c.mm$ATEN.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATC
  ps.cov <- as.matrix( data[,c("X1","X2","X3","X4")])
  out.cov <- as.matrix(cbind(data[,c("X1","X2","X3","X4")]))
  
  result <- ATC(y = data$Y.c, z = data$Z, X = ps.cov, DR = TRUE, X.out = out.cov)
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.c.mm$ATC.c[i] <- PE
  sim.aug.c.mm$ATC.c.var[i] <- Var
  sim.aug.c.mm$ATC.c.lwr[i] <- Lower
  sim.aug.c.mm$ATC.c.upr[i] <- Upper
  sim.aug.c.mm$ATC.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- ATT(y = data$Y.c, z = data$Z, X = ps.cov, DR = TRUE, X.out = out.cov)
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.c.mm$ATT.c[i] <- PE
  sim.aug.c.mm$ATT.c.var[i] <- Var
  sim.aug.c.mm$ATT.c.lwr[i] <- Lower
  sim.aug.c.mm$ATT.c.upr[i] <- Upper
  sim.aug.c.mm$ATT.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # Heterogeneous treatment effect
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.mm$IPW.h[i] <- PE
  sim.aug.h.mm$IPW.h.var[i] <- Var
  sim.aug.h.mm$IPW.h.lwr[i] <- Lower
  sim.aug.h.mm$IPW.h.upr[i] <- Upper
  sim.aug.h.mm$IPW.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.h < Upper & PE.h$True.IPW.h > Lower), 1, 0))
  
  # --- IPW with 0.05 trimming
  trim1 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.05)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = trim1$data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.mm$IPW.5.h[i] <- PE
  sim.aug.h.mm$IPW.5.h.var[i] <- Var
  sim.aug.h.mm$IPW.5.h.lwr[i] <- Lower
  sim.aug.h.mm$IPW.5.h.upr[i] <- Upper
  sim.aug.h.mm$IPW.5.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.5.h < Upper & PE.h$True.IPW.5.h > Lower), 1, 0))
  
  # --- IPW with 0.10 trimming
  trim2 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.10)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = trim2$data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.mm$IPW.10.h[i] <- PE
  sim.aug.h.mm$IPW.10.h.var[i] <- Var
  sim.aug.h.mm$IPW.10.h.lwr[i] <- Lower
  sim.aug.h.mm$IPW.10.h.upr[i] <- Upper
  sim.aug.h.mm$IPW.10.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.10.h < Upper & PE.h$True.IPW.10.h > Lower), 1, 0))
  
  # --- IPW with 0.15 trimming
  trim3 <- PStrim(data = data, ps.formula = ps.mult, delta = 0.15)
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = trim3$data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.mm$IPW.15.h[i] <- PE
  sim.aug.h.mm$IPW.15.h.var[i] <- Var
  sim.aug.h.mm$IPW.15.h.lwr[i] <- Lower
  sim.aug.h.mm$IPW.15.h.upr[i] <- Upper
  sim.aug.h.mm$IPW.15.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.15.h < Upper & PE.h$True.IPW.15.h > Lower), 1, 0))
  
  # --- ATO
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "overlap", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.mm$ATO.h[i] <- PE
  sim.aug.h.mm$ATO.h.var[i] <- Var
  sim.aug.h.mm$ATO.h.lwr[i] <- Lower
  sim.aug.h.mm$ATO.h.upr[i] <- Upper
  sim.aug.h.mm$ATO.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATO.h < Upper & PE.h$True.ATO.h > Lower), 1, 0))
  
  # --- ATM
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "matching", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.mm$ATM.h[i] <- PE
  sim.aug.h.mm$ATM.h.var[i] <- Var
  sim.aug.h.mm$ATM.h.lwr[i] <- Lower
  sim.aug.h.mm$ATM.h.upr[i] <- Upper
  sim.aug.h.mm$ATM.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATM.h < Upper & PE.h$True.ATM.h > Lower), 1, 0))
  
  # --- ATEN
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "entropy", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.mm$ATEN.h[i] <- PE
  sim.aug.h.mm$ATEN.h.var[i] <- Var
  sim.aug.h.mm$ATEN.h.lwr[i] <- Lower
  sim.aug.h.mm$ATEN.h.upr[i] <- Upper
  sim.aug.h.mm$ATEN.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATEN.h < Upper & PE.h$True.ATEN.h > Lower), 1, 0))
  
  # --- ATC
  ps.cov <- as.matrix( data[,c("X1","X2","X3","X4")])
  out.cov <- as.matrix(cbind(data[,c("X1","X2","X3","X4")], X1X3 = data[,"X1"]*data[,"X3"]))
  
  result <- ATC(y = data$Y.h, z = data$Z, X = ps.cov, DR = TRUE, X.out = out.cov)
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.h.mm$ATC.h[i] <- PE
  sim.aug.h.mm$ATC.h.var[i] <- Var
  sim.aug.h.mm$ATC.h.lwr[i] <- Lower
  sim.aug.h.mm$ATC.h.upr[i] <- Upper
  sim.aug.h.mm$ATC.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))
  
  # --- ATT
  result <- ATT(y = data$Y.h, z = data$Z, X = ps.cov, DR = TRUE, X.out = out.cov)
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.h.mm$ATT.h[i] <- PE
  sim.aug.h.mm$ATT.h.var[i] <- Var
  sim.aug.h.mm$ATT.h.lwr[i] <- Lower
  sim.aug.h.mm$ATT.h.upr[i] <- Upper
  sim.aug.h.mm$ATT.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
}

# One way to address NA problem in variance and CIs
num.na <- function(x) { sum(is.na(x)) }
all.sim.dat <- cbind(sim.wt.c, sim.wt.h, sim.aug.c.cc, sim.aug.h.cc, 
                     sim.aug.c.cm, sim.aug.h.cm, sim.aug.c.mc, sim.aug.h.mc, 
                     sim.aug.c.mm, sim.aug.h.mm) 

colnames(all.sim.dat) = as.character(1:460)
all.sim.dat <- all.sim.dat %>% mutate(Num.NA = rep(0, M), Row.Num = 1:M)
all.sim.dat$Num.NA <- apply(all.sim.dat, 1, num.na)

all.sim.dat <- all.sim.dat %>% filter(Num.NA==0)

R <- all.sim.dat$Row.Num
if(length(all.sim.dat$Row.Num) < M0) { warning("The number of rows < M0!") } else {
  R <- sample(all.sim.dat$Row.Num, size = M0, replace = FALSE)
}

sim1.wt.c <- sim.wt.c %>% filter(Row.Num %in% R)
sim1.wt.h <- sim.wt.h %>% filter(Row.Num %in% R)

sim1.aug.c.cc <- sim.aug.c.cc %>% filter(Row.Num %in% R)
sim1.aug.h.cc <- sim.aug.h.cc %>% filter(Row.Num %in% R)

sim1.aug.c.cm <- sim.aug.c.cm %>% filter(Row.Num %in% R)
sim1.aug.h.cm <- sim.aug.h.cm %>% filter(Row.Num %in% R)

sim1.aug.c.mc <- sim.aug.c.mc %>% filter(Row.Num %in% R)
sim1.aug.h.mc <- sim.aug.h.mc %>% filter(Row.Num %in% R)

sim1.aug.c.mm <- sim.aug.c.mm %>% filter(Row.Num %in% R)
sim1.aug.h.mm <- sim.aug.h.mm %>% filter(Row.Num %in% R)

weight.data.reps1 <- weight.data.reps %>% filter(Ite %in% R)

# --- proportion
sim1.p <- length(weight.data.reps1$Z[weight.data.reps1$Z == 1])/(M0*N)

# --- ratio
ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7

ps.var.ctrl <- rep(NA, M0)
ps.var.trt <- rep(NA, M0)

k = 1
for(i in R) {
  ps.dat <- data.frame(ps = rep(NA, N), Z = rep(NA, N))
  data <- weight.data.reps1 %>% filter(Ite == i)
  
  ps <- SumStat(ps.formula = ps.mult, data = data, weight = "overlap")$propensity
  ps.dat$ps <- ps[, "1"]
  ps.dat$Z <- data$Z
  
  ps.var.ctrl[k] = ps.dat %>% filter(Z==0) %>% summarize(var(ps)) %>% as.numeric()
  ps.var.trt[k] = ps.dat %>% filter(Z==1) %>% summarize(var(ps)) %>% as.numeric()
  k = k+1
}

sim1.ratio <- mean(ps.var.trt/ps.var.ctrl)
md1.true.h <- PE.h

# Save the data
save(file = "md1_sims.RData", weight.data.reps1, sim1.p, sim1.ratio, md1.true.h,
     sim1.wt.c, sim1.wt.h, sim1.aug.c.cc, sim1.aug.h.cc, 
     sim1.aug.c.cm, sim1.aug.h.cm, sim1.aug.c.mc, sim1.aug.h.mc, sim1.aug.c.mm, sim1.aug.h.mm )

