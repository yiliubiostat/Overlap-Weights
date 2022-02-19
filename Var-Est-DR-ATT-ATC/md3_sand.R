### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ ATT/ATC DR Variance       ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

####### Sandwich variance (new), Model 3

### Created by Yi Liu, modified by Yunji Zhou
### Oct 9, 2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~ Using M Replicates for Estimands ~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(dplyr)
library(PSweight)
source("newSand_func.R")

load("md3_sims_PSwt.RData")
weight.data.reps <- weight.data.reps3
rm(sim3.aug.c.cc,sim3.aug.c.cm,sim3.aug.c.mc,sim3.aug.c.mm,
   sim3.aug.h.cc,sim3.aug.h.cm,sim3.aug.h.mc,sim3.aug.h.mm)

iteration <- unique(weight.data.reps$Ite)
M <- 2000

results.new <- data.frame(IPW.new = rep(NA, M),
                          IPW.new.var = rep(NA, M),
                          IPW.new.lwr = rep(NA, M),
                          IPW.new.upr = rep(NA, M),
                          IPW.new.ifci = rep(NA, M),
                          
                          ATC.new = rep(NA, M),
                          ATC.new.var = rep(NA, M),
                          ATC.new.lwr = rep(NA, M),
                          ATC.new.upr = rep(NA, M),
                          ATC.new.ifci = rep(NA, M),
                          
                          ATT.new = rep(NA, M),
                          ATT.new.var = rep(NA, M),
                          ATT.new.lwr = rep(NA, M),
                          ATT.new.upr = rep(NA, M),
                          ATT.new.ifci = rep(NA, M),
                          
                          Row.Num = 1:M)

results.old <- data.frame(IPW.old = rep(NA, M),
                          IPW.old.var = rep(NA, M),
                          IPW.old.lwr = rep(NA, M),
                          IPW.old.upr = rep(NA, M),
                          IPW.old.ifci = rep(NA, M),
                          
                          ATC.old = rep(NA, M),
                          ATC.old.var = rep(NA, M),
                          ATC.old.lwr = rep(NA, M),
                          ATC.old.upr = rep(NA, M),
                          ATC.old.ifci = rep(NA, M),
                          
                          ATT.old = rep(NA, M),
                          ATT.old.var = rep(NA, M),
                          ATT.old.lwr = rep(NA, M),
                          ATT.old.upr = rep(NA, M),
                          ATT.old.ifci = rep(NA, M),
                          
                          Row.Num = 1:M)
PE.c <- 4
PE.h <- md3.true.h

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~ Augmented Estimators ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Case 1: Both PS and outcome are correctly specified
sim.aug.new.cc <- results.new
sim.aug.old.cc <- results.old

c.sim.aug.new.cc <- results.new
c.sim.aug.old.cc <- results.old

ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
atc.mult <- ZZ ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
out.form.c <- Y.c ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2)
out.form.h <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2) + I(X1*X3)

for(i in 1:M) {
  
  data <- weight.data.reps %>% filter(Ite == iteration[i])
  ps.cov <- as.matrix( data[,c("X1","X2","X3","X4","X5","X6","X7")])
  out.cov <- as.matrix(cbind(data[,c("X1","X2","X3","X4")],
                             X1X2 = data[,"X1"]*data[,"X2"],X1sq = data[,"X1"]^2,
                             X2sq = data[,"X2"]^2,X1X3 = data[,"X1"]*data[,"X3"]))
  
  # Heterogeneous treatment effect
  # Estimation using new sandwich estimator 
  # --- ATE: IPW without trimming
  result <- ATE(y=data$Y.h, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.new.cc$IPW.new[i] <- PE
  sim.aug.new.cc$IPW.new.var[i] <- Var
  sim.aug.new.cc$IPW.new.lwr[i] <- Lower
  sim.aug.new.cc$IPW.new.upr[i] <- Upper
  sim.aug.new.cc$IPW.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.h < Upper & PE.h$True.IPW.h > Lower), 1, 0))
  
  # --- ATC
  result <- ATC(y=data$Y.h, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.new.cc$ATC.new[i] <- PE
  sim.aug.new.cc$ATC.new.var[i] <- Var
  sim.aug.new.cc$ATC.new.lwr[i] <- Lower
  sim.aug.new.cc$ATC.new.upr[i] <- Upper
  sim.aug.new.cc$ATC.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))
  
  # --- ATT
  result <- ATT(y=data$Y.h, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.new.cc$ATT.new[i] <- PE
  sim.aug.new.cc$ATT.new.var[i] <- Var
  sim.aug.new.cc$ATT.new.lwr[i] <- Lower
  sim.aug.new.cc$ATT.new.upr[i] <- Upper
  sim.aug.new.cc$ATT.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
  
  # Estimation using old sandwich estimator through PSweight
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.old.cc$IPW.old[i] <- PE
  sim.aug.old.cc$IPW.old.var[i] <- Var
  sim.aug.old.cc$IPW.old.lwr[i] <- Lower
  sim.aug.old.cc$IPW.old.upr[i] <- Upper
  sim.aug.old.cc$IPW.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.h < Upper & PE.h$True.IPW.h > Lower), 1, 0))
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.h", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  sim.aug.old.cc$ATC.old[i] <- PE
  sim.aug.old.cc$ATC.old.var[i] <- Var
  sim.aug.old.cc$ATC.old.lwr[i] <- Lower
  sim.aug.old.cc$ATC.old.upr[i] <- Upper
  sim.aug.old.cc$ATC.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.old.cc$ATT.old[i] <- PE
  sim.aug.old.cc$ATT.old.var[i] <- Var
  sim.aug.old.cc$ATT.old.lwr[i] <- Lower
  sim.aug.old.cc$ATT.old.upr[i] <- Upper
  sim.aug.old.cc$ATT.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
  
  # Homogeneous trt effect
  # Estimation using new sandwich estimator 
  # --- ATE: IPW without trimming
  result <- ATE(y=data$Y.c, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  c.sim.aug.new.cc$IPW.new[i] <- PE
  c.sim.aug.new.cc$IPW.new.var[i] <- Var
  c.sim.aug.new.cc$IPW.new.lwr[i] <- Lower
  c.sim.aug.new.cc$IPW.new.upr[i] <- Upper
  c.sim.aug.new.cc$IPW.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATC
  result <- ATC(y=data$Y.c, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  c.sim.aug.new.cc$ATC.new[i] <- PE
  c.sim.aug.new.cc$ATC.new.var[i] <- Var
  c.sim.aug.new.cc$ATC.new.lwr[i] <- Lower
  c.sim.aug.new.cc$ATC.new.upr[i] <- Upper
  c.sim.aug.new.cc$ATC.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- ATT(y=data$Y.c, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  c.sim.aug.new.cc$ATT.new[i] <- PE
  c.sim.aug.new.cc$ATT.new.var[i] <- Var
  c.sim.aug.new.cc$ATT.new.lwr[i] <- Lower
  c.sim.aug.new.cc$ATT.new.upr[i] <- Upper
  c.sim.aug.new.cc$ATT.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # Estimation using old sandwich estimator through PSweight
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  c.sim.aug.old.cc$IPW.old[i] <- PE
  c.sim.aug.old.cc$IPW.old.var[i] <- Var
  c.sim.aug.old.cc$IPW.old.lwr[i] <- Lower
  c.sim.aug.old.cc$IPW.old.upr[i] <- Upper
  c.sim.aug.old.cc$IPW.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.c", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  c.sim.aug.old.cc$ATC.old[i] <- PE
  c.sim.aug.old.cc$ATC.old.var[i] <- Var
  c.sim.aug.old.cc$ATC.old.lwr[i] <- Lower
  c.sim.aug.old.cc$ATC.old.upr[i] <- Upper
  c.sim.aug.old.cc$ATC.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  c.sim.aug.old.cc$ATT.old[i] <- PE
  c.sim.aug.old.cc$ATT.old.var[i] <- Var
  c.sim.aug.old.cc$ATT.old.lwr[i] <- Lower
  c.sim.aug.old.cc$ATT.old.upr[i] <- Upper
  c.sim.aug.old.cc$ATT.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  print(paste0("CC Ite ", i))
}

# Case 2: PS correct but outcome misspecified
sim.aug.new.cm <- results.new
sim.aug.old.cm <- results.old

c.sim.aug.new.cm <- results.new
c.sim.aug.old.cm <- results.old

ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
atc.mult <- ZZ ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
out.form.c <- Y.c ~ X1 + X2 + X3 + X4
out.form.h <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X3)

for(i in 1:M) {
  
  data <- weight.data.reps %>% filter(Ite == iteration[i])
  ps.cov <- as.matrix( data[,c("X1","X2","X3","X4","X5","X6","X7")])
  out.cov <- as.matrix(cbind(data[,c("X1","X2","X3","X4")],
                             X1X3 = data[,"X1"]*data[,"X3"]))
  
  # Heterogeneous trt effect
  # Estimation using new sandwich estimator 
  # --- ATE: IPW without trimming
  result <- ATE(y=data$Y.h, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.new.cm$IPW.new[i] <- PE
  sim.aug.new.cm$IPW.new.var[i] <- Var
  sim.aug.new.cm$IPW.new.lwr[i] <- Lower
  sim.aug.new.cm$IPW.new.upr[i] <- Upper
  sim.aug.new.cm$IPW.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.h < Upper & PE.h$True.IPW.h > Lower), 1, 0))
  
  # --- ATC
  result <- ATC(y=data$Y.h, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.new.cm$ATC.new[i] <- PE
  sim.aug.new.cm$ATC.new.var[i] <- Var
  sim.aug.new.cm$ATC.new.lwr[i] <- Lower
  sim.aug.new.cm$ATC.new.upr[i] <- Upper
  sim.aug.new.cm$ATC.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))
  
  # --- ATT
  result <- ATT(y=data$Y.h, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.new.cm$ATT.new[i] <- PE
  sim.aug.new.cm$ATT.new.var[i] <- Var
  sim.aug.new.cm$ATT.new.lwr[i] <- Lower
  sim.aug.new.cm$ATT.new.upr[i] <- Upper
  sim.aug.new.cm$ATT.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
  
  # Estimation using old sandwich estimator through PSweight
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.old.cm$IPW.old[i] <- PE
  sim.aug.old.cm$IPW.old.var[i] <- Var
  sim.aug.old.cm$IPW.old.lwr[i] <- Lower
  sim.aug.old.cm$IPW.old.upr[i] <- Upper
  sim.aug.old.cm$IPW.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.h < Upper & PE.h$True.IPW.h > Lower), 1, 0))
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.h", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  sim.aug.old.cm$ATC.old[i] <- PE
  sim.aug.old.cm$ATC.old.var[i] <- Var
  sim.aug.old.cm$ATC.old.lwr[i] <- Lower
  sim.aug.old.cm$ATC.old.upr[i] <- Upper
  sim.aug.old.cm$ATC.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.old.cm$ATT.old[i] <- PE
  sim.aug.old.cm$ATT.old.var[i] <- Var
  sim.aug.old.cm$ATT.old.lwr[i] <- Lower
  sim.aug.old.cm$ATT.old.upr[i] <- Upper
  sim.aug.old.cm$ATT.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
  
  # Homogeneous trt effect
  # Estimation using new sandwich estimator 
  # --- ATE: IPW without trimming
  result <- ATE(y=data$Y.c, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  c.sim.aug.new.cm$IPW.new[i] <- PE
  c.sim.aug.new.cm$IPW.new.var[i] <- Var
  c.sim.aug.new.cm$IPW.new.lwr[i] <- Lower
  c.sim.aug.new.cm$IPW.new.upr[i] <- Upper
  c.sim.aug.new.cm$IPW.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATC
  result <- ATC(y=data$Y.c, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  c.sim.aug.new.cm$ATC.new[i] <- PE
  c.sim.aug.new.cm$ATC.new.var[i] <- Var
  c.sim.aug.new.cm$ATC.new.lwr[i] <- Lower
  c.sim.aug.new.cm$ATC.new.upr[i] <- Upper
  c.sim.aug.new.cm$ATC.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- ATT(y=data$Y.c, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  c.sim.aug.new.cm$ATT.new[i] <- PE
  c.sim.aug.new.cm$ATT.new.var[i] <- Var
  c.sim.aug.new.cm$ATT.new.lwr[i] <- Lower
  c.sim.aug.new.cm$ATT.new.upr[i] <- Upper
  c.sim.aug.new.cm$ATT.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # Estimation using old sandwich estimator through PSweight
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  c.sim.aug.old.cm$IPW.old[i] <- PE
  c.sim.aug.old.cm$IPW.old.var[i] <- Var
  c.sim.aug.old.cm$IPW.old.lwr[i] <- Lower
  c.sim.aug.old.cm$IPW.old.upr[i] <- Upper
  c.sim.aug.old.cm$IPW.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.c", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  c.sim.aug.old.cm$ATC.old[i] <- PE
  c.sim.aug.old.cm$ATC.old.var[i] <- Var
  c.sim.aug.old.cm$ATC.old.lwr[i] <- Lower
  c.sim.aug.old.cm$ATC.old.upr[i] <- Upper
  c.sim.aug.old.cm$ATC.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  c.sim.aug.old.cm$ATT.old[i] <- PE
  c.sim.aug.old.cm$ATT.old.var[i] <- Var
  c.sim.aug.old.cm$ATT.old.lwr[i] <- Lower
  c.sim.aug.old.cm$ATT.old.upr[i] <- Upper
  c.sim.aug.old.cm$ATT.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  print(paste0("CM Ite ", i))
}

# Case 3: PS misspecified but outcome correct
sim.aug.new.mc <- results.new
sim.aug.old.mc <- results.old

c.sim.aug.new.mc <- results.new
c.sim.aug.old.mc <- results.old

ps.mult <- Z ~ X1 + X2 + X3 + X4
atc.mult <- ZZ ~ X1 + X2 + X3 + X4
out.form.c <- Y.c ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2)
out.form.h <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2) + I(X1*X3)

for(i in 1:M) {
  
  data <- weight.data.reps %>% filter(Ite == iteration[i])
  ps.cov <- as.matrix( data[,c("X1","X2","X3","X4")])
  out.cov <- as.matrix(cbind(data[,c("X1","X2","X3","X4")],
                             X1X2 = data[,"X1"]*data[,"X2"],X1sq = data[,"X1"]^2,
                             X2sq = data[,"X2"]^2,X1X3 = data[,"X1"]*data[,"X3"]))
  
  # Heterogeneous trt effect
  # Estimation using new sandwich estimator 
  # --- ATE: IPW without trimming
  result <- ATE(y=data$Y.h, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.new.mc$IPW.new[i] <- PE
  sim.aug.new.mc$IPW.new.var[i] <- Var
  sim.aug.new.mc$IPW.new.lwr[i] <- Lower
  sim.aug.new.mc$IPW.new.upr[i] <- Upper
  sim.aug.new.mc$IPW.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.h < Upper & PE.h$True.IPW.h > Lower), 1, 0))
  
  # --- ATC
  result <- ATC(y=data$Y.h, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.new.mc$ATC.new[i] <- PE
  sim.aug.new.mc$ATC.new.var[i] <- Var
  sim.aug.new.mc$ATC.new.lwr[i] <- Lower
  sim.aug.new.mc$ATC.new.upr[i] <- Upper
  sim.aug.new.mc$ATC.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))
  
  # --- ATT
  result <- ATT(y=data$Y.h, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.new.mc$ATT.new[i] <- PE
  sim.aug.new.mc$ATT.new.var[i] <- Var
  sim.aug.new.mc$ATT.new.lwr[i] <- Lower
  sim.aug.new.mc$ATT.new.upr[i] <- Upper
  sim.aug.new.mc$ATT.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
  
  # Estimation using old sandwich estimator through PSweight
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.old.mc$IPW.old[i] <- PE
  sim.aug.old.mc$IPW.old.var[i] <- Var
  sim.aug.old.mc$IPW.old.lwr[i] <- Lower
  sim.aug.old.mc$IPW.old.upr[i] <- Upper
  sim.aug.old.mc$IPW.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.h < Upper & PE.h$True.IPW.h > Lower), 1, 0))
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.h", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  sim.aug.old.mc$ATC.old[i] <- PE
  sim.aug.old.mc$ATC.old.var[i] <- Var
  sim.aug.old.mc$ATC.old.lwr[i] <- Lower
  sim.aug.old.mc$ATC.old.upr[i] <- Upper
  sim.aug.old.mc$ATC.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.old.mc$ATT.old[i] <- PE
  sim.aug.old.mc$ATT.old.var[i] <- Var
  sim.aug.old.mc$ATT.old.lwr[i] <- Lower
  sim.aug.old.mc$ATT.old.upr[i] <- Upper
  sim.aug.old.mc$ATT.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
  
  # Homogeneous trt effect
  # Estimation using new sandwich estimator 
  # --- ATE: IPW without trimming
  result <- ATE(y=data$Y.c, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  c.sim.aug.new.mc$IPW.new[i] <- PE
  c.sim.aug.new.mc$IPW.new.var[i] <- Var
  c.sim.aug.new.mc$IPW.new.lwr[i] <- Lower
  c.sim.aug.new.mc$IPW.new.upr[i] <- Upper
  c.sim.aug.new.mc$IPW.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATC
  result <- ATC(y=data$Y.c, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  c.sim.aug.new.mc$ATC.new[i] <- PE
  c.sim.aug.new.mc$ATC.new.var[i] <- Var
  c.sim.aug.new.mc$ATC.new.lwr[i] <- Lower
  c.sim.aug.new.mc$ATC.new.upr[i] <- Upper
  c.sim.aug.new.mc$ATC.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- ATT(y=data$Y.c, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  c.sim.aug.new.mc$ATT.new[i] <- PE
  c.sim.aug.new.mc$ATT.new.var[i] <- Var
  c.sim.aug.new.mc$ATT.new.lwr[i] <- Lower
  c.sim.aug.new.mc$ATT.new.upr[i] <- Upper
  c.sim.aug.new.mc$ATT.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # Estimation using old sandwich estimator through PSweight
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  c.sim.aug.old.mc$IPW.old[i] <- PE
  c.sim.aug.old.mc$IPW.old.var[i] <- Var
  c.sim.aug.old.mc$IPW.old.lwr[i] <- Lower
  c.sim.aug.old.mc$IPW.old.upr[i] <- Upper
  c.sim.aug.old.mc$IPW.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.c", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  c.sim.aug.old.mc$ATC.old[i] <- PE
  c.sim.aug.old.mc$ATC.old.var[i] <- Var
  c.sim.aug.old.mc$ATC.old.lwr[i] <- Lower
  c.sim.aug.old.mc$ATC.old.upr[i] <- Upper
  c.sim.aug.old.mc$ATC.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  c.sim.aug.old.mc$ATT.old[i] <- PE
  c.sim.aug.old.mc$ATT.old.var[i] <- Var
  c.sim.aug.old.mc$ATT.old.lwr[i] <- Lower
  c.sim.aug.old.mc$ATT.old.upr[i] <- Upper
  c.sim.aug.old.mc$ATT.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  print(paste0("MC Ite ", i))
}

# Case 4: Both PS and outcome are misspecified

sim.aug.new.mm <- results.new
sim.aug.old.mm <- results.old

c.sim.aug.new.mm <- results.new
c.sim.aug.old.mm <- results.old

ps.mult <- Z ~ X1 + X2 + X3 + X4
atc.mult <- ZZ ~ X1 + X2 + X3 + X4
out.form.c <- Y.c ~ X1 + X2 + X3 + X4 
out.form.h <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X3)

for(i in 1:M) {
  
  data <- weight.data.reps %>% filter(Ite == iteration[i])
  ps.cov <- as.matrix( data[,c("X1","X2","X3","X4")])
  out.cov <- as.matrix(cbind(data[,c("X1","X2","X3","X4")],
                             X1X3 = data[,"X1"]*data[,"X3"]))
  
  # Heterogeneous trt effect
  # Estimation using new sandwich estimator 
  # --- ATE: IPW without trimming
  result <- ATE(y=data$Y.h, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.new.mm$IPW.new[i] <- PE
  sim.aug.new.mm$IPW.new.var[i] <- Var
  sim.aug.new.mm$IPW.new.lwr[i] <- Lower
  sim.aug.new.mm$IPW.new.upr[i] <- Upper
  sim.aug.new.mm$IPW.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.h < Upper & PE.h$True.IPW.h > Lower), 1, 0))
  
  # --- ATC
  result <- ATC(y=data$Y.h, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.new.mm$ATC.new[i] <- PE
  sim.aug.new.mm$ATC.new.var[i] <- Var
  sim.aug.new.mm$ATC.new.lwr[i] <- Lower
  sim.aug.new.mm$ATC.new.upr[i] <- Upper
  sim.aug.new.mm$ATC.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))
  
  # --- ATT
  result <- ATT(y=data$Y.h, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  sim.aug.new.mm$ATT.new[i] <- PE
  sim.aug.new.mm$ATT.new.var[i] <- Var
  sim.aug.new.mm$ATT.new.lwr[i] <- Lower
  sim.aug.new.mm$ATT.new.upr[i] <- Upper
  sim.aug.new.mm$ATT.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
  
  # Estimation using old sandwich estimator through PSweight
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "IPW", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.old.mm$IPW.old[i] <- PE
  sim.aug.old.mm$IPW.old.var[i] <- Var
  sim.aug.old.mm$IPW.old.lwr[i] <- Lower
  sim.aug.old.mm$IPW.old.upr[i] <- Upper
  sim.aug.old.mm$IPW.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.IPW.h < Upper & PE.h$True.IPW.h > Lower), 1, 0))
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.h", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  sim.aug.old.mm$ATC.old[i] <- PE
  sim.aug.old.mm$ATC.old.var[i] <- Var
  sim.aug.old.mm$ATC.old.lwr[i] <- Lower
  sim.aug.old.mm$ATC.old.upr[i] <- Upper
  sim.aug.old.mm$ATC.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.old.mm$ATT.old[i] <- PE
  sim.aug.old.mm$ATT.old.var[i] <- Var
  sim.aug.old.mm$ATT.old.lwr[i] <- Lower
  sim.aug.old.mm$ATT.old.upr[i] <- Upper
  sim.aug.old.mm$ATT.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
  
  # Homogeneous trt effect
  # Estimation using new sandwich estimator 
  # --- ATE: IPW without trimming
  result <- ATE(y=data$Y.c, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  c.sim.aug.new.mm$IPW.new[i] <- PE
  c.sim.aug.new.mm$IPW.new.var[i] <- Var
  c.sim.aug.new.mm$IPW.new.lwr[i] <- Lower
  c.sim.aug.new.mm$IPW.new.upr[i] <- Upper
  c.sim.aug.new.mm$IPW.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATC
  result <- ATC(y=data$Y.c, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  c.sim.aug.new.mm$ATC.new[i] <- PE
  c.sim.aug.new.mm$ATC.new.var[i] <- Var
  c.sim.aug.new.mm$ATC.new.lwr[i] <- Lower
  c.sim.aug.new.mm$ATC.new.upr[i] <- Upper
  c.sim.aug.new.mm$ATC.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- ATT(y=data$Y.c, z=data$Z, X=ps.cov, DR=TRUE, X.out = out.cov) 
  
  PE <- result$tau 
  Var <- (result$se)^2
  Lower <- result$tau - qnorm(0.975)*result$se
  Upper <- result$tau + qnorm(0.975)*result$se
  
  c.sim.aug.new.mm$ATT.new[i] <- PE
  c.sim.aug.new.mm$ATT.new.var[i] <- Var
  c.sim.aug.new.mm$ATT.new.lwr[i] <- Lower
  c.sim.aug.new.mm$ATT.new.upr[i] <- Upper
  c.sim.aug.new.mm$ATT.new.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # Estimation using old sandwich estimator through PSweight
  # --- ATE: IPW without trimming
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "IPW", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  c.sim.aug.old.mm$IPW.old[i] <- PE
  c.sim.aug.old.mm$IPW.old.var[i] <- Var
  c.sim.aug.old.mm$IPW.old.lwr[i] <- Lower
  c.sim.aug.old.mm$IPW.old.upr[i] <- Upper
  c.sim.aug.old.mm$IPW.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.c", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  c.sim.aug.old.mm$ATC.old[i] <- PE
  c.sim.aug.old.mm$ATC.old.var[i] <- Var
  c.sim.aug.old.mm$ATC.old.lwr[i] <- Lower
  c.sim.aug.old.mm$ATC.old.upr[i] <- Upper
  c.sim.aug.old.mm$ATC.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  c.sim.aug.old.mm$ATT.old[i] <- PE
  c.sim.aug.old.mm$ATT.old.var[i] <- Var
  c.sim.aug.old.mm$ATT.old.lwr[i] <- Lower
  c.sim.aug.old.mm$ATT.old.upr[i] <- Upper
  c.sim.aug.old.mm$ATT.old.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
}

# Save data
save(file = "md3_sims.RData", 
     weight.data.reps3, sim3.p, sim3.ratio, md3.true.h,
     sim.aug.new.cc, sim.aug.old.cc, 
     sim.aug.new.cm, sim.aug.old.cm, 
     sim.aug.new.mc, sim.aug.old.mc, 
     sim.aug.new.mm, sim.aug.old.mm,
     c.sim.aug.new.cc, c.sim.aug.old.cc, 
     c.sim.aug.new.cm, c.sim.aug.old.cm, 
     c.sim.aug.new.mc, c.sim.aug.old.mc, 
     c.sim.aug.new.mm, c.sim.aug.old.mm)
