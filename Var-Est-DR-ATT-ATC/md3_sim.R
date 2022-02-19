### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Variance estimations ATE ATT ATC ~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study                 ~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### by Yi Liu
### Oct 24, 2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~ Using M Replicates for Estimands ~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load("model_data_sim.Rdata")
load("truth.Rdata")

weight.data.reps <- md3.simdat
M <- length(unique(weight.data.reps$Ite))

results.c <- data.frame(IPW.c = rep(NA, M),
                        IPW.c.var = rep(NA, M),
                        IPW.c.lwr = rep(NA, M),
                        IPW.c.upr = rep(NA, M),
                        IPW.c.ifci = rep(NA, M),
                        
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
# ~~~~~~~~~~~~~~~~~~~~~~ Augmented Estimators ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Case 1: Both PS and outcome are correct
PE.c <- 4
PE.h <- md3.true_eff

sim.aug.c.cc <- results.c
sim.aug.h.cc <- results.h

ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
atc.mult <- ZZ ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
out.form.c <- Y.c ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2)
out.form.h <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2) + I(X1*X3)

for(i in 1:M) {
  
  data <- weight.data.reps %>% filter(Ite == i)
  
  # Constant model
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
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.c", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  sim.aug.c.cc$ATC.c[i] <- PE
  sim.aug.c.cc$ATC.c.var[i] <- Var
  sim.aug.c.cc$ATC.c.lwr[i] <- Lower
  sim.aug.c.cc$ATC.c.upr[i] <- Upper
  sim.aug.c.cc$ATC.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.cc$ATT.c[i] <- PE
  sim.aug.c.cc$ATT.c.var[i] <- Var
  sim.aug.c.cc$ATT.c.lwr[i] <- Lower
  sim.aug.c.cc$ATT.c.upr[i] <- Upper
  sim.aug.c.cc$ATT.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # Hetergeneous model
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
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.h", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  sim.aug.h.cc$ATC.h[i] <- PE
  sim.aug.h.cc$ATC.h.var[i] <- Var
  sim.aug.h.cc$ATC.h.lwr[i] <- Lower
  sim.aug.h.cc$ATC.h.upr[i] <- Upper
  sim.aug.h.cc$ATC.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))

  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.cc$ATT.h[i] <- PE
  sim.aug.h.cc$ATT.h.var[i] <- Var
  sim.aug.h.cc$ATT.h.lwr[i] <- Lower
  sim.aug.h.cc$ATT.h.upr[i] <- Upper
  sim.aug.h.cc$ATT.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
}

# Case 2: PS correct but outcome misspecified
sim.aug.c.cm <- results.c
sim.aug.h.cm <- results.h

ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
atc.mult <- ZZ ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
out.form.c <- Y.c ~ X1 + X2 + X3 + X4
out.form.h <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X3)

for(i in 1:M) {
  
  data <- weight.data.reps %>% filter(Ite == i)
  
  # Constant model
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
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.c", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  sim.aug.c.cm$ATC.c[i] <- PE
  sim.aug.c.cm$ATC.c.var[i] <- Var
  sim.aug.c.cm$ATC.c.lwr[i] <- Lower
  sim.aug.c.cm$ATC.c.upr[i] <- Upper
  sim.aug.c.cm$ATC.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.cm$ATT.c[i] <- PE
  sim.aug.c.cm$ATT.c.var[i] <- Var
  sim.aug.c.cm$ATT.c.lwr[i] <- Lower
  sim.aug.c.cm$ATT.c.upr[i] <- Upper
  sim.aug.c.cm$ATT.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # Hetergeneous model
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
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.h", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  sim.aug.h.cm$ATC.h[i] <- PE
  sim.aug.h.cm$ATC.h.var[i] <- Var
  sim.aug.h.cm$ATC.h.lwr[i] <- Lower
  sim.aug.h.cm$ATC.h.upr[i] <- Upper
  sim.aug.h.cm$ATC.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.cm$ATT.h[i] <- PE
  sim.aug.h.cm$ATT.h.var[i] <- Var
  sim.aug.h.cm$ATT.h.lwr[i] <- Lower
  sim.aug.h.cm$ATT.h.upr[i] <- Upper
  sim.aug.h.cm$ATT.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
}

# Case 3: PS misspecified but outcome correct
sim.aug.c.mc <- results.c
sim.aug.h.mc <- results.h

ps.mult <- Z ~ X1 + X2 + X3 + X4
atc.mult <- ZZ ~ X1 + X2 + X3 + X4
out.form.c <- Y.c ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2)
out.form.h <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X2) + I(X1^2) + I(X2^2) + I(X1*X3)

for(i in 1:M) {
  
  data <- weight.data.reps %>% filter(Ite == i)
  
  # Constant model
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
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.c", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  sim.aug.c.mc$ATC.c[i] <- PE
  sim.aug.c.mc$ATC.c.var[i] <- Var
  sim.aug.c.mc$ATC.c.lwr[i] <- Lower
  sim.aug.c.mc$ATC.c.upr[i] <- Upper
  sim.aug.c.mc$ATC.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.mc$ATT.c[i] <- PE
  sim.aug.c.mc$ATT.c.var[i] <- Var
  sim.aug.c.mc$ATT.c.lwr[i] <- Lower
  sim.aug.c.mc$ATT.c.upr[i] <- Upper
  sim.aug.c.mc$ATT.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # Hetergeneous model
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
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.h", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  sim.aug.h.mc$ATC.h[i] <- PE
  sim.aug.h.mc$ATC.h.var[i] <- Var
  sim.aug.h.mc$ATC.h.lwr[i] <- Lower
  sim.aug.h.mc$ATC.h.upr[i] <- Upper
  sim.aug.h.mc$ATC.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.mc$ATT.h[i] <- PE
  sim.aug.h.mc$ATT.h.var[i] <- Var
  sim.aug.h.mc$ATT.h.lwr[i] <- Lower
  sim.aug.h.mc$ATT.h.upr[i] <- Upper
  sim.aug.h.mc$ATT.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
}

# Case 4: Both PS and outcome are misspecified
sim.aug.c.mm <- results.c
sim.aug.h.mm <- results.h

ps.mult <- Z ~ X1 + X2 + X3 + X4
atc.mult <- ZZ ~ X1 + X2 + X3 + X4
out.form.c <- Y.c ~ X1 + X2 + X3 + X4 
out.form.h <- Y.h ~ X1 + X2 + X3 + X4 + I(X1*X3)

for(i in 1:M) {
  
  data <- weight.data.reps %>% filter(Ite == i)
  
  # Constant model
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
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.c", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  sim.aug.c.mm$ATC.c[i] <- PE
  sim.aug.c.mm$ATC.c.var[i] <- Var
  sim.aug.c.mm$ATC.c.lwr[i] <- Lower
  sim.aug.c.mm$ATC.c.upr[i] <- Upper
  sim.aug.c.mm$ATC.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.c", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.c, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.c.mm$ATT.c[i] <- PE
  sim.aug.c.mm$ATT.c.var[i] <- Var
  sim.aug.c.mm$ATT.c.lwr[i] <- Lower
  sim.aug.c.mm$ATT.c.upr[i] <- Upper
  sim.aug.c.mm$ATT.c.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.c < Upper & PE.c > Lower), 1, 0))
  
  # Hetergeneous model
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
  
  # --- ATC
  result <- summary(PSweight(ps.formula = atc.mult, yname = "Y.h", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- -(result$estimates[1])
  Var <- (result$estimates[2])^2
  Lower <- -(result$estimates[5])
  Upper <- -(result$estimates[4])
  
  sim.aug.h.mm$ATC.h[i] <- PE
  sim.aug.h.mm$ATC.h.var[i] <- Var
  sim.aug.h.mm$ATC.h.lwr[i] <- Lower
  sim.aug.h.mm$ATC.h.upr[i] <- Upper
  sim.aug.h.mm$ATC.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATC.h < Upper & PE.h$True.ATC.h > Lower), 1, 0))
  
  # --- ATT
  result <- summary(PSweight(ps.formula = ps.mult, yname = "Y.h", weight = "treated", data = data,
                             augmentation = TRUE, out.formula = out.form.h, family = "gaussian"), type = "DIF")
  PE <- result$estimates[1]
  Var <- (result$estimates[2])^2
  Lower <- result$estimates[4]
  Upper <- result$estimates[5]
  
  sim.aug.h.mm$ATT.h[i] <- PE
  sim.aug.h.mm$ATT.h.var[i] <- Var
  sim.aug.h.mm$ATT.h.lwr[i] <- Lower
  sim.aug.h.mm$ATT.h.upr[i] <- Upper
  sim.aug.h.mm$ATT.h.ifci[i] <- ifelse(is.na(Var), NA, ifelse((PE.h$True.ATT.h < Upper & PE.h$True.ATT.h > Lower), 1, 0))
}

# One way to address NA problem in variance and CIs
num.na <- function(x) { sum(is.na(x)) }
M0 <- round(M/1.5)

all.sim.dat <- cbind(sim.aug.c.cc, sim.aug.h.cc, 
                     sim.aug.c.cm, sim.aug.h.cm, 
                     sim.aug.c.mc, sim.aug.h.mc, 
                     sim.aug.c.mm, sim.aug.h.mm) 

colnames(all.sim.dat) = as.character(1:128)
all.sim.dat <- all.sim.dat %>% mutate(Num.NA = rep(0, M), Row.Num = 1:M)
all.sim.dat$Num.NA <- apply(all.sim.dat, 1, num.na)

all.sim.dat <- all.sim.dat %>% filter(Num.NA==0)

R <- all.sim.dat$Row.Num
if(length(all.sim.dat$Row.Num) < M0) { warning("The number of rows < M0!") } else {
  R <- sample(all.sim.dat$Row.Num, size = M0, replace = FALSE)
}

sim3.aug.c.cc <- sim.aug.c.cc %>% filter(Row.Num %in% R)
sim3.aug.h.cc <- sim.aug.h.cc %>% filter(Row.Num %in% R)

sim3.aug.c.cm <- sim.aug.c.cm %>% filter(Row.Num %in% R)
sim3.aug.h.cm <- sim.aug.h.cm %>% filter(Row.Num %in% R)

sim3.aug.c.mc <- sim.aug.c.mc %>% filter(Row.Num %in% R)
sim3.aug.h.mc <- sim.aug.h.mc %>% filter(Row.Num %in% R)

sim3.aug.c.mm <- sim.aug.c.mm %>% filter(Row.Num %in% R)
sim3.aug.h.mm <- sim.aug.h.mm %>% filter(Row.Num %in% R)

# Design module
weight.data.reps3 <- weight.data.reps %>% filter(Ite %in% R)

# --- Propensity score plot will be draw in SumStat file
# --- p
sim3.p <- length(weight.data.reps3$Z[weight.data.reps3$Z == 1])/(M0*N)

# --- ratio
ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7

ps.var.ctrl <- rep(NA, M0)
ps.var.trt <- rep(NA, M0)

k = 1
for(i in R) {
  ps.dat <- data.frame(ps = rep(NA, N), Z = rep(NA, N))
  data <- weight.data.reps3 %>% filter(Ite == i)
  
  ps <- SumStat(ps.formula = ps.mult, data = data, weight = "overlap")$propensity
  ps.dat$ps <- ps[, "1"]
  ps.dat$Z <- data$Z
  
  ps.var.ctrl[k] = ps.dat %>% filter(Z==0) %>% summarize(var(ps)) %>% as.numeric()
  ps.var.trt[k] = ps.dat %>% filter(Z==1) %>% summarize(var(ps)) %>% as.numeric()
  k = k+1
}

sim3.ratio <- mean(ps.var.trt/ps.var.ctrl)

md3.true.h <- PE.h

# Save data
save(file = "md3_sims_PSwt.RData", 
     weight.data.reps3, sim3.p, sim3.ratio, md3.true.h,
     sim3.aug.c.cc, sim3.aug.h.cc, sim3.aug.c.cm, sim3.aug.h.cm, 
     sim3.aug.c.mc, sim3.aug.h.mc, sim3.aug.c.mm, sim3.aug.h.mm )
