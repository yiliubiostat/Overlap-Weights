### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Variance estimations ATE ATT ATC ~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### by Yi Liu
### Nov 17, 2021

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~ Wild Bootstrap Simulation        ~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~ For constant treatment effect    ~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load("md1_sims_PSwt.RData")
source("wild_boot_func.R")
library(dplyr)

weight.data.reps <- weight.data.reps1 %>% mutate(H = X1*X3)
M <- length(unique(weight.data.reps$Ite))
N <- nrow(weight.data.reps)/length(unique(weight.data.reps$Ite))
weight.data.reps <- weight.data.reps %>% select(-Ite,-ZZ) %>% mutate(Ite = rep(1:M, rep(N,M)))

results <- data.frame(ATE.wbvar = rep(NA, M),
                      ATE = rep(NA, M),
                      ATE.lwr = rep(NA, M),
                      ATE.upr = rep(NA, M),
                      ATE.ifci = rep(NA, M),
                      
                      ATE.biasC = rep(NA, M),
                      ATE.biasC.lwr = rep(NA, M),
                      ATE.biasC.upr = rep(NA, M),
                      ATE.biasC.ifci = rep(NA, M),
                      
                      ATC.wbvar = rep(NA, M),
                      ATC = rep(NA, M),
                      ATC.lwr = rep(NA, M),
                      ATC.upr = rep(NA, M),
                      ATC.ifci = rep(NA, M),
                      
                      ATC.biasC = rep(NA, M),
                      ATC.biasC.lwr = rep(NA, M),
                      ATC.biasC.upr = rep(NA, M),
                      ATC.biasC.ifci = rep(NA, M),
                      
                      ATT.wbvar = rep(NA, M),
                      ATT = rep(NA, M),
                      ATT.lwr = rep(NA, M),
                      ATT.upr = rep(NA, M),
                      ATT.ifci = rep(NA, M),
                      
                      ATT.biasC = rep(NA, M),
                      ATT.biasC.lwr = rep(NA, M),
                      ATT.biasC.upr = rep(NA, M),
                      ATT.biasC.ifci = rep(NA, M) )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~ Augmented Estimators ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Case 1: Both PS and outcome are correct
sim.rad.h.cc <- results
sim.exp.h.cc <- results

for(i in 1:M) {
  
  data <- weight.data.reps %>% filter(Ite == i)
  
  y = data$Y.c
  z = data$Z
  X = data[,c("X1","X2","X3","X4","X5","X6","X7")] %>% as.matrix()
  X.out = data[,c("X1","X2","X3","X4","X5","X6","X7","H")] %>% as.matrix()
  
  # --- ATE: IPW without trimming
  # --- Using Rademacher distribution
  result <- ATE.boot(y=y, z=z, X=X, X.out=X.out)
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.rad.h.cc$ATE.wbvar[i] <- Var
  
  sim.rad.h.cc$ATE[i] <- PE
  sim.rad.h.cc$ATE.lwr[i] <- Lower
  sim.rad.h.cc$ATE.upr[i] <- Upper
  sim.rad.h.cc$ATE.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.rad.h.cc$ATE.biasC[i] <- PE.bc
  sim.rad.h.cc$ATE.biasC.lwr[i] <- Lower.bc
  sim.rad.h.cc$ATE.biasC.upr[i] <- Upper.bc
  sim.rad.h.cc$ATE.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- Using exponential distribution
  result <- ATE.boot(y=y, z=z, X=X, X.out=X.out, RV="Exp")
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.exp.h.cc$ATE.wbvar[i] <- Var
  
  sim.exp.h.cc$ATE[i] <- PE
  sim.exp.h.cc$ATE.lwr[i] <- Lower
  sim.exp.h.cc$ATE.upr[i] <- Upper
  sim.exp.h.cc$ATE.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.exp.h.cc$ATE.biasC[i] <- PE.bc
  sim.exp.h.cc$ATE.biasC.lwr[i] <- Lower.bc
  sim.exp.h.cc$ATE.biasC.upr[i] <- Upper.bc
  sim.exp.h.cc$ATE.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- ATC
  # --- Using Rademacher distribution
  result <- ATC.boot(y=y, z=z, X=X, X.out=X.out)
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.rad.h.cc$ATC.wbvar[i] <- Var
  
  sim.rad.h.cc$ATC[i] <- PE
  sim.rad.h.cc$ATC.lwr[i] <- Lower
  sim.rad.h.cc$ATC.upr[i] <- Upper
  sim.rad.h.cc$ATC.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.rad.h.cc$ATC.biasC[i] <- PE.bc
  sim.rad.h.cc$ATC.biasC.lwr[i] <- Lower.bc
  sim.rad.h.cc$ATC.biasC.upr[i] <- Upper.bc
  sim.rad.h.cc$ATC.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- Using exponential distribution
  result <- ATC.boot(y=y, z=z, X=X, X.out=X.out, RV="Exp")
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.exp.h.cc$ATC.wbvar[i] <- Var
  
  sim.exp.h.cc$ATC[i] <- PE
  sim.exp.h.cc$ATC.lwr[i] <- Lower
  sim.exp.h.cc$ATC.upr[i] <- Upper
  sim.exp.h.cc$ATC.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.exp.h.cc$ATC.biasC[i] <- PE.bc
  sim.exp.h.cc$ATC.biasC.lwr[i] <- Lower.bc
  sim.exp.h.cc$ATC.biasC.upr[i] <- Upper.bc
  sim.exp.h.cc$ATC.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- ATT
  # --- Using Rademacher distribution
  result <- ATT.boot(y=y, z=z, X=X, X.out=X.out)
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.rad.h.cc$ATT.wbvar[i] <- Var
  
  sim.rad.h.cc$ATT[i] <- PE
  sim.rad.h.cc$ATT.lwr[i] <- Lower
  sim.rad.h.cc$ATT.upr[i] <- Upper
  sim.rad.h.cc$ATT.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.rad.h.cc$ATT.biasC[i] <- PE.bc
  sim.rad.h.cc$ATT.biasC.lwr[i] <- Lower.bc
  sim.rad.h.cc$ATT.biasC.upr[i] <- Upper.bc
  sim.rad.h.cc$ATT.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- Using exponential distribution
  result <- ATT.boot(y=y, z=z, X=X, X.out=X.out, RV="Exp")
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.exp.h.cc$ATT.wbvar[i] <- Var
  
  sim.exp.h.cc$ATT[i] <- PE
  sim.exp.h.cc$ATT.lwr[i] <- Lower
  sim.exp.h.cc$ATT.upr[i] <- Upper
  sim.exp.h.cc$ATT.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.exp.h.cc$ATT.biasC[i] <- PE.bc
  sim.exp.h.cc$ATT.biasC.lwr[i] <- Lower.bc
  sim.exp.h.cc$ATT.biasC.upr[i] <- Upper.bc
  sim.exp.h.cc$ATT.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
}

# Case 2: PS correct but outcome misspecified
sim.rad.h.cm <- results
sim.exp.h.cm <- results

for(i in 1:M) {
  
  data <- weight.data.reps %>% filter(Ite == i)
  
  y = data$Y.c
  z = data$Z
  X = data[,c("X1","X2","X3","X4","X5","X6","X7")] %>% as.matrix()
  X.out = data[,c("X1","X2","X3","X4","H")] %>% as.matrix()
  
  # --- ATE: IPW without trimming
  # --- Using Rademacher distribution
  result <- ATE.boot(y=y, z=z, X=X, X.out=X.out)
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.rad.h.cm$ATE.wbvar[i] <- Var
  
  sim.rad.h.cm$ATE[i] <- PE
  sim.rad.h.cm$ATE.lwr[i] <- Lower
  sim.rad.h.cm$ATE.upr[i] <- Upper
  sim.rad.h.cm$ATE.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.rad.h.cm$ATE.biasC[i] <- PE.bc
  sim.rad.h.cm$ATE.biasC.lwr[i] <- Lower.bc
  sim.rad.h.cm$ATE.biasC.upr[i] <- Upper.bc
  sim.rad.h.cm$ATE.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- Using exponential distribution
  result <- ATE.boot(y=y, z=z, X=X, X.out=X.out, RV="Exp")
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.exp.h.cm$ATE.wbvar[i] <- Var
  
  sim.exp.h.cm$ATE[i] <- PE
  sim.exp.h.cm$ATE.lwr[i] <- Lower
  sim.exp.h.cm$ATE.upr[i] <- Upper
  sim.exp.h.cm$ATE.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.exp.h.cm$ATE.biasC[i] <- PE.bc
  sim.exp.h.cm$ATE.biasC.lwr[i] <- Lower.bc
  sim.exp.h.cm$ATE.biasC.upr[i] <- Upper.bc
  sim.exp.h.cm$ATE.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- ATC
  # --- Using Rademacher distribution
  result <- ATC.boot(y=y, z=z, X=X, X.out=X.out)
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.rad.h.cm$ATC.wbvar[i] <- Var
  
  sim.rad.h.cm$ATC[i] <- PE
  sim.rad.h.cm$ATC.lwr[i] <- Lower
  sim.rad.h.cm$ATC.upr[i] <- Upper
  sim.rad.h.cm$ATC.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.rad.h.cm$ATC.biasC[i] <- PE.bc
  sim.rad.h.cm$ATC.biasC.lwr[i] <- Lower.bc
  sim.rad.h.cm$ATC.biasC.upr[i] <- Upper.bc
  sim.rad.h.cm$ATC.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- Using exponential distribution
  result <- ATC.boot(y=y, z=z, X=X, X.out=X.out, RV="Exp")
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.exp.h.cm$ATC.wbvar[i] <- Var
  
  sim.exp.h.cm$ATC[i] <- PE
  sim.exp.h.cm$ATC.lwr[i] <- Lower
  sim.exp.h.cm$ATC.upr[i] <- Upper
  sim.exp.h.cm$ATC.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.exp.h.cm$ATC.biasC[i] <- PE.bc
  sim.exp.h.cm$ATC.biasC.lwr[i] <- Lower.bc
  sim.exp.h.cm$ATC.biasC.upr[i] <- Upper.bc
  sim.exp.h.cm$ATC.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- ATT
  # --- Using Rademacher distribution
  result <- ATT.boot(y=y, z=z, X=X, X.out=X.out)
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.rad.h.cm$ATT.wbvar[i] <- Var
  
  sim.rad.h.cm$ATT[i] <- PE
  sim.rad.h.cm$ATT.lwr[i] <- Lower
  sim.rad.h.cm$ATT.upr[i] <- Upper
  sim.rad.h.cm$ATT.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.rad.h.cm$ATT.biasC[i] <- PE.bc
  sim.rad.h.cm$ATT.biasC.lwr[i] <- Lower.bc
  sim.rad.h.cm$ATT.biasC.upr[i] <- Upper.bc
  sim.rad.h.cm$ATT.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- Using exponential distribution
  result <- ATT.boot(y=y, z=z, X=X, X.out=X.out, RV="Exp")
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.exp.h.cm$ATT.wbvar[i] <- Var
  
  sim.exp.h.cm$ATT[i] <- PE
  sim.exp.h.cm$ATT.lwr[i] <- Lower
  sim.exp.h.cm$ATT.upr[i] <- Upper
  sim.exp.h.cm$ATT.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.exp.h.cm$ATT.biasC[i] <- PE.bc
  sim.exp.h.cm$ATT.biasC.lwr[i] <- Lower.bc
  sim.exp.h.cm$ATT.biasC.upr[i] <- Upper.bc
  sim.exp.h.cm$ATT.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
}

# Case 3: PS misspecified but outcome correct
sim.rad.h.mc <- results
sim.exp.h.mc <- results

for(i in 1:M) {
  
  data <- weight.data.reps %>% filter(Ite == i)
  
  y = data$Y.c
  z = data$Z
  X = data[,c("X1","X2","X3","X4")] %>% as.matrix()
  X.out = data[,c("X1","X2","X3","X4","X5","X6","X7","H")] %>% as.matrix()
  
  # --- ATE: IPW without trimming
  # --- Using Rademacher distribution
  result <- ATE.boot(y=y, z=z, X=X, X.out=X.out)
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.rad.h.mc$ATE.wbvar[i] <- Var
  
  sim.rad.h.mc$ATE[i] <- PE
  sim.rad.h.mc$ATE.lwr[i] <- Lower
  sim.rad.h.mc$ATE.upr[i] <- Upper
  sim.rad.h.mc$ATE.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.rad.h.mc$ATE.biasC[i] <- PE.bc
  sim.rad.h.mc$ATE.biasC.lwr[i] <- Lower.bc
  sim.rad.h.mc$ATE.biasC.upr[i] <- Upper.bc
  sim.rad.h.mc$ATE.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- Using exponential distribution
  result <- ATE.boot(y=y, z=z, X=X, X.out=X.out, RV="Exp")
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.exp.h.mc$ATE.wbvar[i] <- Var
  
  sim.exp.h.mc$ATE[i] <- PE
  sim.exp.h.mc$ATE.lwr[i] <- Lower
  sim.exp.h.mc$ATE.upr[i] <- Upper
  sim.exp.h.mc$ATE.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.exp.h.mc$ATE.biasC[i] <- PE.bc
  sim.exp.h.mc$ATE.biasC.lwr[i] <- Lower.bc
  sim.exp.h.mc$ATE.biasC.upr[i] <- Upper.bc
  sim.exp.h.mc$ATE.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- ATC
  # --- Using Rademacher distribution
  result <- ATC.boot(y=y, z=z, X=X, X.out=X.out)
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.rad.h.mc$ATC.wbvar[i] <- Var
  
  sim.rad.h.mc$ATC[i] <- PE
  sim.rad.h.mc$ATC.lwr[i] <- Lower
  sim.rad.h.mc$ATC.upr[i] <- Upper
  sim.rad.h.mc$ATC.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.rad.h.mc$ATC.biasC[i] <- PE.bc
  sim.rad.h.mc$ATC.biasC.lwr[i] <- Lower.bc
  sim.rad.h.mc$ATC.biasC.upr[i] <- Upper.bc
  sim.rad.h.mc$ATC.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- Using exponential distribution
  result <- ATC.boot(y=y, z=z, X=X, X.out=X.out, RV="Exp")
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.exp.h.mc$ATC.wbvar[i] <- Var
  
  sim.exp.h.mc$ATC[i] <- PE
  sim.exp.h.mc$ATC.lwr[i] <- Lower
  sim.exp.h.mc$ATC.upr[i] <- Upper
  sim.exp.h.mc$ATC.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.exp.h.mc$ATC.biasC[i] <- PE.bc
  sim.exp.h.mc$ATC.biasC.lwr[i] <- Lower.bc
  sim.exp.h.mc$ATC.biasC.upr[i] <- Upper.bc
  sim.exp.h.mc$ATC.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- ATT
  # --- Using Rademacher distribution
  result <- ATT.boot(y=y, z=z, X=X, X.out=X.out)
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.rad.h.mc$ATT.wbvar[i] <- Var
  
  sim.rad.h.mc$ATT[i] <- PE
  sim.rad.h.mc$ATT.lwr[i] <- Lower
  sim.rad.h.mc$ATT.upr[i] <- Upper
  sim.rad.h.mc$ATT.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.rad.h.mc$ATT.biasC[i] <- PE.bc
  sim.rad.h.mc$ATT.biasC.lwr[i] <- Lower.bc
  sim.rad.h.mc$ATT.biasC.upr[i] <- Upper.bc
  sim.rad.h.mc$ATT.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- Using exponential distribution
  result <- ATT.boot(y=y, z=z, X=X, X.out=X.out, RV="Exp")
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.exp.h.mc$ATT.wbvar[i] <- Var
  
  sim.exp.h.mc$ATT[i] <- PE
  sim.exp.h.mc$ATT.lwr[i] <- Lower
  sim.exp.h.mc$ATT.upr[i] <- Upper
  sim.exp.h.mc$ATT.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.exp.h.mc$ATT.biasC[i] <- PE.bc
  sim.exp.h.mc$ATT.biasC.lwr[i] <- Lower.bc
  sim.exp.h.mc$ATT.biasC.upr[i] <- Upper.bc
  sim.exp.h.mc$ATT.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
}

# Case 4: Both PS and outcome are misspecified
sim.rad.h.mm <- results
sim.exp.h.mm <- results

for(i in 1:M) {
  
  data <- weight.data.reps %>% filter(Ite == i)
  
  y = data$Y.c
  z = data$Z
  X = data[,c("X1","X2","X3","X4")] %>% as.matrix()
  X.out = data[,c("X1","X2","X3","X4","H")] %>% as.matrix()
  
  # --- ATE: IPW without trimming
  # --- Using Rademacher distribution
  result <- ATE.boot(y=y, z=z, X=X, X.out=X.out)
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.rad.h.mm$ATE.wbvar[i] <- Var
  
  sim.rad.h.mm$ATE[i] <- PE
  sim.rad.h.mm$ATE.lwr[i] <- Lower
  sim.rad.h.mm$ATE.upr[i] <- Upper
  sim.rad.h.mm$ATE.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.rad.h.mm$ATE.biasC[i] <- PE.bc
  sim.rad.h.mm$ATE.biasC.lwr[i] <- Lower.bc
  sim.rad.h.mm$ATE.biasC.upr[i] <- Upper.bc
  sim.rad.h.mm$ATE.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- Using exponential distribution
  result <- ATE.boot(y=y, z=z, X=X, X.out=X.out, RV="Exp")
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.exp.h.mm$ATE.wbvar[i] <- Var
  
  sim.exp.h.mm$ATE[i] <- PE
  sim.exp.h.mm$ATE.lwr[i] <- Lower
  sim.exp.h.mm$ATE.upr[i] <- Upper
  sim.exp.h.mm$ATE.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.exp.h.mm$ATE.biasC[i] <- PE.bc
  sim.exp.h.mm$ATE.biasC.lwr[i] <- Lower.bc
  sim.exp.h.mm$ATE.biasC.upr[i] <- Upper.bc
  sim.exp.h.mm$ATE.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- ATC
  # --- Using Rademacher distribution
  result <- ATC.boot(y=y, z=z, X=X, X.out=X.out)
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.rad.h.mm$ATC.wbvar[i] <- Var
  
  sim.rad.h.mm$ATC[i] <- PE
  sim.rad.h.mm$ATC.lwr[i] <- Lower
  sim.rad.h.mm$ATC.upr[i] <- Upper
  sim.rad.h.mm$ATC.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.rad.h.mm$ATC.biasC[i] <- PE.bc
  sim.rad.h.mm$ATC.biasC.lwr[i] <- Lower.bc
  sim.rad.h.mm$ATC.biasC.upr[i] <- Upper.bc
  sim.rad.h.mm$ATC.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- Using exponential distribution
  result <- ATC.boot(y=y, z=z, X=X, X.out=X.out, RV="Exp")
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.exp.h.mm$ATC.wbvar[i] <- Var
  
  sim.exp.h.mm$ATC[i] <- PE
  sim.exp.h.mm$ATC.lwr[i] <- Lower
  sim.exp.h.mm$ATC.upr[i] <- Upper
  sim.exp.h.mm$ATC.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.exp.h.mm$ATC.biasC[i] <- PE.bc
  sim.exp.h.mm$ATC.biasC.lwr[i] <- Lower.bc
  sim.exp.h.mm$ATC.biasC.upr[i] <- Upper.bc
  sim.exp.h.mm$ATC.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- ATT
  # --- Using Rademacher distribution
  result <- ATT.boot(y=y, z=z, X=X, X.out=X.out)
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.rad.h.mm$ATT.wbvar[i] <- Var
  
  sim.rad.h.mm$ATT[i] <- PE
  sim.rad.h.mm$ATT.lwr[i] <- Lower
  sim.rad.h.mm$ATT.upr[i] <- Upper
  sim.rad.h.mm$ATT.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.rad.h.mm$ATT.biasC[i] <- PE.bc
  sim.rad.h.mm$ATT.biasC.lwr[i] <- Lower.bc
  sim.rad.h.mm$ATT.biasC.upr[i] <- Upper.bc
  sim.rad.h.mm$ATT.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
  
  # --- Using exponential distribution
  result <- ATT.boot(y=y, z=z, X=X, X.out=X.out, RV="Exp")
  Var <- (result$Std.Boot)^2
  
  PE <- result$PE
  Lower <- result$PE.lwr
  Upper <- result$PE.upr
  
  PE.bc <- result$biasC.Est
  Lower.bc <- result$biasC.lwr
  Upper.bc <- result$biasC.upr
  
  sim.exp.h.mm$ATT.wbvar[i] <- Var
  
  sim.exp.h.mm$ATT[i] <- PE
  sim.exp.h.mm$ATT.lwr[i] <- Lower
  sim.exp.h.mm$ATT.upr[i] <- Upper
  sim.exp.h.mm$ATT.ifci[i] <- ifelse((4 < Upper & 4 > Lower), 1, 0)
  
  sim.exp.h.mm$ATT.biasC[i] <- PE.bc
  sim.exp.h.mm$ATT.biasC.lwr[i] <- Lower.bc
  sim.exp.h.mm$ATT.biasC.upr[i] <- Upper.bc
  sim.exp.h.mm$ATT.biasC.ifci[i] <- ifelse((4 < Upper.bc & 4 > Lower.bc), 1, 0)
}

# Summarize results
wb1.rad.c.cc <- sim.rad.h.cc
wb1.rad.c.cm <- sim.rad.h.cm
wb1.rad.c.mc <- sim.rad.h.mc
wb1.rad.c.mm <- sim.rad.h.mm

wb1.exp.c.cc <- sim.exp.h.cc
wb1.exp.c.cm <- sim.exp.h.cm
wb1.exp.c.mc <- sim.exp.h.mc
wb1.exp.c.mm <- sim.exp.h.mm

# Save data
save(file = "md1_sims_WBoot_cons.RData", 
     wb1.rad.c.cc, wb1.rad.c.cm, wb1.rad.c.mc, wb1.rad.c.mm, 
     wb1.exp.c.cc, wb1.exp.c.cm, wb1.exp.c.mc, wb1.exp.c.mm )
