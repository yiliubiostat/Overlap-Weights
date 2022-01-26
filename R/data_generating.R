### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Variance estimations ATE ATT ATC ~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### ~~~~ Generating PS and Outcome models

### by Yi Liu
### Create date: Oct 24, 2021


# Model data for simulation
source("ps_out_md_func.R")
alpha1 <- c(0.3, 0.4, 0.4, 0.4, -0.1, -0.1, 0.1)

### Data for model 1: alpha0 = -2.17
md1.simdat <- PS.model.reps(alpha0 = -2.17, alpha1)

### Data for model 2: alpha0 = -0.78
md2.simdat <- PS.model.reps(alpha0 = -0.78, alpha1)

### Data for model 3: alpha0 = 0.98
md3.simdat <- PS.model.reps(alpha0 = 0.98, alpha1)

save(file = "model_data_sim.Rdata", md1.simdat, md2.simdat, md3.simdat)


# True estimands data
### Constant effect: We use tau = 4 by model assumption

### Heterogeneous effect:
load("md_true.Rdata")
source("trueEst_func.R")

md1.true_eff <- Heter_eff(md2.dat)
md2.true_eff <- Heter_eff(md3.dat)
md3.true_eff <- Heter_eff(md4.dat)

save(file = "truth.Rdata", md1.true_eff, md2.true_eff, md3.true_eff)
