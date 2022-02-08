### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

####### Model 1 results plotting

### by Yi Liu
### Sept 20, 2021

load("md1_sims.RData")
source("SumStat_func.R")

# --- Plotting path
here <- getwd()
path <- paste(here, "/Results", sep="")

# Propensity score plot
ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
data <- weight.data.reps1 %>% filter(Ite == sample(Ite, 1))
ps <- SumStat(ps.formula = ps.mult, data = data, weight = "IPW")

png(filename = paste(path, "/md1_ps.png", sep=""), 
    width = 1800, height = 1300, res = 200)
plot(ps, type = "hist", breaks = 40)
dev.off()

# Point Estimation Plot
png(filename = paste(path, "/md1_PEst.png", sep=""), 
    width = 1600, height = 800, res = 200)
PEst_plot(sim1.wt.c, sim1.wt.h, sim1.aug.c.cc, sim1.aug.h.cc,
          sim1.aug.c.cm, sim1.aug.h.cm, sim1.aug.c.mc, sim1.aug.h.mc,
          sim1.aug.c.mm, sim1.aug.h.mm,
          plot.title = "")
dev.off()

# ARBias Plot
png(filename = paste(path, "/md1_ARBias.png", sep=""), 
    width = 1600, height = 800, res = 200)
abias_plot(md1.true.h, sim1.wt.c, sim1.wt.h, sim1.aug.c.cc, sim1.aug.h.cc,
          sim1.aug.c.cm, sim1.aug.h.cm, sim1.aug.c.mc, sim1.aug.h.mc,
          sim1.aug.c.mm, sim1.aug.h.mm,
          plot.title = "")
dev.off()
