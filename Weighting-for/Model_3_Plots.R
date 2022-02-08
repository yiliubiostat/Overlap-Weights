### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

####### Model 3 results plotting

### by Yi Liu
### Sept 20, 2021

load("md3_sims.RData")
source("SumStat_func.R")

# --- Plotting path
here <- getwd()
path <- paste(here, "/Results", sep="")

# Propensity score plot
ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
data <- weight.data.reps3 %>% filter(Ite == sample(Ite, 1))
ps <- SumStat(ps.formula = ps.mult, data = data, weight = "IPW")

png(filename = paste(path, "/md3_ps.png", sep=""), 
    width = 1800, height = 1300, res = 200)
plot(ps, type = "hist", breaks = 40)
dev.off()

# Point Estimation Plot
png(filename = paste(path, "/md3_PEst.png", sep=""), 
    width = 1600, height = 800, res = 200)
PEst_plot(sim3.wt.c, sim3.wt.h, sim3.aug.c.cc, sim3.aug.h.cc,
          sim3.aug.c.cm, sim3.aug.h.cm, sim3.aug.c.mc, sim3.aug.h.mc,
          sim3.aug.c.mm, sim3.aug.h.mm,
          plot.title = "")
dev.off()

# ARBias Plot
png(filename = paste(path, "/md3_ARBias.png", sep=""), 
    width = 1600, height = 800, res = 200)
abias_plot(md3.true.h, sim3.wt.c, sim3.wt.h, sim3.aug.c.cc, sim3.aug.h.cc,
           sim3.aug.c.cm, sim3.aug.h.cm, sim3.aug.c.mc, sim3.aug.h.mc,
           sim3.aug.c.mm, sim3.aug.h.mm,
           plot.title = "")
dev.off()
