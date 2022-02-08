### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

####### Model 4 results plotting

### by Yi Liu
### Sept 20, 2021

load("md4_sims.RData")
source("SumStat_func.R")

# --- Plotting path
here <- getwd()
path <- paste(here, "/Results", sep="")

# Propensity score plot
ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
data <- weight.data.reps4 %>% filter(Ite == sample(Ite, 1))
ps <- SumStat(ps.formula = ps.mult, data = data, weight = "IPW")

png(filename = paste(path, "/md4_ps.png", sep=""), 
    width = 1800, height = 1300, res = 200)
plot(ps, type = "hist", breaks = 40)
dev.off()

# Point Estimation Plot
png(filename = paste(path, "/md4_PEst.png", sep=""), 
    width = 1600, height = 800, res = 200)
PEst_plot(sim4.wt.c, sim4.wt.h, sim4.aug.c.cc, sim4.aug.h.cc,
          sim4.aug.c.cm, sim4.aug.h.cm, sim4.aug.c.mc, sim4.aug.h.mc,
          sim4.aug.c.mm, sim4.aug.h.mm,
          plot.title = "")
dev.off()

# ARBias Plot
png(filename = paste(path, "/md4_ARBias.png", sep=""), 
    width = 1600, height = 800, res = 200)
abias_plot(md4.true.h, sim4.wt.c, sim4.wt.h, sim4.aug.c.cc, sim4.aug.h.cc,
           sim4.aug.c.cm, sim4.aug.h.cm, sim4.aug.c.mc, sim4.aug.h.mc,
           sim4.aug.c.mm, sim4.aug.h.mm,
           plot.title = "")
dev.off()
