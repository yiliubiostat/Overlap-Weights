### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

####### Model 5 results plotting

### by Yi Liu
### Sept 20, 2021

load("md5_sims.RData")
source("SumStat_func.R")

# --- Plotting path
here <- getwd()
path <- paste(here, "/Results", sep="")

# Propensity score plot
ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
data <- weight.data.reps5 %>% filter(Ite == sample(Ite, 1))
ps <- SumStat(ps.formula = ps.mult, data = data, weight = "IPW")

png(filename = paste(path, "/md5_ps.png", sep=""), 
    width = 1800, height = 1300, res = 200)
plot(ps, type = "hist", breaks = 40)
dev.off()

# Point Estimation Plot
png(filename = paste(path, "/md5_PEst.png", sep=""), 
    width = 1600, height = 800, res = 200)
PEst_plot(sim5.wt.c, sim5.wt.h, sim5.aug.c.cc, sim5.aug.h.cc,
          sim5.aug.c.cm, sim5.aug.h.cm, sim5.aug.c.mc, sim5.aug.h.mc,
          sim5.aug.c.mm, sim5.aug.h.mm,
          plot.title = "")
dev.off()

# ARBias Plot
png(filename = paste(path, "/md5_ARBias.png", sep=""), 
    width = 1600, height = 800, res = 200)
abias_plot(md5.true.h, sim5.wt.c, sim5.wt.h, sim5.aug.c.cc, sim5.aug.h.cc,
           sim5.aug.c.cm, sim5.aug.h.cm, sim5.aug.c.mc, sim5.aug.h.mc,
           sim5.aug.c.mm, sim5.aug.h.mm,
           plot.title = "")
dev.off()
