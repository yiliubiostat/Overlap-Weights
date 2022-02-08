### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

####### Model 2 results plotting

### by Yi Liu
### Oct 10, 2021

load("md2_sims.RData")
source("SumStat_func.R")

# --- Plotting path
here <- getwd()
path <- paste(here, "/Results", sep="")

# Propensity score plot
ps.mult <- Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7
data <- weight.data.reps2 %>% filter(Ite == sample(Ite, 1))
ps <- SumStat(ps.formula = ps.mult, data = data, weight = "IPW")

png(filename = paste(path, "/md2_ps.png", sep=""), 
    width = 1800, height = 1300, res = 200)
plot(ps, type = "hist", breaks = 40)
dev.off()

# Point Estimation Plot
png(filename = paste(path, "/md2_PEst.png", sep=""), 
    width = 1600, height = 800, res = 200)
PEst_plot(sim2.wt.c, sim2.wt.h, sim2.aug.c.cc, sim2.aug.h.cc,
          sim2.aug.c.cm, sim2.aug.h.cm, sim2.aug.c.mc, sim2.aug.h.mc,
          sim2.aug.c.mm, sim2.aug.h.mm,
          plot.title = "")
dev.off()

# ARBias Plot
png(filename = paste(path, "/md2_ARBias.png", sep=""), 
    width = 1600, height = 800, res = 200)
abias_plot(md2.true.h, sim2.wt.c, sim2.wt.h, sim2.aug.c.cc, sim2.aug.h.cc,
           sim2.aug.c.cm, sim2.aug.h.cm, sim2.aug.c.mc, sim2.aug.h.mc,
           sim2.aug.c.mm, sim2.aug.h.mm,
           plot.title = "")
dev.off()
