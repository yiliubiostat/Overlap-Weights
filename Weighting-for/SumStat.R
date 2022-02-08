### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### Summary statistics and LaTeX output

### by Yi Liu
### Oct 25, 2021

library(dplyr)
library(xtable)

source("SumStat_func.R")
source("ESS_func.R")

# Load the data
load("md1_sims.RData")
load("md2_sims.RData")
load("md3_sims.RData")
load("md4_sims.RData")
load("md5_sims.RData")
load("md_true.RData")

###### Truth table

# --- Constant effect: 4
# --- Heterogeneous effect:
cols = c("p (%)", "ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", 
              "ATO", "ATM", "ATEN", "ATC", "ATT")

p <- round(100*c(mean(md1.dat$Z), mean(md2.dat$Z), mean(md3.dat$Z), mean(md4.dat$Z), mean(md5.dat$Z)), 2)
true.df <- cbind(p,
                 round(rbind(md1.true.h, md2.true.h, md3.true.h, md4.true.h, md5.true.h), 2))
colnames(true.df) <- cols
rownames(true.df) <- NULL
xtable(true.df)

###### Point estimations, Bias, RMSE, RE, CP

rows = c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", "ATO", "ATM", "ATEN", "ATC", "ATT")
cols = c("Point Est", "ARBias (%)", "RMSE", "RE", "CP")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Model 1
sim1.ratio
sim1.p

######### ESS
ess <- ESS(weight.data.reps1)
xtable(ess)

######### PE, ARBias, RMSE, RE, CP
ac <- PE(sim1.wt.c)
bc <- PE(sim1.aug.c.cc)
cc <- PE(sim1.aug.c.cm)
dc <- PE(sim1.aug.c.mc)
ec <- PE(sim1.aug.c.mm)

ah <- PE(sim1.wt.h)
bh <- PE(sim1.aug.h.cc)
ch <- PE(sim1.aug.h.cm)
dh <- PE(sim1.aug.h.mc)
eh <- PE(sim1.aug.h.mm)

a1 <- avg.rbias(true.h = md1.true.h, simdata.c = sim1.wt.c, simdata.h = sim1.wt.h)
b1 <- avg.rbias(true.h = md1.true.h, simdata.c = sim1.aug.c.cc, simdata.h = sim1.aug.h.cc)
c1 <- avg.rbias(true.h = md1.true.h, simdata.c = sim1.aug.c.cm, simdata.h = sim1.aug.h.cm)
d1 <- avg.rbias(true.h = md1.true.h, simdata.c = sim1.aug.c.mc, simdata.h = sim1.aug.h.mc)
e1 <- avg.rbias(true.h = md1.true.h, simdata.c = sim1.aug.c.mm, simdata.h = sim1.aug.h.mm)

a2 <- RMSE(true.h = md1.true.h, simdata.c = sim1.wt.c, simdata.h = sim1.wt.h)
b2 <- RMSE(true.h = md1.true.h, simdata.c = sim1.aug.c.cc, simdata.h = sim1.aug.h.cc)
c2 <- RMSE(true.h = md1.true.h, simdata.c = sim1.aug.c.cm, simdata.h = sim1.aug.h.cm)
d2 <- RMSE(true.h = md1.true.h, simdata.c = sim1.aug.c.mc, simdata.h = sim1.aug.h.mc)
e2 <- RMSE(true.h = md1.true.h, simdata.c = sim1.aug.c.mm, simdata.h = sim1.aug.h.mm)

a3 <- RE(simdata.c = sim1.wt.c, simdata.h = sim1.wt.h)
b3 <- RE(simdata.c = sim1.aug.c.cc, simdata.h = sim1.aug.h.cc)
c3 <- RE(simdata.c = sim1.aug.c.cm, simdata.h = sim1.aug.h.cm)
d3 <- RE(simdata.c = sim1.aug.c.mc, simdata.h = sim1.aug.h.mc)
e3 <- RE(simdata.c = sim1.aug.c.mm, simdata.h = sim1.aug.h.mm)

a4 <- CP(simdata.c = sim1.wt.c, simdata.h = sim1.wt.h)
b4 <- CP(simdata.c = sim1.aug.c.cc, simdata.h = sim1.aug.h.cc)
c4 <- CP(simdata.c = sim1.aug.c.cm, simdata.h = sim1.aug.h.cm)
d4 <- CP(simdata.c = sim1.aug.c.mc, simdata.h = sim1.aug.h.mc)
e4 <- CP(simdata.c = sim1.aug.c.mm, simdata.h = sim1.aug.h.mm)

######### Hajek like (Weighted) estimator
cons <- cbind(ac, a1$Rbias.avg.c, a2$RMSE.c, a3$RE.c, a4$CP.c)
hete <- cbind(ah, as.numeric(a1$Rbias.avg.h), a2$RMSE.h, a3$RE.h, a4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: both models are correctly specified
cons <- cbind(bc, b1$Rbias.avg.c, b2$RMSE.c, b3$RE.c, b4$CP.c)
hete <- cbind(bh, as.numeric(b1$Rbias.avg.h), b2$RMSE.h, b3$RE.h, b4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: only OR model is misspecified
cons <- cbind(cc, c1$Rbias.avg.c, c2$RMSE.c, c3$RE.c, c4$CP.c)
hete <- cbind(ch, as.numeric(c1$Rbias.avg.h), c2$RMSE.h, c3$RE.h, c4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: only PS model is misspecified
cons <- cbind(dc, d1$Rbias.avg.c, d2$RMSE.c, d3$RE.c, d4$CP.c)
hete <- cbind(dh, as.numeric(d1$Rbias.avg.h), d2$RMSE.h, d3$RE.h, d4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: both models are misspecified
cons <- cbind(ec, e1$Rbias.avg.c, e2$RMSE.c, e3$RE.c, e4$CP.c)
hete <- cbind(eh, as.numeric(e1$Rbias.avg.h), e2$RMSE.h, e3$RE.h, e4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

###### For table used in paper
### Table part I
cons <- cbind(ac, a1$Rbias.avg.c, a2$RMSE.c, a3$RE.c, a4$CP.c,
              bc, b1$Rbias.avg.c, b2$RMSE.c, b3$RE.c, b4$CP.c,
              cc, c1$Rbias.avg.c, c2$RMSE.c, c3$RE.c, c4$CP.c)
rownames(cons) <- rows
xtable(cons, digits = 2)

### Table part II
cons <- cbind(dc, d1$Rbias.avg.c, d2$RMSE.c, d3$RE.c, d4$CP.c,
              ec, e1$Rbias.avg.c, e2$RMSE.c, e3$RE.c, e4$CP.c)
rownames(cons) <- rows
xtable(cons, digits = 2)

### Table part III
cons <- cbind(ah, as.numeric(a1$Rbias.avg.h), a2$RMSE.h, a3$RE.h, a4$CP.h,
              bh, as.numeric(b1$Rbias.avg.h), b2$RMSE.h, b3$RE.h, b4$CP.h,
              ch, as.numeric(c1$Rbias.avg.h), c2$RMSE.h, c3$RE.h, c4$CP.h)
rownames(cons) <- rows
xtable(cons, digits = 2)

### Table part IV
cons <- cbind(dh, as.numeric(d1$Rbias.avg.h), d2$RMSE.h, d3$RE.h, d4$CP.h,
              eh, as.numeric(e1$Rbias.avg.h), e2$RMSE.h, e3$RE.h, e4$CP.h)
rownames(cons) <- rows
xtable(cons, digits = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Model 2
sim2.ratio
sim2.p

######## ESS
ess <- ESS(weight.data.reps2)
xtable(ess)

######### PE, ARBias, RMSE, RE, CP
ac <- PE(sim2.wt.c)
bc <- PE(sim2.aug.c.cc)
cc <- PE(sim2.aug.c.cm)
dc <- PE(sim2.aug.c.mc)
ec <- PE(sim2.aug.c.mm)

ah <- PE(sim2.wt.h)
bh <- PE(sim2.aug.h.cc)
ch <- PE(sim2.aug.h.cm)
dh <- PE(sim2.aug.h.mc)
eh <- PE(sim2.aug.h.mm)

a1 <- avg.rbias(true.h = md2.true.h, simdata.c = sim2.wt.c, simdata.h = sim2.wt.h)
b1 <- avg.rbias(true.h = md2.true.h, simdata.c = sim2.aug.c.cc, simdata.h = sim2.aug.h.cc)
c1 <- avg.rbias(true.h = md2.true.h, simdata.c = sim2.aug.c.cm, simdata.h = sim2.aug.h.cm)
d1 <- avg.rbias(true.h = md2.true.h, simdata.c = sim2.aug.c.mc, simdata.h = sim2.aug.h.mc)
e1 <- avg.rbias(true.h = md2.true.h, simdata.c = sim2.aug.c.mm, simdata.h = sim2.aug.h.mm)

a2 <- RMSE(true.h = md2.true.h, simdata.c = sim2.wt.c, simdata.h = sim2.wt.h)
b2 <- RMSE(true.h = md2.true.h, simdata.c = sim2.aug.c.cc, simdata.h = sim2.aug.h.cc)
c2 <- RMSE(true.h = md2.true.h, simdata.c = sim2.aug.c.cm, simdata.h = sim2.aug.h.cm)
d2 <- RMSE(true.h = md2.true.h, simdata.c = sim2.aug.c.mc, simdata.h = sim2.aug.h.mc)
e2 <- RMSE(true.h = md2.true.h, simdata.c = sim2.aug.c.mm, simdata.h = sim2.aug.h.mm)

a3 <- RE(simdata.c = sim2.wt.c, simdata.h = sim2.wt.h)
b3 <- RE(simdata.c = sim2.aug.c.cc, simdata.h = sim2.aug.h.cc)
c3 <- RE(simdata.c = sim2.aug.c.cm, simdata.h = sim2.aug.h.cm)
d3 <- RE(simdata.c = sim2.aug.c.mc, simdata.h = sim2.aug.h.mc)
e3 <- RE(simdata.c = sim2.aug.c.mm, simdata.h = sim2.aug.h.mm)

a4 <- CP(simdata.c = sim2.wt.c, simdata.h = sim2.wt.h)
b4 <- CP(simdata.c = sim2.aug.c.cc, simdata.h = sim2.aug.h.cc)
c4 <- CP(simdata.c = sim2.aug.c.cm, simdata.h = sim2.aug.h.cm)
d4 <- CP(simdata.c = sim2.aug.c.mc, simdata.h = sim2.aug.h.mc)
e4 <- CP(simdata.c = sim2.aug.c.mm, simdata.h = sim2.aug.h.mm)

######### Hajek like (Weighted) estimator
cons <- cbind(ac, a1$Rbias.avg.c, a2$RMSE.c, a3$RE.c, a4$CP.c)
hete <- cbind(ah, as.numeric(a1$Rbias.avg.h), a2$RMSE.h, a3$RE.h, a4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: both models are correctly specified
cons <- cbind(bc, b1$Rbias.avg.c, b2$RMSE.c, b3$RE.c, b4$CP.c)
hete <- cbind(bh, as.numeric(b1$Rbias.avg.h), b2$RMSE.h, b3$RE.h, b4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: only OR model is misspecified
cons <- cbind(cc, c1$Rbias.avg.c, c2$RMSE.c, c3$RE.c, c4$CP.c)
hete <- cbind(ch, as.numeric(c1$Rbias.avg.h), c2$RMSE.h, c3$RE.h, c4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: only PS model is misspecified
cons <- cbind(dc, d1$Rbias.avg.c, d2$RMSE.c, d3$RE.c, d4$CP.c)
hete <- cbind(dh, as.numeric(d1$Rbias.avg.h), d2$RMSE.h, d3$RE.h, d4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: both models are misspecified
cons <- cbind(ec, e1$Rbias.avg.c, e2$RMSE.c, e3$RE.c, e4$CP.c)
hete <- cbind(eh, as.numeric(e1$Rbias.avg.h), e2$RMSE.h, e3$RE.h, e4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

###### For table used in paper
### Table part I
cons <- cbind(ac, a1$Rbias.avg.c, a2$RMSE.c, a3$RE.c, a4$CP.c,
              bc, b1$Rbias.avg.c, b2$RMSE.c, b3$RE.c, b4$CP.c,
              cc, c1$Rbias.avg.c, c2$RMSE.c, c3$RE.c, c4$CP.c)
rownames(cons) <- rows
xtable(cons, digits = 2)

### Table part II
cons <- cbind(dc, d1$Rbias.avg.c, d2$RMSE.c, d3$RE.c, d4$CP.c,
              ec, e1$Rbias.avg.c, e2$RMSE.c, e3$RE.c, e4$CP.c)
rownames(cons) <- rows
xtable(cons, digits = 2)

### Table part III
cons <- cbind(ah, as.numeric(a1$Rbias.avg.h), a2$RMSE.h, a3$RE.h, a4$CP.h,
              bh, as.numeric(b1$Rbias.avg.h), b2$RMSE.h, b3$RE.h, b4$CP.h,
              ch, as.numeric(c1$Rbias.avg.h), c2$RMSE.h, c3$RE.h, c4$CP.h)
rownames(cons) <- rows
xtable(cons, digits = 2)

### Table part IV
cons <- cbind(dh, as.numeric(d1$Rbias.avg.h), d2$RMSE.h, d3$RE.h, d4$CP.h,
              eh, as.numeric(e1$Rbias.avg.h), e2$RMSE.h, e3$RE.h, e4$CP.h)
rownames(cons) <- rows
xtable(cons, digits = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model 3
sim3.ratio
sim3.p

### ESS
ess <- ESS(weight.data.reps3)
xtable(ess)

######### PE, ARBias, RMSE, RE, CP
ac <- PE(sim3.wt.c)
bc <- PE(sim3.aug.c.cc)
cc <- PE(sim3.aug.c.cm)
dc <- PE(sim3.aug.c.mc)
ec <- PE(sim3.aug.c.mm)

ah <- PE(sim3.wt.h)
bh <- PE(sim3.aug.h.cc)
ch <- PE(sim3.aug.h.cm)
dh <- PE(sim3.aug.h.mc)
eh <- PE(sim3.aug.h.mm)

a1 <- avg.rbias(true.h = md3.true.h, simdata.c = sim3.wt.c, simdata.h = sim3.wt.h)
b1 <- avg.rbias(true.h = md3.true.h, simdata.c = sim3.aug.c.cc, simdata.h = sim3.aug.h.cc)
c1 <- avg.rbias(true.h = md3.true.h, simdata.c = sim3.aug.c.cm, simdata.h = sim3.aug.h.cm)
d1 <- avg.rbias(true.h = md3.true.h, simdata.c = sim3.aug.c.mc, simdata.h = sim3.aug.h.mc)
e1 <- avg.rbias(true.h = md3.true.h, simdata.c = sim3.aug.c.mm, simdata.h = sim3.aug.h.mm)

a2 <- RMSE(true.h = md3.true.h, simdata.c = sim3.wt.c, simdata.h = sim3.wt.h)
b2 <- RMSE(true.h = md3.true.h, simdata.c = sim3.aug.c.cc, simdata.h = sim3.aug.h.cc)
c2 <- RMSE(true.h = md3.true.h, simdata.c = sim3.aug.c.cm, simdata.h = sim3.aug.h.cm)
d2 <- RMSE(true.h = md3.true.h, simdata.c = sim3.aug.c.mc, simdata.h = sim3.aug.h.mc)
e2 <- RMSE(true.h = md3.true.h, simdata.c = sim3.aug.c.mm, simdata.h = sim3.aug.h.mm)

a3 <- RE(simdata.c = sim3.wt.c, simdata.h = sim3.wt.h)
b3 <- RE(simdata.c = sim3.aug.c.cc, simdata.h = sim3.aug.h.cc)
c3 <- RE(simdata.c = sim3.aug.c.cm, simdata.h = sim3.aug.h.cm)
d3 <- RE(simdata.c = sim3.aug.c.mc, simdata.h = sim3.aug.h.mc)
e3 <- RE(simdata.c = sim3.aug.c.mm, simdata.h = sim3.aug.h.mm)

a4 <- CP(simdata.c = sim3.wt.c, simdata.h = sim3.wt.h)
b4 <- CP(simdata.c = sim3.aug.c.cc, simdata.h = sim3.aug.h.cc)
c4 <- CP(simdata.c = sim3.aug.c.cm, simdata.h = sim3.aug.h.cm)
d4 <- CP(simdata.c = sim3.aug.c.mc, simdata.h = sim3.aug.h.mc)
e4 <- CP(simdata.c = sim3.aug.c.mm, simdata.h = sim3.aug.h.mm)

######### Hajek like (Weighted) estimator
cons <- cbind(ac, a1$Rbias.avg.c, a2$RMSE.c, a3$RE.c, a4$CP.c)
hete <- cbind(ah, as.numeric(a1$Rbias.avg.h), a2$RMSE.h, a3$RE.h, a4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: both models are correctly specified
cons <- cbind(bc, b1$Rbias.avg.c, b2$RMSE.c, b3$RE.c, b4$CP.c)
hete <- cbind(bh, as.numeric(b1$Rbias.avg.h), b2$RMSE.h, b3$RE.h, b4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: only OR model is misspecified
cons <- cbind(cc, c1$Rbias.avg.c, c2$RMSE.c, c3$RE.c, c4$CP.c)
hete <- cbind(ch, as.numeric(c1$Rbias.avg.h), c2$RMSE.h, c3$RE.h, c4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: only PS model is misspecified
cons <- cbind(dc, d1$Rbias.avg.c, d2$RMSE.c, d3$RE.c, d4$CP.c)
hete <- cbind(dh, as.numeric(d1$Rbias.avg.h), d2$RMSE.h, d3$RE.h, d4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: both models are misspecified
cons <- cbind(ec, e1$Rbias.avg.c, e2$RMSE.c, e3$RE.c, e4$CP.c)
hete <- cbind(eh, as.numeric(e1$Rbias.avg.h), e2$RMSE.h, e3$RE.h, e4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

###### For table used in paper
### Table part I
cons <- cbind(ac, a1$Rbias.avg.c, a2$RMSE.c, a3$RE.c, a4$CP.c,
              bc, b1$Rbias.avg.c, b2$RMSE.c, b3$RE.c, b4$CP.c,
              cc, c1$Rbias.avg.c, c2$RMSE.c, c3$RE.c, c4$CP.c)
rownames(cons) <- rows
xtable(cons, digits = 2)

### Table part II
cons <- cbind(dc, d1$Rbias.avg.c, d2$RMSE.c, d3$RE.c, d4$CP.c,
              ec, e1$Rbias.avg.c, e2$RMSE.c, e3$RE.c, e4$CP.c)
rownames(cons) <- rows
xtable(cons, digits = 2)

### Table part III
cons <- cbind(ah, as.numeric(a1$Rbias.avg.h), a2$RMSE.h, a3$RE.h, a4$CP.h,
              bh, as.numeric(b1$Rbias.avg.h), b2$RMSE.h, b3$RE.h, b4$CP.h,
              ch, as.numeric(c1$Rbias.avg.h), c2$RMSE.h, c3$RE.h, c4$CP.h)
rownames(cons) <- rows
xtable(cons, digits = 2)

### Table part IV
cons <- cbind(dh, as.numeric(d1$Rbias.avg.h), d2$RMSE.h, d3$RE.h, d4$CP.h,
              eh, as.numeric(e1$Rbias.avg.h), e2$RMSE.h, e3$RE.h, e4$CP.h)
rownames(cons) <- rows
xtable(cons, digits = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model 4

sim4.ratio
sim4.p

######### ESS
ess <- ESS(weight.data.reps4)
xtable(ess)

######### PE, ARBias, RMSE, RE, CP
ac <- PE(sim4.wt.c)
bc <- PE(sim4.aug.c.cc)
cc <- PE(sim4.aug.c.cm)
dc <- PE(sim4.aug.c.mc)
ec <- PE(sim4.aug.c.mm)

ah <- PE(sim4.wt.h)
bh <- PE(sim4.aug.h.cc)
ch <- PE(sim4.aug.h.cm)
dh <- PE(sim4.aug.h.mc)
eh <- PE(sim4.aug.h.mm)

a1 <- avg.rbias(true.h = md4.true.h, simdata.c = sim4.wt.c, simdata.h = sim4.wt.h)
b1 <- avg.rbias(true.h = md4.true.h, simdata.c = sim4.aug.c.cc, simdata.h = sim4.aug.h.cc)
c1 <- avg.rbias(true.h = md4.true.h, simdata.c = sim4.aug.c.cm, simdata.h = sim4.aug.h.cm)
d1 <- avg.rbias(true.h = md4.true.h, simdata.c = sim4.aug.c.mc, simdata.h = sim4.aug.h.mc)
e1 <- avg.rbias(true.h = md4.true.h, simdata.c = sim4.aug.c.mm, simdata.h = sim4.aug.h.mm)

a2 <- RMSE(true.h = md4.true.h, simdata.c = sim4.wt.c, simdata.h = sim4.wt.h)
b2 <- RMSE(true.h = md4.true.h, simdata.c = sim4.aug.c.cc, simdata.h = sim4.aug.h.cc)
c2 <- RMSE(true.h = md4.true.h, simdata.c = sim4.aug.c.cm, simdata.h = sim4.aug.h.cm)
d2 <- RMSE(true.h = md4.true.h, simdata.c = sim4.aug.c.mc, simdata.h = sim4.aug.h.mc)
e2 <- RMSE(true.h = md4.true.h, simdata.c = sim4.aug.c.mm, simdata.h = sim4.aug.h.mm)

a3 <- RE(simdata.c = sim4.wt.c, simdata.h = sim4.wt.h)
b3 <- RE(simdata.c = sim4.aug.c.cc, simdata.h = sim4.aug.h.cc)
c3 <- RE(simdata.c = sim4.aug.c.cm, simdata.h = sim4.aug.h.cm)
d3 <- RE(simdata.c = sim4.aug.c.mc, simdata.h = sim4.aug.h.mc)
e3 <- RE(simdata.c = sim4.aug.c.mm, simdata.h = sim4.aug.h.mm)

a4 <- CP(simdata.c = sim4.wt.c, simdata.h = sim4.wt.h)
b4 <- CP(simdata.c = sim4.aug.c.cc, simdata.h = sim4.aug.h.cc)
c4 <- CP(simdata.c = sim4.aug.c.cm, simdata.h = sim4.aug.h.cm)
d4 <- CP(simdata.c = sim4.aug.c.mc, simdata.h = sim4.aug.h.mc)
e4 <- CP(simdata.c = sim4.aug.c.mm, simdata.h = sim4.aug.h.mm)

######### Hajek like (Weighted) estimator
cons <- cbind(ac, a1$Rbias.avg.c, a2$RMSE.c, a3$RE.c, a4$CP.c)
hete <- cbind(ah, as.numeric(a1$Rbias.avg.h), a2$RMSE.h, a3$RE.h, a4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: both models are correctly specified
cons <- cbind(bc, b1$Rbias.avg.c, b2$RMSE.c, b3$RE.c, b4$CP.c)
hete <- cbind(bh, as.numeric(b1$Rbias.avg.h), b2$RMSE.h, b3$RE.h, b4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: only OR model is misspecified
cons <- cbind(cc, c1$Rbias.avg.c, c2$RMSE.c, c3$RE.c, c4$CP.c)
hete <- cbind(ch, as.numeric(c1$Rbias.avg.h), c2$RMSE.h, c3$RE.h, c4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: only PS model is misspecified
cons <- cbind(dc, d1$Rbias.avg.c, d2$RMSE.c, d3$RE.c, d4$CP.c)
hete <- cbind(dh, as.numeric(d1$Rbias.avg.h), d2$RMSE.h, d3$RE.h, d4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: both models are misspecified
cons <- cbind(ec, e1$Rbias.avg.c, e2$RMSE.c, e3$RE.c, e4$CP.c)
hete <- cbind(eh, as.numeric(e1$Rbias.avg.h), e2$RMSE.h, e3$RE.h, e4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

###### For table used in paper
### Table part I
cons <- cbind(ac, a1$Rbias.avg.c, a2$RMSE.c, a3$RE.c, a4$CP.c,
              bc, b1$Rbias.avg.c, b2$RMSE.c, b3$RE.c, b4$CP.c,
              cc, c1$Rbias.avg.c, c2$RMSE.c, c3$RE.c, c4$CP.c)
rownames(cons) <- rows
xtable(cons, digits = 2)

### Table part II
cons <- cbind(dc, d1$Rbias.avg.c, d2$RMSE.c, d3$RE.c, d4$CP.c,
              ec, e1$Rbias.avg.c, e2$RMSE.c, e3$RE.c, e4$CP.c)
rownames(cons) <- rows
xtable(cons, digits = 2)

### Table part III
cons <- cbind(ah, as.numeric(a1$Rbias.avg.h), a2$RMSE.h, a3$RE.h, a4$CP.h,
              bh, as.numeric(b1$Rbias.avg.h), b2$RMSE.h, b3$RE.h, b4$CP.h,
              ch, as.numeric(c1$Rbias.avg.h), c2$RMSE.h, c3$RE.h, c4$CP.h)
rownames(cons) <- rows
xtable(cons, digits = 2)

### Table part IV
cons <- cbind(dh, as.numeric(d1$Rbias.avg.h), d2$RMSE.h, d3$RE.h, d4$CP.h,
              eh, as.numeric(e1$Rbias.avg.h), e2$RMSE.h, e3$RE.h, e4$CP.h)
rownames(cons) <- rows
xtable(cons, digits = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model 5

sim5.ratio
sim5.p

######### ESS
ess <- ESS(weight.data.reps5)
xtable(ess)

######### PE, ARBias, RMSE, RE, CP
ac <- PE(sim5.wt.c)
bc <- PE(sim5.aug.c.cc)
cc <- PE(sim5.aug.c.cm)
dc <- PE(sim5.aug.c.mc)
ec <- PE(sim5.aug.c.mm)

ah <- PE(sim5.wt.h)
bh <- PE(sim5.aug.h.cc)
ch <- PE(sim5.aug.h.cm)
dh <- PE(sim5.aug.h.mc)
eh <- PE(sim5.aug.h.mm)

a1 <- avg.rbias(true.h = md5.true.h, simdata.c = sim5.wt.c, simdata.h = sim5.wt.h)
b1 <- avg.rbias(true.h = md5.true.h, simdata.c = sim5.aug.c.cc, simdata.h = sim5.aug.h.cc)
c1 <- avg.rbias(true.h = md5.true.h, simdata.c = sim5.aug.c.cm, simdata.h = sim5.aug.h.cm)
d1 <- avg.rbias(true.h = md5.true.h, simdata.c = sim5.aug.c.mc, simdata.h = sim5.aug.h.mc)
e1 <- avg.rbias(true.h = md5.true.h, simdata.c = sim5.aug.c.mm, simdata.h = sim5.aug.h.mm)

a2 <- RMSE(true.h = md5.true.h, simdata.c = sim5.wt.c, simdata.h = sim5.wt.h)
b2 <- RMSE(true.h = md5.true.h, simdata.c = sim5.aug.c.cc, simdata.h = sim5.aug.h.cc)
c2 <- RMSE(true.h = md5.true.h, simdata.c = sim5.aug.c.cm, simdata.h = sim5.aug.h.cm)
d2 <- RMSE(true.h = md5.true.h, simdata.c = sim5.aug.c.mc, simdata.h = sim5.aug.h.mc)
e2 <- RMSE(true.h = md5.true.h, simdata.c = sim5.aug.c.mm, simdata.h = sim5.aug.h.mm)

a3 <- RE(simdata.c = sim5.wt.c, simdata.h = sim5.wt.h)
b3 <- RE(simdata.c = sim5.aug.c.cc, simdata.h = sim5.aug.h.cc)
c3 <- RE(simdata.c = sim5.aug.c.cm, simdata.h = sim5.aug.h.cm)
d3 <- RE(simdata.c = sim5.aug.c.mc, simdata.h = sim5.aug.h.mc)
e3 <- RE(simdata.c = sim5.aug.c.mm, simdata.h = sim5.aug.h.mm)

a4 <- CP(simdata.c = sim5.wt.c, simdata.h = sim5.wt.h)
b4 <- CP(simdata.c = sim5.aug.c.cc, simdata.h = sim5.aug.h.cc)
c4 <- CP(simdata.c = sim5.aug.c.cm, simdata.h = sim5.aug.h.cm)
d4 <- CP(simdata.c = sim5.aug.c.mc, simdata.h = sim5.aug.h.mc)
e4 <- CP(simdata.c = sim5.aug.c.mm, simdata.h = sim5.aug.h.mm)

######### Hajek like (Weighted) estimator
cons <- cbind(ac, a1$Rbias.avg.c, a2$RMSE.c, a3$RE.c, a4$CP.c)
hete <- cbind(ah, as.numeric(a1$Rbias.avg.h), a2$RMSE.h, a3$RE.h, a4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: both models are correctly specified
cons <- cbind(bc, b1$Rbias.avg.c, b2$RMSE.c, b3$RE.c, b4$CP.c)
hete <- cbind(bh, as.numeric(b1$Rbias.avg.h), b2$RMSE.h, b3$RE.h, b4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: only OR model is misspecified
cons <- cbind(cc, c1$Rbias.avg.c, c2$RMSE.c, c3$RE.c, c4$CP.c)
hete <- cbind(ch, as.numeric(c1$Rbias.avg.h), c2$RMSE.h, c3$RE.h, c4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: only PS model is misspecified
cons <- cbind(dc, d1$Rbias.avg.c, d2$RMSE.c, d3$RE.c, d4$CP.c)
hete <- cbind(dh, as.numeric(d1$Rbias.avg.h), d2$RMSE.h, d3$RE.h, d4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

######## Augmented estimator: both models are misspecified
cons <- cbind(ec, e1$Rbias.avg.c, e2$RMSE.c, e3$RE.c, e4$CP.c)
hete <- cbind(eh, as.numeric(e1$Rbias.avg.h), e2$RMSE.h, e3$RE.h, e4$CP.h)

colnames(cons) <- cols
rownames(cons) <- rows
colnames(hete) <- cols
rownames(hete) <- rows

xtable(cons, digits = 2)
xtable(hete, digits = 2)

###### For table used in paper
### Table part I
cons <- cbind(ac, a1$Rbias.avg.c, a2$RMSE.c, a3$RE.c, a4$CP.c,
              bc, b1$Rbias.avg.c, b2$RMSE.c, b3$RE.c, b4$CP.c,
              cc, c1$Rbias.avg.c, c2$RMSE.c, c3$RE.c, c4$CP.c)
rownames(cons) <- rows
xtable(cons, digits = 2)

### Table part II
cons <- cbind(dc, d1$Rbias.avg.c, d2$RMSE.c, d3$RE.c, d4$CP.c,
              ec, e1$Rbias.avg.c, e2$RMSE.c, e3$RE.c, e4$CP.c)
rownames(cons) <- rows
xtable(cons, digits = 2)

### Table part III
cons <- cbind(ah, as.numeric(a1$Rbias.avg.h), a2$RMSE.h, a3$RE.h, a4$CP.h,
              bh, as.numeric(b1$Rbias.avg.h), b2$RMSE.h, b3$RE.h, b4$CP.h,
              ch, as.numeric(c1$Rbias.avg.h), c2$RMSE.h, c3$RE.h, c4$CP.h)
rownames(cons) <- rows
xtable(cons, digits = 2)

### Table part IV
cons <- cbind(dh, as.numeric(d1$Rbias.avg.h), d2$RMSE.h, d3$RE.h, d4$CP.h,
              eh, as.numeric(e1$Rbias.avg.h), e2$RMSE.h, e3$RE.h, e4$CP.h)
rownames(cons) <- rows
xtable(cons, digits = 2)
