### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### XYplots for presentation

### by Yi Liu
### Jan 8, 2022
                      
library(lattice)

source("SumStat_func.R")
source("ESS_func.R")

# --- Plotting path
here <- getwd()
path <- paste(here, "/Results", sep="")

# Load the data
load("md1_sims.RData")
load("md2_sims.RData")
load("md3_sims.RData")
load("md4_sims.RData")
load("md5_sims.RData")
load("md_true.RData")

#### Point Estimation ####
p11c <- PE(sim1.wt.c)
p12c <- PE(sim1.aug.c.cc)
p13c <- PE(sim1.aug.c.cm)
p14c <- PE(sim1.aug.c.mc)
p15c <- PE(sim1.aug.c.mm)

p11h <- PE(sim1.wt.h)
p12h <- PE(sim1.aug.h.cc)
p13h <- PE(sim1.aug.h.cm)
p14h <- PE(sim1.aug.h.mc)
p15h <- PE(sim1.aug.h.mm)

p21c <- PE(sim3.wt.c)
p22c <- PE(sim3.aug.c.cc)
p23c <- PE(sim3.aug.c.cm)
p24c <- PE(sim3.aug.c.mc)
p25c <- PE(sim3.aug.c.mm)

p21h <- PE(sim3.wt.h)
p22h <- PE(sim3.aug.h.cc)
p23h <- PE(sim3.aug.h.cm)
p24h <- PE(sim3.aug.h.mc)
p25h <- PE(sim3.aug.h.mm)

p31c <- PE(sim5.wt.c)
p32c <- PE(sim5.aug.c.cc)
p33c <- PE(sim5.aug.c.cm)
p34c <- PE(sim5.aug.c.mc)
p35c <- PE(sim5.aug.c.mm)

p31h <- PE(sim5.wt.h)
p32h <- PE(sim5.aug.h.cc)
p33h <- PE(sim5.aug.h.cm)
p34h <- PE(sim5.aug.h.mc)
p35h <- PE(sim5.aug.h.mm)

XY.Est.c <- data.frame(Est = rep(0, 135), 
                       Model = c(rep("p = 10.05%", 45), rep("p = 45.94%", 45), rep("p = 89.18%", 45)),
                       Case = rep(c(rep("A", 9), rep("B", 9), rep("C", 9), rep("D", 9), rep("E", 9)), 3),
                       Estimand = rep(c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", 
                                        "ATO", "ATM", "ATEN", "ATC", "ATT"), 15))
XY.Est.h <- data.frame(Est = rep(0, 135), 
                       Model = c(rep("p = 10.05%", 45), rep("p = 45.94%", 45), rep("p = 89.18%", 45)),
                       Case = rep(c(rep("A", 9), rep("B", 9), rep("C", 9), rep("D", 9), rep("E", 9)), 3),
                       Estimand = rep(c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", 
                                        "ATO", "ATM", "ATEN", "ATC", "ATT"), 15))

XY.Est.c$Est <- c(p11c, p12c, p13c, p14c, p15c, 
                  p21c, p22c, p23c, p24c, p25c,
                  p31c, p32c, p33c, p34c, p35c)
XY.Est.h$Est <- c(p11h, p12h, p13h, p14h, p15h, 
                  p21h, p22h, p23h, p24h, p25h,
                  p31h, p32h, p33h, p34h, p35h)

XY.Est.c$Estimand <- factor(XY.Est.c$Estimand, 
                            levels = c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", "ATO", "ATM", "ATEN", "ATC", "ATT"))
XY.Est.h$Estimand <- factor(XY.Est.h$Estimand, 
                            levels = c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", "ATO", "ATM", "ATEN", "ATC", "ATT"))

colors <- c("blue", "dodgerblue", "green", "forestgreen", "brown", "red", "mediumslateblue", "goldenrod2", "gray27")
dotplot(Est~Case | Model, data = XY.Est.c, 
        ylab = "Point Estimation", ylim = c(-3,5),
        groups = Estimand, layout = c(3,1), 
        par.settings = list(superpose.symbol = list(pch = c(0,1,2,5,3,4,8,12,12), 
                                                    cex = 1.15, 
                                                    col = colors)),
        auto.key = list(space = "top", points = "TRUE", columns=3)) -> p.est.c
dotplot(Est~Case | Model, data = XY.Est.h, 
        ylab = "Point Estimation",
        groups = Estimand, layout = c(3,1),
        par.settings = list(superpose.symbol = list(pch = c(0,1,2,5,3,4,8,12,12), 
                                                    cex = 1.15, 
                                                    col = colors)),
        auto.key = list(space = "top", points = "TRUE", columns=3)) -> p.est.h

png(filename = paste(path, "/consXY.png", sep=""),
    width = 2200, height = 1800, res = 312)
p.est.c
dev.off()

png(filename = paste(path, "/heteXY.png", sep=""), 
    width = 2500, height = 1300, res = 312)
p.est.h
dev.off()

#### ARBias ####
b11 <- avg.rbias(true.h = md1.true.h, simdata.c = sim1.wt.c, simdata.h = sim1.wt.h)
b12 <- avg.rbias(true.h = md1.true.h, simdata.c = sim1.aug.c.cc, simdata.h = sim1.aug.h.cc)
b13 <- avg.rbias(true.h = md1.true.h, simdata.c = sim1.aug.c.cm, simdata.h = sim1.aug.h.cm)
b14 <- avg.rbias(true.h = md1.true.h, simdata.c = sim1.aug.c.mc, simdata.h = sim1.aug.h.mc)
b15 <- avg.rbias(true.h = md1.true.h, simdata.c = sim1.aug.c.mm, simdata.h = sim1.aug.h.mm)

b21 <- avg.rbias(true.h = md3.true.h, simdata.c = sim3.wt.c, simdata.h = sim3.wt.h)
b22 <- avg.rbias(true.h = md3.true.h, simdata.c = sim3.aug.c.cc, simdata.h = sim3.aug.h.cc)
b23 <- avg.rbias(true.h = md3.true.h, simdata.c = sim3.aug.c.cm, simdata.h = sim3.aug.h.cm)
b24 <- avg.rbias(true.h = md3.true.h, simdata.c = sim3.aug.c.mc, simdata.h = sim3.aug.h.mc)
b25 <- avg.rbias(true.h = md3.true.h, simdata.c = sim3.aug.c.mm, simdata.h = sim3.aug.h.mm)

b31 <- avg.rbias(true.h = md5.true.h, simdata.c = sim5.wt.c, simdata.h = sim5.wt.h)
b32 <- avg.rbias(true.h = md5.true.h, simdata.c = sim5.aug.c.cc, simdata.h = sim5.aug.h.cc)
b33 <- avg.rbias(true.h = md5.true.h, simdata.c = sim5.aug.c.cm, simdata.h = sim5.aug.h.cm)
b34 <- avg.rbias(true.h = md5.true.h, simdata.c = sim5.aug.c.mc, simdata.h = sim5.aug.h.mc)
b35 <- avg.rbias(true.h = md5.true.h, simdata.c = sim5.aug.c.mm, simdata.h = sim5.aug.h.mm)

XY.bias <- data.frame(Bias = rep(0, 270), 
                      Model = c(rep("p = 10.05%", 90), rep("p = 45.94%", 90), rep("p = 89.18%", 90)),
                      Case = rep(c(rep("A", 18), rep("B", 18), rep("C", 18), rep("D", 18), rep("E", 18)), 3),
                      Effect = rep(c(rep("Constant", 9), rep("Heterogeneous", 9)), 15),
                      Estimand = rep(c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", 
                                       "ATO", "ATM", "ATEN", "ATC", "ATT"), 30))

XY.bias$Bias <- c(b11$Rbias.avg.c, as.numeric(b11$Rbias.avg.h),
                  b12$Rbias.avg.c, as.numeric(b12$Rbias.avg.h),
                  b13$Rbias.avg.c, as.numeric(b13$Rbias.avg.h),
                  b14$Rbias.avg.c, as.numeric(b14$Rbias.avg.h),
                  b15$Rbias.avg.c, as.numeric(b15$Rbias.avg.h),
                  
                  b21$Rbias.avg.c, as.numeric(b21$Rbias.avg.h),
                  b22$Rbias.avg.c, as.numeric(b22$Rbias.avg.h),
                  b23$Rbias.avg.c, as.numeric(b23$Rbias.avg.h),
                  b24$Rbias.avg.c, as.numeric(b24$Rbias.avg.h),
                  b25$Rbias.avg.c, as.numeric(b25$Rbias.avg.h),
                  
                  b31$Rbias.avg.c, as.numeric(b31$Rbias.avg.h),
                  b32$Rbias.avg.c, as.numeric(b32$Rbias.avg.h),
                  b33$Rbias.avg.c, as.numeric(b33$Rbias.avg.h),
                  b34$Rbias.avg.c, as.numeric(b34$Rbias.avg.h),
                  b35$Rbias.avg.c, as.numeric(b35$Rbias.avg.h))

XY.bias$Estimand <- factor(XY.bias$Estimand, 
                           levels = c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", "ATO", "ATM", "ATEN", "ATC", "ATT"))

colors <- c("blue", "dodgerblue", "green", "forestgreen",  "red", "brown","mediumslateblue", "goldenrod2", "gray27")
dotplot(Bias~Case | Model*Effect, data = XY.bias, 
        ylab = "Absolute Relative Percent Bias (%)",
        groups = Estimand, ylim = c(-5,55),
        par.settings = list(superpose.symbol = list(pch = c(0,1,2,5,4,3,8,12,12), 
                                                    cex = 1.15, 
                                                    col = colors)),
        auto.key = list(space = "right", points = "TRUE"),
        panel = function(...) {
          panel.abline(h = 0, lty = "dotted", col = "gray48")
          panel.dotplot(...)
        }) -> p.bias

png(filename = paste(path, "/bias_XY.png", sep=""),
    width = 2500, height = 1800, res = 288)
p.bias
dev.off()

#### RMSE ####
rm11 <- RMSE(true.h = md1.true.h, simdata.c = sim1.wt.c, simdata.h = sim1.wt.h)
rm12 <- RMSE(true.h = md1.true.h, simdata.c = sim1.aug.c.cc, simdata.h = sim1.aug.h.cc)
rm13 <- RMSE(true.h = md1.true.h, simdata.c = sim1.aug.c.cm, simdata.h = sim1.aug.h.cm)
rm14 <- RMSE(true.h = md1.true.h, simdata.c = sim1.aug.c.mc, simdata.h = sim1.aug.h.mc)
rm15 <- RMSE(true.h = md1.true.h, simdata.c = sim1.aug.c.mm, simdata.h = sim1.aug.h.mm)

rm21 <- RMSE(true.h = md3.true.h, simdata.c = sim3.wt.c, simdata.h = sim3.wt.h)
rm22 <- RMSE(true.h = md3.true.h, simdata.c = sim3.aug.c.cc, simdata.h = sim3.aug.h.cc)
rm23 <- RMSE(true.h = md3.true.h, simdata.c = sim3.aug.c.cm, simdata.h = sim3.aug.h.cm)
rm24 <- RMSE(true.h = md3.true.h, simdata.c = sim3.aug.c.mc, simdata.h = sim3.aug.h.mc)
rm25 <- RMSE(true.h = md3.true.h, simdata.c = sim3.aug.c.mm, simdata.h = sim3.aug.h.mm)

rm31 <- RMSE(true.h = md5.true.h, simdata.c = sim5.wt.c, simdata.h = sim5.wt.h)
rm32 <- RMSE(true.h = md5.true.h, simdata.c = sim5.aug.c.cc, simdata.h = sim5.aug.h.cc)
rm33 <- RMSE(true.h = md5.true.h, simdata.c = sim5.aug.c.cm, simdata.h = sim5.aug.h.cm)
rm34 <- RMSE(true.h = md5.true.h, simdata.c = sim5.aug.c.mc, simdata.h = sim5.aug.h.mc)
rm35 <- RMSE(true.h = md5.true.h, simdata.c = sim5.aug.c.mm, simdata.h = sim5.aug.h.mm)

XY.RMSE <- data.frame(RMSE = rep(0, 270), 
                      Model = c(rep("p = 10.05%", 90), rep("p = 45.94%", 90), rep("p = 89.18%", 90)),
                      Case = rep(c(rep("A", 18), rep("B", 18), rep("C", 18), rep("D", 18), rep("E", 18)), 3),
                      Effect = rep(c(rep("Constant", 9), rep("Heterogeneous", 9)), 15),
                      Estimand = rep(c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", 
                                       "ATO", "ATM", "ATEN", "ATC", "ATT"), 30))

XY.RMSE$RMSE <- c(rm11$RMSE.c, rm11$RMSE.h,
                  rm12$RMSE.c, rm12$RMSE.h,
                  rm13$RMSE.c, rm13$RMSE.h,
                  rm14$RMSE.c, rm14$RMSE.h,
                  rm15$RMSE.c, rm15$RMSE.h,
                  
                  rm21$RMSE.c, rm21$RMSE.h,
                  rm22$RMSE.c, rm22$RMSE.h,
                  rm23$RMSE.c, rm23$RMSE.h,
                  rm24$RMSE.c, rm24$RMSE.h,
                  rm25$RMSE.c, rm25$RMSE.h,
                  
                  rm31$RMSE.c, rm31$RMSE.h,
                  rm32$RMSE.c, rm32$RMSE.h,
                  rm33$RMSE.c, rm33$RMSE.h,
                  rm34$RMSE.c, rm34$RMSE.h,
                  rm35$RMSE.c, rm35$RMSE.h)

XY.RMSE$Estimand <- factor(XY.RMSE$Estimand, 
                           levels = c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", "ATO", "ATM", "ATEN", "ATC", "ATT"))

colors <- c("blue", "dodgerblue", "green", "forestgreen", "red", "brown", "mediumslateblue", "goldenrod2", "gray27")
dotplot(RMSE~Case | Model*Effect, data = XY.RMSE, 
        ylab = "Rooted Mean Squared Error",
        groups = Estimand, 
        par.settings = list(superpose.symbol = list(pch = c(0,1,2,5,4,3,8,12,12),  
                                                    cex = 1.15, 
                                                    col = colors)),
        auto.key = list(space = "right", points = "TRUE"),
        panel = function(...) {
          panel.abline(h = 0, lty = "dotted", col = "gray48")
          panel.dotplot(...)
        }) -> p.rmse

png(filename = paste(path, "/rmse_XY.png", sep=""), 
    width = 2500, height = 1800, res = 288)
p.rmse
dev.off()

#### Relative Efficiency ####
r11 <- RE(simdata.c = sim1.wt.c, simdata.h = sim1.wt.h)
r12 <- RE(simdata.c = sim1.aug.c.cc, simdata.h = sim1.aug.h.cc)
r13 <- RE(simdata.c = sim1.aug.c.cm, simdata.h = sim1.aug.h.cm)
r14 <- RE(simdata.c = sim1.aug.c.mc, simdata.h = sim1.aug.h.mc)
r15 <- RE(simdata.c = sim1.aug.c.mm, simdata.h = sim1.aug.h.mm)

r21 <- RE(simdata.c = sim3.wt.c, simdata.h = sim3.wt.h)
r22 <- RE(simdata.c = sim3.aug.c.cc, simdata.h = sim3.aug.h.cc)
r23 <- RE(simdata.c = sim3.aug.c.cm, simdata.h = sim3.aug.h.cm)
r24 <- RE(simdata.c = sim3.aug.c.mc, simdata.h = sim3.aug.h.mc)
r25 <- RE(simdata.c = sim3.aug.c.mm, simdata.h = sim3.aug.h.mm)

r31 <- RE(simdata.c = sim5.wt.c, simdata.h = sim5.wt.h)
r32 <- RE(simdata.c = sim5.aug.c.cc, simdata.h = sim5.aug.h.cc)
r33 <- RE(simdata.c = sim5.aug.c.cm, simdata.h = sim5.aug.h.cm)
r34 <- RE(simdata.c = sim5.aug.c.mc, simdata.h = sim5.aug.h.mc)
r35 <- RE(simdata.c = sim5.aug.c.mm, simdata.h = sim5.aug.h.mm)

XY.RE <- data.frame(RE = rep(0, 270), 
                    Model = c(rep("p = 10.05%", 90), rep("p = 45.94%", 90), rep("p = 89.18%", 90)),
                    Case = rep(c(rep("A", 18), rep("B", 18), rep("C", 18), rep("D", 18), rep("E", 18)), 3),
                    Effect = rep(c(rep("Constant", 9), rep("Heterogeneous", 9)), 15),
                    Estimand = rep(c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", 
                                     "ATO", "ATM", "ATEN", "ATC", "ATT"), 30))

XY.RE$RE <- c(r11$RE.c, r11$RE.h, r12$RE.c, r12$RE.h,
              r13$RE.c, r13$RE.h, r14$RE.c, r14$RE.h, r15$RE.c, r15$RE.h,
              
              r21$RE.c, r21$RE.h, r22$RE.c, r22$RE.h,
              r23$RE.c, r23$RE.h, r24$RE.c, r24$RE.h, r25$RE.c, r25$RE.h,
              
              r31$RE.c, r31$RE.h, r32$RE.c, r32$RE.h,
              r33$RE.c, r33$RE.h, r34$RE.c, r34$RE.h, r35$RE.c, r35$RE.h)

XY.RE$Estimand <- factor(XY.RE$Estimand, 
                         levels = c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", "ATO", "ATM", "ATEN", "ATC", "ATT"))

colors <- c("blue", "dodgerblue", "green", "forestgreen", "red", "brown", "mediumslateblue", "goldenrod2", "gray27")
dotplot(RE~Case | Model*Effect, data = XY.RE, 
        ylab = "Relative Efficiency",
        groups = Estimand, 
        par.settings = list(superpose.symbol = list(pch = c(0,1,2,5,4,3,8,12,12), 
                                                    cex = 1.15, 
                                                    col = colors)),
        auto.key = list(space = "right", points = "TRUE"),
        panel = function(...) {
          panel.abline(h = 1, lty = "dotted", col = "gray48")
          panel.dotplot(...)
        }) -> p.re

png(filename = paste(path, "/re_XY.png", sep=""), 
    width = 2500, height = 1800, res = 288)
p.re
dev.off()

#### Coverage Rate ####
c11 <- CP(simdata.c = sim1.wt.c, simdata.h = sim1.wt.h)
c12 <- CP(simdata.c = sim1.aug.c.cc, simdata.h = sim1.aug.h.cc)
c13 <- CP(simdata.c = sim1.aug.c.cm, simdata.h = sim1.aug.h.cm)
c14 <- CP(simdata.c = sim1.aug.c.mc, simdata.h = sim1.aug.h.mc)
c15 <- CP(simdata.c = sim1.aug.c.mm, simdata.h = sim1.aug.h.mm)

c21 <- CP(simdata.c = sim3.wt.c, simdata.h = sim3.wt.h)
c22 <- CP(simdata.c = sim3.aug.c.cc, simdata.h = sim3.aug.h.cc)
c23 <- CP(simdata.c = sim3.aug.c.cm, simdata.h = sim3.aug.h.cm)
c24 <- CP(simdata.c = sim3.aug.c.mc, simdata.h = sim3.aug.h.mc)
c25 <- CP(simdata.c = sim3.aug.c.mm, simdata.h = sim3.aug.h.mm)

c31 <- CP(simdata.c = sim5.wt.c, simdata.h = sim5.wt.h)
c32 <- CP(simdata.c = sim5.aug.c.cc, simdata.h = sim5.aug.h.cc)
c33 <- CP(simdata.c = sim5.aug.c.cm, simdata.h = sim5.aug.h.cm)
c34 <- CP(simdata.c = sim5.aug.c.mc, simdata.h = sim5.aug.h.mc)
c35 <- CP(simdata.c = sim5.aug.c.mm, simdata.h = sim5.aug.h.mm)

XY.CP <- data.frame(CP = rep(0, 270), 
                    Model = c(rep("p = 10.05%", 90), rep("p = 45.94%", 90), rep("p = 89.18%", 90)),
                    Case = rep(c(rep("A", 18), rep("B", 18), rep("C", 18), rep("D", 18), rep("E", 18)), 3),
                    Effect = rep(c(rep("Constant", 9), rep("Heterogeneous", 9)), 15),
                    Estimand = rep(c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", 
                                     "ATO", "ATM", "ATEN", "ATC", "ATT"), 30))

XY.CP$CP <- c(c11$CP.c, c11$CP.h, c12$CP.c, c12$CP.h,
              c13$CP.c, c13$CP.h, c14$CP.c, c14$CP.h, c15$CP.c, c15$CP.h,
              
              c21$CP.c, c21$CP.h, c22$CP.c, c22$CP.h,
              c23$CP.c, c23$CP.h, c24$CP.c, c24$CP.h, c25$CP.c, c25$CP.h,
              
              c31$CP.c, c31$CP.h, c32$CP.c, c32$CP.h,
              c33$CP.c, c33$CP.h, c34$CP.c, c34$CP.h, c35$CP.c, c35$CP.h)

XY.CP$Estimand <- factor(XY.CP$Estimand, 
                         levels = c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", "ATO", "ATM", "ATEN", "ATC", "ATT"))

colors <- c("blue", "dodgerblue", "green", "forestgreen", "red", "brown", "mediumslateblue", "goldenrod2", "gray27")
dotplot(CP~Case | Model*Effect, data = XY.CP, 
        ylab = "Coverage Rate",
        groups = Estimand, 
        par.settings = list(superpose.symbol = list(pch = c(0,1,2,5,4,3,8,12,12), 
                                                    cex = 1.15, 
                                                    col = colors)),
        auto.key = list(space = "right", points = "TRUE"),
        panel = function(...) {
          panel.abline(h = 0.95, lty = "dotted", col = "gray48")
          panel.dotplot(...)
        }) -> p.cp

png(filename = paste(path, "/cp_XY.png", sep=""), 
    width = 2500, height = 1800, res = 288)
p.cp
dev.off()

