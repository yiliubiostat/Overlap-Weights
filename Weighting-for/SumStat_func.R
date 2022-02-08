### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

####### Summary statistics and plotting functions

### by Yi Liu
### Sept 29, 2021

library(PSweight)
library(ggplot2)

# Function for Point Estimation plot
PEst_plot <- function(sim.wt.c, sim.wt.h, sim.aug.c.cc, sim.aug.h.cc,
                      sim.aug.c.cm, sim.aug.h.cm, sim.aug.c.mc, sim.aug.h.mc,
                      sim.aug.c.mm, sim.aug.h.mm,
                      plot.title = "Point Estimation Distribution", plot.subtitle = "") {
  
  # Only allow inputting dataset generated in this simulation project
  # plot.title: because each plot has its own case
  
  # Hajek-like (Weighted)
  df1 <- data.frame(PE = sim.wt.c[,"IPW.c"],
                    method = "ATE", type = "Hajek-like", effect = "Constant Treatment Effect")
  df2 <- data.frame(PE = sim.wt.c[,"IPW.5.c"],
                    method = "ATE 0.05", type = "Hajek-like", effect = "Constant Treatment Effect")
  df3 <- data.frame(PE = sim.wt.c[,"IPW.10.c"],
                    method = "ATE 0.1", type = "Hajek-like", effect = "Constant Treatment Effect")
  df4 <- data.frame(PE = sim.wt.c[,"IPW.15.c"],
                    method = "ATE 0.15", type = "Hajek-like", effect = "Constant Treatment Effect")
  df5 <- data.frame(PE = sim.wt.c[,"ATO.c"],
                    method = "ATO", type = "Hajek-like", effect = "Constant Treatment Effect")
  df6 <- data.frame(PE = sim.wt.c[,"ATM.c"],
                    method = "ATM", type = "Hajek-like", effect = "Constant Treatment Effect")
  df7 <- data.frame(PE = sim.wt.c[,"ATEN.c"],
                    method = "ATEN", type = "Hajek-like", effect = "Constant Treatment Effect")
  df8 <- data.frame(PE = sim.wt.c[,"ATC.c"],
                    method = "ATC", type = "Hajek-like", effect = "Constant Treatment Effect")
  df9 <- data.frame(PE = sim.wt.c[,"ATT.c"],
                    method = "ATT", type = "Hajek-like", effect = "Constant Treatment Effect")
  
  PEst.wt.c <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)

  df1 <- data.frame(PE = sim.wt.h[,"IPW.h"],
                    method = "ATE", type = "Hajek-like", effect = "Heterogenous Treatment Effect")
  df2 <- data.frame(PE = sim.wt.h[,"IPW.5.h"],
                    method = "ATE 0.05", type = "Hajek-like", effect = "Heterogenous Treatment Effect")
  df3 <- data.frame(PE = sim.wt.h[,"IPW.10.h"],
                    method = "ATE 0.1", type = "Hajek-like", effect = "Heterogenous Treatment Effect")
  df4 <- data.frame(PE = sim.wt.h[,"IPW.15.h"],
                    method = "ATE 0.15", type = "Hajek-like", effect = "Heterogenous Treatment Effect")
  df5 <- data.frame(PE = sim.wt.h[,"ATO.h"],
                    method = "ATO", type = "Hajek-like", effect = "Heterogenous Treatment Effect")
  df6 <- data.frame(PE = sim.wt.h[,"ATM.h"],
                    method = "ATM", type = "Hajek-like", effect = "Heterogenous Treatment Effect")
  df7 <- data.frame(PE = sim.wt.h[,"ATEN.h"],
                    method = "ATEN", type = "Hajek-like", effect = "Heterogenous Treatment Effect")
  df8 <- data.frame(PE = sim.wt.h[,"ATC.h"],
                    method = "ATC", type = "Hajek-like", effect = "Heterogenous Treatment Effect")
  df9 <- data.frame(PE = sim.wt.h[,"ATT.h"],
                    method = "ATT", type = "Hajek-like", effect = "Heterogenous Treatment Effect")
  
  PEst.wt.h <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  # Augmented estimator case 1
  df1 <- data.frame(PE = sim.aug.c.cc[,"IPW.c"],
                    method = "ATE", type = "Aug Case 1", effect = "Constant Treatment Effect")
  df2 <- data.frame(PE = sim.aug.c.cc[,"IPW.5.c"],
                    method = "ATE 0.05", type = "Aug Case 1", effect = "Constant Treatment Effect")
  df3 <- data.frame(PE = sim.aug.c.cc[,"IPW.10.c"],
                    method = "ATE 0.1", type = "Aug Case 1", effect = "Constant Treatment Effect")
  df4 <- data.frame(PE = sim.aug.c.cc[,"IPW.15.c"],
                    method = "ATE 0.15", type = "Aug Case 1", effect = "Constant Treatment Effect")
  df5 <- data.frame(PE = sim.aug.c.cc[,"ATO.c"],
                    method = "ATO", type = "Aug Case 1", effect = "Constant Treatment Effect")
  df6 <- data.frame(PE = sim.aug.c.cc[,"ATM.c"],
                    method = "ATM", type = "Aug Case 1", effect = "Constant Treatment Effect")
  df7 <- data.frame(PE = sim.aug.c.cc[,"ATEN.c"],
                    method = "ATEN", type = "Aug Case 1", effect = "Constant Treatment Effect")
  df8 <- data.frame(PE = sim.aug.c.cc[,"ATC.c"],
                    method = "ATC", type = "Aug Case 1", effect = "Constant Treatment Effect")
  df9 <- data.frame(PE = sim.aug.c.cc[,"ATT.c"],
                    method = "ATT", type = "Aug Case 1", effect = "Constant Treatment Effect")
  
  PEst.aug.c.cc <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  df1 <- data.frame(PE = sim.aug.h.cc[,"IPW.h"],
                    method = "ATE", type = "Aug Case 1", effect = "Heterogenous Treatment Effect")
  df2 <- data.frame(PE = sim.aug.h.cc[,"IPW.5.h"],
                    method = "ATE 0.05", type = "Aug Case 1", effect = "Heterogenous Treatment Effect")
  df3 <- data.frame(PE = sim.aug.h.cc[,"IPW.10.h"],
                    method = "ATE 0.1", type = "Aug Case 1", effect = "Heterogenous Treatment Effect")
  df4 <- data.frame(PE = sim.aug.h.cc[,"IPW.15.h"],
                    method = "ATE 0.15", type = "Aug Case 1", effect = "Heterogenous Treatment Effect")
  df5 <- data.frame(PE = sim.aug.h.cc[,"ATO.h"],
                    method = "ATO", type = "Aug Case 1", effect = "Heterogenous Treatment Effect")
  df6 <- data.frame(PE = sim.aug.h.cc[,"ATM.h"],
                    method = "ATM", type = "Aug Case 1", effect = "Heterogenous Treatment Effect")
  df7 <- data.frame(PE = sim.aug.h.cc[,"ATEN.h"],
                    method = "ATEN", type = "Aug Case 1", effect = "Heterogenous Treatment Effect")
  df8 <- data.frame(PE = sim.aug.h.cc[,"ATC.h"],
                    method = "ATC", type = "Aug Case 1", effect = "Heterogenous Treatment Effect")
  df9 <- data.frame(PE = sim.aug.h.cc[,"ATT.h"],
                    method = "ATT", type = "Aug Case 1", effect = "Heterogenous Treatment Effect")
  
  PEst.aug.h.cc <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  # Augmented estimator case 2
  df1 <- data.frame(PE = sim.aug.c.cm[,"IPW.c"],
                    method = "ATE", type = "Aug Case 2", effect = "Constant Treatment Effect")
  df2 <- data.frame(PE = sim.aug.c.cm[,"IPW.5.c"],
                    method = "ATE 0.05", type = "Aug Case 2", effect = "Constant Treatment Effect")
  df3 <- data.frame(PE = sim.aug.c.cm[,"IPW.10.c"],
                    method = "ATE 0.1", type = "Aug Case 2", effect = "Constant Treatment Effect")
  df4 <- data.frame(PE = sim.aug.c.cm[,"IPW.15.c"],
                    method = "ATE 0.15", type = "Aug Case 2", effect = "Constant Treatment Effect")
  df5 <- data.frame(PE = sim.aug.c.cm[,"ATO.c"],
                    method = "ATO", type = "Aug Case 2", effect = "Constant Treatment Effect")
  df6 <- data.frame(PE = sim.aug.c.cm[,"ATM.c"],
                    method = "ATM", type = "Aug Case 2", effect = "Constant Treatment Effect")
  df7 <- data.frame(PE = sim.aug.c.cm[,"ATEN.c"],
                    method = "ATEN", type = "Aug Case 2", effect = "Constant Treatment Effect")
  df8 <- data.frame(PE = sim.aug.c.cm[,"ATC.c"],
                    method = "ATC", type = "Aug Case 2", effect = "Constant Treatment Effect")
  df9 <- data.frame(PE = sim.aug.c.cm[,"ATT.c"],
                    method = "ATT", type = "Aug Case 2", effect = "Constant Treatment Effect")
  
  PEst.aug.c.cm <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  df1 <- data.frame(PE = sim.aug.h.cm[,"IPW.h"],
                    method = "ATE", type = "Aug Case 2", effect = "Heterogenous Treatment Effect")
  df2 <- data.frame(PE = sim.aug.h.cm[,"IPW.5.h"],
                    method = "ATE 0.05", type = "Aug Case 2", effect = "Heterogenous Treatment Effect")
  df3 <- data.frame(PE = sim.aug.h.cm[,"IPW.10.h"],
                    method = "ATE 0.1", type = "Aug Case 2", effect = "Heterogenous Treatment Effect")
  df4 <- data.frame(PE = sim.aug.h.cm[,"IPW.15.h"],
                    method = "ATE 0.15", type = "Aug Case 2", effect = "Heterogenous Treatment Effect")
  df5 <- data.frame(PE = sim.aug.h.cm[,"ATO.h"],
                    method = "ATO", type = "Aug Case 2", effect = "Heterogenous Treatment Effect")
  df6 <- data.frame(PE = sim.aug.h.cm[,"ATM.h"],
                    method = "ATM", type = "Aug Case 2", effect = "Heterogenous Treatment Effect")
  df7 <- data.frame(PE = sim.aug.h.cm[,"ATEN.h"],
                    method = "ATEN", type = "Aug Case 2", effect = "Heterogenous Treatment Effect")
  df8 <- data.frame(PE = sim.aug.h.cm[,"ATC.h"],
                    method = "ATC", type = "Aug Case 2", effect = "Heterogenous Treatment Effect")
  df9 <- data.frame(PE = sim.aug.h.cm[,"ATT.h"],
                    method = "ATT", type = "Aug Case 2", effect = "Heterogenous Treatment Effect")
  
  PEst.aug.h.cm <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)

  # Augmented estimator case 3
  df1 <- data.frame(PE = sim.aug.c.mc[,"IPW.c"],
                    method = "ATE", type = "Aug Case 3", effect = "Constant Treatment Effect")
  df2 <- data.frame(PE = sim.aug.c.mc[,"IPW.5.c"],
                    method = "ATE 0.05", type = "Aug Case 3", effect = "Constant Treatment Effect")
  df3 <- data.frame(PE = sim.aug.c.mc[,"IPW.10.c"],
                    method = "ATE 0.1", type = "Aug Case 3", effect = "Constant Treatment Effect")
  df4 <- data.frame(PE = sim.aug.c.mc[,"IPW.15.c"],
                    method = "ATE 0.15", type = "Aug Case 3", effect = "Constant Treatment Effect")
  df5 <- data.frame(PE = sim.aug.c.mc[,"ATO.c"],
                    method = "ATO", type = "Aug Case 3", effect = "Constant Treatment Effect")
  df6 <- data.frame(PE = sim.aug.c.mc[,"ATM.c"],
                    method = "ATM", type = "Aug Case 3", effect = "Constant Treatment Effect")
  df7 <- data.frame(PE = sim.aug.c.mc[,"ATEN.c"],
                    method = "ATEN", type = "Aug Case 3", effect = "Constant Treatment Effect")
  df8 <- data.frame(PE = sim.aug.c.mc[,"ATC.c"],
                    method = "ATC", type = "Aug Case 3", effect = "Constant Treatment Effect")
  df9 <- data.frame(PE = sim.aug.c.mc[,"ATT.c"],
                    method = "ATT", type = "Aug Case 3", effect = "Constant Treatment Effect")
  
  PEst.aug.c.mc <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  df1 <- data.frame(PE = sim.aug.h.mc[,"IPW.h"],
                    method = "ATE", type = "Aug Case 3", effect = "Heterogenous Treatment Effect")
  df2 <- data.frame(PE = sim.aug.h.mc[,"IPW.5.h"],
                    method = "ATE 0.05", type = "Aug Case 3", effect = "Heterogenous Treatment Effect")
  df3 <- data.frame(PE = sim.aug.h.mc[,"IPW.10.h"],
                    method = "ATE 0.1", type = "Aug Case 3", effect = "Heterogenous Treatment Effect")
  df4 <- data.frame(PE = sim.aug.h.mc[,"IPW.15.h"],
                    method = "ATE 0.15", type = "Aug Case 3", effect = "Heterogenous Treatment Effect")
  df5 <- data.frame(PE = sim.aug.h.mc[,"ATO.h"],
                    method = "ATO", type = "Aug Case 3", effect = "Heterogenous Treatment Effect")
  df6 <- data.frame(PE = sim.aug.h.mc[,"ATM.h"],
                    method = "ATM", type = "Aug Case 3", effect = "Heterogenous Treatment Effect")
  df7 <- data.frame(PE = sim.aug.h.mc[,"ATEN.h"],
                    method = "ATEN", type = "Aug Case 3", effect = "Heterogenous Treatment Effect")
  df8 <- data.frame(PE = sim.aug.h.mc[,"ATC.h"],
                    method = "ATC", type = "Aug Case 3", effect = "Heterogenous Treatment Effect")
  df9 <- data.frame(PE = sim.aug.h.mc[,"ATT.h"],
                    method = "ATT", type = "Aug Case 3", effect = "Heterogenous Treatment Effect")
  
  PEst.aug.h.mc <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  # Augmented estimator case 4
  df1 <- data.frame(PE = sim.aug.c.mm[,"IPW.c"],
                    method = "ATE", type = "Aug Case 4", effect = "Constant Treatment Effect")
  df2 <- data.frame(PE = sim.aug.c.mm[,"IPW.5.c"],
                    method = "ATE 0.05", type = "Aug Case 4", effect = "Constant Treatment Effect")
  df3 <- data.frame(PE = sim.aug.c.mm[,"IPW.10.c"],
                    method = "ATE 0.1", type = "Aug Case 4", effect = "Constant Treatment Effect")
  df4 <- data.frame(PE = sim.aug.c.mm[,"IPW.15.c"],
                    method = "ATE 0.15", type = "Aug Case 4", effect = "Constant Treatment Effect")
  df5 <- data.frame(PE = sim.aug.c.mm[,"ATO.c"],
                    method = "ATO", type = "Aug Case 4", effect = "Constant Treatment Effect")
  df6 <- data.frame(PE = sim.aug.c.mm[,"ATM.c"],
                    method = "ATM", type = "Aug Case 4", effect = "Constant Treatment Effect")
  df7 <- data.frame(PE = sim.aug.c.mm[,"ATEN.c"],
                    method = "ATEN", type = "Aug Case 4", effect = "Constant Treatment Effect")
  df8 <- data.frame(PE = sim.aug.c.mm[,"ATC.c"],
                    method = "ATC", type = "Aug Case 4", effect = "Constant Treatment Effect")
  df9 <- data.frame(PE = sim.aug.c.mm[,"ATT.c"],
                    method = "ATT", type = "Aug Case 4", effect = "Constant Treatment Effect")
  
  PEst.aug.c.mm <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  df1 <- data.frame(PE = sim.aug.h.mm[,"IPW.h"],
                    method = "ATE", type = "Aug Case 4", effect = "Heterogenous Treatment Effect")
  df2 <- data.frame(PE = sim.aug.h.mm[,"IPW.5.h"],
                    method = "ATE 0.05", type = "Aug Case 4", effect = "Heterogenous Treatment Effect")
  df3 <- data.frame(PE = sim.aug.h.mm[,"IPW.10.h"],
                    method = "ATE 0.1", type = "Aug Case 4", effect = "Heterogenous Treatment Effect")
  df4 <- data.frame(PE = sim.aug.h.mm[,"IPW.15.h"],
                    method = "ATE 0.15", type = "Aug Case 4", effect = "Heterogenous Treatment Effect")
  df5 <- data.frame(PE = sim.aug.h.mm[,"ATO.h"],
                    method = "ATO", type = "Aug Case 4", effect = "Heterogenous Treatment Effect")
  df6 <- data.frame(PE = sim.aug.h.mm[,"ATM.h"],
                    method = "ATM", type = "Aug Case 4", effect = "Heterogenous Treatment Effect")
  df7 <- data.frame(PE = sim.aug.h.mm[,"ATEN.h"],
                    method = "ATEN", type = "Aug Case 4", effect = "Heterogenous Treatment Effect")
  df8 <- data.frame(PE = sim.aug.h.mm[,"ATC.h"],
                    method = "ATC", type = "Aug Case 4", effect = "Heterogenous Treatment Effect")
  df9 <- data.frame(PE = sim.aug.h.mm[,"ATT.h"],
                    method = "ATT", type = "Aug Case 4", effect = "Heterogenous Treatment Effect")
  
  PEst.aug.h.mm <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  PEst.df  <- rbind(PEst.wt.c, PEst.wt.h, PEst.aug.c.cc, PEst.aug.h.cc,
                    PEst.aug.c.cm, PEst.aug.h.cm, PEst.aug.c.mc, PEst.aug.h.mc,
                    PEst.aug.c.mm, PEst.aug.h.mm)
  
  PEst.df$method <- factor(PEst.df$method,
                          levels = c("ATE", "ATE 0.05", "ATE 0.1", "ATE 0.15", "ATO", "ATM", "ATEN", "ATC", "ATT"))
  
  PEst.df$type <- factor(PEst.df$type,
                         levels = c("Hajek-like", "Aug Case 1", "Aug Case 2", "Aug Case 3", "Aug Case 4"))
  
  # --- Plot them in the same graph
  ggplot(PEst.df, aes(x = type, y = PE)) + 
    geom_boxplot(aes(fill = method), outlier.size = 0.5) +
    labs(y = "Point Estimation", x = "Estimator and Model Specification", title = plot.title, subtitle = plot.subtitle) + 
    scale_fill_manual(values = c ("royalblue", "deepskyblue", "darkseagreen1", "forestgreen",   
                                  "sienna", "indianred2", "plum2", "goldenrod2", "dimgrey"),
                      name = "Method") + 
    facet_wrap( ~ effect, ncol=1) + ylim(-5, 40) +
    theme_bw()
}

# Function for ARBias plot
abias_plot <- function(md.true.h,
                       sim.wt.c, sim.wt.h, sim.aug.c.cc, sim.aug.h.cc,
                       sim.aug.c.cm, sim.aug.h.cm, sim.aug.c.mc, sim.aug.h.mc,
                       sim.aug.c.mm, sim.aug.h.mm,
                       plot.title = "ARBias Distribution", plot.subtitle = "") {
  
  # Only allowinputting dataset generated in this project (True and Sims files of each model)
  # plot.title: because each plot is for its own case
  
  # Hajek-like (Weighted) estimator
  md.true.c <- 4
  
  df1 <- data.frame(abias = 100*abs(1 - sim.wt.c[,"IPW.c"]/md.true.c),
                    method = "ATE", type = "Hajek-like", effect = "Constant Treatment Effect")
  df2 <- data.frame(abias = 100*abs(1 - sim.wt.c[,"IPW.5.c"]/md.true.c),
                    method = "ATE 0.05", type = "Hajek-like", effect = "Constant Treatment Effect")
  df3 <- data.frame(abias = 100*abs(1 - sim.wt.c[,"IPW.10.c"]/md.true.c),
                    method = "ATE 0.1", type = "Hajek-like", effect = "Constant Treatment Effect")
  df4 <- data.frame(abias = 100*abs(1 - sim.wt.c[,"IPW.15.c"]/md.true.c),
                    method = "ATE 0.15", type = "Hajek-like", effect = "Constant Treatment Effect")
  df5 <- data.frame(abias = 100*abs(1 - sim.wt.c[,"ATO.c"]/md.true.c),
                    method = "ATO", type = "Hajek-like", effect = "Constant Treatment Effect")
  df6 <- data.frame(abias = 100*abs(1 - sim.wt.c[,"ATM.c"]/md.true.c),
                    method = "ATM", type = "Hajek-like", effect = "Constant Treatment Effect")
  df7 <- data.frame(abias = 100*abs(1 - sim.wt.c[,"ATEN.c"]/md.true.c),
                    method = "ATEN", type = "Hajek-like", effect = "Constant Treatment Effect")
  df8 <- data.frame(abias = 100*abs(1 - sim.wt.c[,"ATC.c"]/md.true.c),
                    method = "ATC", type = "Hajek-like", effect = "Constant Treatment Effect")
  df9 <- data.frame(abias = 100*abs(1 - sim.wt.c[,"ATT.c"]/md.true.c),
                    method = "ATT", type = "Hajek-like", effect = "Constant Treatment Effect")
  
  abias.wt.c <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  df1 <- data.frame(abias = 100*abs(1 - sim.wt.h[,"IPW.h"]/md.true.h[,1]),
                    method = "ATE", type = "Hajek-like", effect = "Heterogeneous Treatment Effect")
  df2 <- data.frame(abias = 100*abs(1 - sim.wt.h[,"IPW.5.h"]/md.true.h[,2]),
                    method = "ATE 0.05", type = "Hajek-like", effect = "Heterogeneous Treatment Effect")
  df3 <- data.frame(abias = 100*abs(1 - sim.wt.h[,"IPW.10.h"]/md.true.h[,3]),
                    method = "ATE 0.1", type = "Hajek-like", effect = "Heterogeneous Treatment Effect")
  df4 <- data.frame(abias = 100*abs(1 - sim.wt.h[,"IPW.15.h"]/md.true.h[,4]),
                    method = "ATE 0.15", type = "Hajek-like", effect = "Heterogeneous Treatment Effect")
  df5 <- data.frame(abias = 100*abs(1 - sim.wt.h[,"ATO.h"]/md.true.h[,5]),
                    method = "ATO", type = "Hajek-like", effect = "Heterogeneous Treatment Effect")
  df6 <- data.frame(abias = 100*abs(1 - sim.wt.h[,"ATM.h"]/md.true.h[,6]),
                    method = "ATM", type = "Hajek-like", effect = "Heterogeneous Treatment Effect")
  df7 <- data.frame(abias = 100*abs(1 - sim.wt.h[,"ATEN.h"]/md.true.h[,7]),
                    method = "ATEN", type = "Hajek-like", effect = "Heterogeneous Treatment Effect")
  df8 <- data.frame(abias = 100*abs(1 - sim.wt.h[,"ATC.h"]/md.true.h[,8]),
                    method = "ATC", type = "Hajek-like", effect = "Heterogeneous Treatment Effect")
  df9 <- data.frame(abias = 100*abs(1 - sim.wt.h[,"ATT.h"]/md.true.h[,9]),
                    method = "ATT", type = "Hajek-like", effect = "Heterogeneous Treatment Effect")
  
  abias.wt.h <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  # Augmented estimator case 1
  df1 <- data.frame(abias = 100*abs(1 - sim.aug.c.cc[,"IPW.c"]/md.true.c),
                    method = "ATE", type = "Aug Case 1", effect = "Constant Treatment Effect")
  df2 <- data.frame(abias = 100*abs(1 - sim.aug.c.cc[,"IPW.5.c"]/md.true.c),
                    method = "ATE 0.05", type = "Aug Case 1", effect = "Constant Treatment Effect")
  df3 <- data.frame(abias = 100*abs(1 - sim.aug.c.cc[,"IPW.10.c"]/md.true.c),
                    method = "ATE 0.1", type = "Aug Case 1", effect = "Constant Treatment Effect")
  df4 <- data.frame(abias = 100*abs(1 - sim.aug.c.cc[,"IPW.15.c"]/md.true.c),
                    method = "ATE 0.15", type = "Aug Case 1", effect = "Constant Treatment Effect")
  df5 <- data.frame(abias = 100*abs(1 - sim.aug.c.cc[,"ATO.c"]/md.true.c),
                    method = "ATO", type = "Aug Case 1", effect = "Constant Treatment Effect")
  df6 <- data.frame(abias = 100*abs(1 - sim.aug.c.cc[,"ATM.c"]/md.true.c),
                    method = "ATM", type = "Aug Case 1", effect = "Constant Treatment Effect")
  df7 <- data.frame(abias = 100*abs(1 - sim.aug.c.cc[,"ATEN.c"]/md.true.c),
                    method = "ATEN", type = "Aug Case 1", effect = "Constant Treatment Effect")
  df8 <- data.frame(abias = 100*abs(1 - sim.aug.c.cc[,"ATC.c"]/md.true.c),
                    method = "ATC", type = "Aug Case 1", effect = "Constant Treatment Effect")
  df9 <- data.frame(abias = 100*abs(1 - sim.aug.c.cc[,"ATT.c"]/md.true.c),
                    method = "ATT", type = "Aug Case 1", effect = "Constant Treatment Effect")
  
  abias.aug.c.cc <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  df1 <- data.frame(abias = 100*abs(1 - sim.aug.h.cc[,"IPW.h"]/md.true.h[,1]),
                    method = "ATE", type = "Aug Case 1", effect = "Heterogeneous Treatment Effect")
  df2 <- data.frame(abias = 100*abs(1 - sim.aug.h.cc[,"IPW.5.h"]/md.true.h[,2]),
                    method = "ATE 0.05", type = "Aug Case 1", effect = "Heterogeneous Treatment Effect")
  df3 <- data.frame(abias = 100*abs(1 - sim.aug.h.cc[,"IPW.10.h"]/md.true.h[,3]),
                    method = "ATE 0.1", type = "Aug Case 1", effect = "Heterogeneous Treatment Effect")
  df4 <- data.frame(abias = 100*abs(1 - sim.aug.h.cc[,"IPW.15.h"]/md.true.h[,4]),
                    method = "ATE 0.15", type = "Aug Case 1", effect = "Heterogeneous Treatment Effect")
  df5 <- data.frame(abias = 100*abs(1 - sim.aug.h.cc[,"ATO.h"]/md.true.h[,5]),
                    method = "ATO", type = "Aug Case 1", effect = "Heterogeneous Treatment Effect")
  df6 <- data.frame(abias = 100*abs(1 - sim.aug.h.cc[,"ATM.h"]/md.true.h[,6]),
                    method = "ATM", type = "Aug Case 1", effect = "Heterogeneous Treatment Effect")
  df7 <- data.frame(abias = 100*abs(1 - sim.aug.h.cc[,"ATEN.h"]/md.true.h[,7]),
                    method = "ATEN", type = "Aug Case 1", effect = "Heterogeneous Treatment Effect")
  df8 <- data.frame(abias = 100*abs(1 - sim.aug.h.cc[,"ATC.h"]/md.true.h[,8]),
                    method = "ATC", type = "Aug Case 1", effect = "Heterogeneous Treatment Effect")
  df9 <- data.frame(abias = 100*abs(1 - sim.aug.h.cc[,"ATT.h"]/md.true.h[,9]),
                    method = "ATT", type = "Aug Case 1", effect = "Heterogeneous Treatment Effect")
  
  abias.aug.h.cc <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  # Augmented estimator case 2
  df1 <- data.frame(abias = 100*abs(1 - sim.aug.c.cm[,"IPW.c"]/md.true.c),
                    method = "ATE", type = "Aug Case 2", effect = "Constant Treatment Effect")
  df2 <- data.frame(abias = 100*abs(1 - sim.aug.c.cm[,"IPW.5.c"]/md.true.c),
                    method = "ATE 0.05", type = "Aug Case 2", effect = "Constant Treatment Effect")
  df3 <- data.frame(abias = 100*abs(1 - sim.aug.c.cm[,"IPW.10.c"]/md.true.c),
                    method = "ATE 0.1", type = "Aug Case 2", effect = "Constant Treatment Effect")
  df4 <- data.frame(abias = 100*abs(1 - sim.aug.c.cm[,"IPW.15.c"]/md.true.c),
                    method = "ATE 0.15", type = "Aug Case 2", effect = "Constant Treatment Effect")
  df5 <- data.frame(abias = 100*abs(1 - sim.aug.c.cm[,"ATO.c"]/md.true.c),
                    method = "ATO", type = "Aug Case 2", effect = "Constant Treatment Effect")
  df6 <- data.frame(abias = 100*abs(1 - sim.aug.c.cm[,"ATM.c"]/md.true.c),
                    method = "ATM", type = "Aug Case 2", effect = "Constant Treatment Effect")
  df7 <- data.frame(abias = 100*abs(1 - sim.aug.c.cm[,"ATEN.c"]/md.true.c),
                    method = "ATEN", type = "Aug Case 2", effect = "Constant Treatment Effect")
  df8 <- data.frame(abias = 100*abs(1 - sim.aug.c.cm[,"ATC.c"]/md.true.c),
                    method = "ATC", type = "Aug Case 2", effect = "Constant Treatment Effect")
  df9 <- data.frame(abias = 100*abs(1 - sim.aug.c.cm[,"ATT.c"]/md.true.c),
                    method = "ATT", type = "Aug Case 2", effect = "Constant Treatment Effect")
  
  abias.aug.c.cm <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  df1 <- data.frame(abias = 100*abs(1 - sim.aug.h.cm[,"IPW.h"]/md.true.h[,1]),
                    method = "ATE", type = "Aug Case 2", effect = "Heterogeneous Treatment Effect")
  df2 <- data.frame(abias = 100*abs(1 - sim.aug.h.cm[,"IPW.5.h"]/md.true.h[,2]),
                    method = "ATE 0.05", type = "Aug Case 2", effect = "Heterogeneous Treatment Effect")
  df3 <- data.frame(abias = 100*abs(1 - sim.aug.h.cm[,"IPW.10.h"]/md.true.h[,3]),
                    method = "ATE 0.1", type = "Aug Case 2", effect = "Heterogeneous Treatment Effect")
  df4 <- data.frame(abias = 100*abs(1 - sim.aug.h.cm[,"IPW.15.h"]/md.true.h[,4]),
                    method = "ATE 0.15", type = "Aug Case 2", effect = "Heterogeneous Treatment Effect")
  df5 <- data.frame(abias = 100*abs(1 - sim.aug.h.cm[,"ATO.h"]/md.true.h[,5]),
                    method = "ATO", type = "Aug Case 2", effect = "Heterogeneous Treatment Effect")
  df6 <- data.frame(abias = 100*abs(1 - sim.aug.h.cm[,"ATM.h"]/md.true.h[,6]),
                    method = "ATM", type = "Aug Case 2", effect = "Heterogeneous Treatment Effect")
  df7 <- data.frame(abias = 100*abs(1 - sim.aug.h.cm[,"ATEN.h"]/md.true.h[,7]),
                    method = "ATEN", type = "Aug Case 2", effect = "Heterogeneous Treatment Effect")
  df8 <- data.frame(abias = 100*abs(1 - sim.aug.h.cm[,"ATC.h"]/md.true.h[,8]),
                    method = "ATC", type = "Aug Case 2", effect = "Heterogeneous Treatment Effect")
  df9 <- data.frame(abias = 100*abs(1 - sim.aug.h.cm[,"ATT.h"]/md.true.h[,9]),
                    method = "ATT", type = "Aug Case 2", effect = "Heterogeneous Treatment Effect")
  
  abias.aug.h.cm <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  # Augmented estimator case 3
  df1 <- data.frame(abias = 100*abs(1 - sim.aug.c.mc[,"IPW.c"]/md.true.c),
                    method = "ATE", type = "Aug Case 3", effect = "Constant Treatment Effect")
  df2 <- data.frame(abias = 100*abs(1 - sim.aug.c.mc[,"IPW.5.c"]/md.true.c),
                    method = "ATE 0.05", type = "Aug Case 3", effect = "Constant Treatment Effect")
  df3 <- data.frame(abias = 100*abs(1 - sim.aug.c.mc[,"IPW.10.c"]/md.true.c),
                    method = "ATE 0.1", type = "Aug Case 3", effect = "Constant Treatment Effect")
  df4 <- data.frame(abias = 100*abs(1 - sim.aug.c.mc[,"IPW.15.c"]/md.true.c),
                    method = "ATE 0.15", type = "Aug Case 3", effect = "Constant Treatment Effect")
  df5 <- data.frame(abias = 100*abs(1 - sim.aug.c.mc[,"ATO.c"]/md.true.c),
                    method = "ATO", type = "Aug Case 3", effect = "Constant Treatment Effect")
  df6 <- data.frame(abias = 100*abs(1 - sim.aug.c.mc[,"ATM.c"]/md.true.c), 
                    method = "ATM", type = "Aug Case 3", effect = "Constant Treatment Effect")
  df7 <- data.frame(abias = 100*abs(1 - sim.aug.c.mc[,"ATEN.c"]/md.true.c),
                    method = "ATEN", type = "Aug Case 3", effect = "Constant Treatment Effect")
  df8 <- data.frame(abias = 100*abs(1 - sim.aug.c.mc[,"ATC.c"]/md.true.c),
                    method = "ATC", type = "Aug Case 3", effect = "Constant Treatment Effect")
  df9 <- data.frame(abias = 100*abs(1 - sim.aug.c.mc[,"ATT.c"]/md.true.c),
                    method = "ATT", type = "Aug Case 3", effect = "Constant Treatment Effect")
  
  abias.aug.c.mc <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  df1 <- data.frame(abias = 100*abs(1 - sim.aug.h.mc[,"IPW.h"]/md.true.h[,1]),
                    method = "ATE", type = "Aug Case 3", effect = "Heterogeneous Treatment Effect")
  df2 <- data.frame(abias = 100*abs(1 - sim.aug.h.mc[,"IPW.5.h"]/md.true.h[,2]),
                    method = "ATE 0.05", type = "Aug Case 3", effect = "Heterogeneous Treatment Effect")
  df3 <- data.frame(abias = 100*abs(1 - sim.aug.h.mc[,"IPW.10.h"]/md.true.h[,3]),
                    method = "ATE 0.1", type = "Aug Case 3", effect = "Heterogeneous Treatment Effect")
  df4 <- data.frame(abias = 100*abs(1 - sim.aug.h.mc[,"IPW.15.h"]/md.true.h[,4]),
                    method = "ATE 0.15", type = "Aug Case 3", effect = "Heterogeneous Treatment Effect")
  df5 <- data.frame(abias = 100*abs(1 - sim.aug.h.mc[,"ATO.h"]/md.true.h[,5]),
                    method = "ATO", type = "Aug Case 3", effect = "Heterogeneous Treatment Effect")
  df6 <- data.frame(abias = 100*abs(1 - sim.aug.h.mc[,"ATM.h"]/md.true.h[,6]),
                    method = "ATM", type = "Aug Case 3", effect = "Heterogeneous Treatment Effect")
  df7 <- data.frame(abias = 100*abs(1 - sim.aug.h.mc[,"ATEN.h"]/md.true.h[,7]),
                    method = "ATEN", type = "Aug Case 3", effect = "Heterogeneous Treatment Effect")
  df8 <- data.frame(abias = 100*abs(1 - sim.aug.h.mc[,"ATC.h"]/md.true.h[,8]),
                    method = "ATC", type = "Aug Case 3", effect = "Heterogeneous Treatment Effect")
  df9 <- data.frame(abias = 100*abs(1 - sim.aug.h.mc[,"ATT.h"]/md.true.h[,9]),
                    method = "ATT", type = "Aug Case 3", effect = "Heterogeneous Treatment Effect")
  
  abias.aug.h.mc <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  # Augmented estimator case 4
  df1 <- data.frame(abias = 100*abs(1 - sim.aug.c.mm[,"IPW.c"]/md.true.c),
                    method = "ATE", type = "Aug Case 4", effect = "Constant Treatment Effect")
  df2 <- data.frame(abias = 100*abs(1 - sim.aug.c.mm[,"IPW.5.c"]/md.true.c),
                    method = "ATE 0.05", type = "Aug Case 4", effect = "Constant Treatment Effect")
  df3 <- data.frame(abias = 100*abs(1 - sim.aug.c.mm[,"IPW.10.c"]/md.true.c),
                    method = "ATE 0.1", type = "Aug Case 4", effect = "Constant Treatment Effect")
  df4 <- data.frame(abias = 100*abs(1 - sim.aug.c.mm[,"IPW.15.c"]/md.true.c),
                    method = "ATE 0.15", type = "Aug Case 4", effect = "Constant Treatment Effect")
  df5 <- data.frame(abias = 100*abs(1 - sim.aug.c.mm[,"ATO.c"]/md.true.c),
                    method = "ATO", type = "Aug Case 4", effect = "Constant Treatment Effect")
  df6 <- data.frame(abias = 100*abs(1 - sim.aug.c.mm[,"ATM.c"]/md.true.c),
                    method = "ATM", type = "Aug Case 4", effect = "Constant Treatment Effect")
  df7 <- data.frame(abias = 100*abs(1 - sim.aug.c.mm[,"ATEN.c"]/md.true.c),
                    method = "ATEN", type = "Aug Case 4", effect = "Constant Treatment Effect")
  df8 <- data.frame(abias = 100*abs(1 - sim.aug.c.mm[,"ATC.c"]/md.true.c),
                    method = "ATC", type = "Aug Case 4", effect = "Constant Treatment Effect")
  df9 <- data.frame(abias = 100*abs(1 - sim.aug.c.mm[,"ATT.c"]/md.true.c),
                    method = "ATT", type = "Aug Case 4", effect = "Constant Treatment Effect")
  
  abias.aug.c.mm <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  df1 <- data.frame(abias = 100*abs(1 - sim.aug.h.mm[,"IPW.h"]/md.true.h[,1]),
                    method = "ATE", type = "Aug Case 4", effect = "Heterogeneous Treatment Effect")
  df2 <- data.frame(abias = 100*abs(1 - sim.aug.h.mm[,"IPW.5.h"]/md.true.h[,2]),
                    method = "ATE 0.05", type = "Aug Case 4", effect = "Heterogeneous Treatment Effect")
  df3 <- data.frame(abias = 100*abs(1 - sim.aug.h.mm[,"IPW.10.h"]/md.true.h[,3]),
                    method = "ATE 0.1", type = "Aug Case 4", effect = "Heterogeneous Treatment Effect")
  df4 <- data.frame(abias = 100*abs(1 - sim.aug.h.mm[,"IPW.15.h"]/md.true.h[,4]),
                    method = "ATE 0.15", type = "Aug Case 4", effect = "Heterogeneous Treatment Effect")
  df5 <- data.frame(abias = 100*abs(1 - sim.aug.h.mm[,"ATO.h"]/md.true.h[,5]),
                    method = "ATO", type = "Aug Case 4", effect = "Heterogeneous Treatment Effect")
  df6 <- data.frame(abias = 100*abs(1 - sim.aug.h.mm[,"ATM.h"]/md.true.h[,6]),
                    method = "ATM", type = "Aug Case 4", effect = "Heterogeneous Treatment Effect")
  df7 <- data.frame(abias = 100*abs(1 - sim.aug.h.mm[,"ATEN.h"]/md.true.h[,7]),
                    method = "ATEN", type = "Aug Case 4", effect = "Heterogeneous Treatment Effect")
  df8 <- data.frame(abias = 100*abs(1 - sim.aug.h.mm[,"ATC.h"]/md.true.h[,8]),
                    method = "ATC", type = "Aug Case 4", effect = "Heterogeneous Treatment Effect")
  df9 <- data.frame(abias = 100*abs(1 - sim.aug.h.mm[,"ATT.h"]/md.true.h[,9]),
                    method = "ATT", type = "Aug Case 4", effect = "Heterogeneous Treatment Effect")
  
  abias.aug.h.mm <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  abias.df <- rbind(abias.wt.c, abias.wt.h, abias.aug.c.cc, abias.aug.h.cc,
                    abias.aug.c.cm, abias.aug.h.cm, abias.aug.c.mc, abias.aug.h.mc,
                    abias.aug.c.mm, abias.aug.h.mm)
  
  abias.df$method <- factor(abias.df$method,
                           levels = c("ATE", "ATE 0.05", "ATE 0.1", "ATE 0.15", "ATO", "ATM", "ATEN", "ATC", "ATT"))
  
  abias.df$type <- factor(abias.df$type,
                         levels = c("Hajek-like", "Aug Case 1", "Aug Case 2", "Aug Case 3", "Aug Case 4"))
  
  # --- Plot them in the same graph
  ggplot(abias.df, aes(x = type, y = abias)) + 
    geom_boxplot(aes(fill = method), outlier.size = 0.5) +
    labs(y = "Absolute Relative Percent Bias (%)", x = "Estimator and Model Specification", 
         title = plot.title, subtitle = plot.subtitle) + 
    scale_fill_manual(values = c ("royalblue", "deepskyblue", "darkseagreen1", "forestgreen",   
                                  "sienna", "indianred2", "plum2", "goldenrod2", "dimgrey"),
                      name = "Method") + 
    facet_wrap( ~ effect, ncol=1) + ylim(0, 250) +
    theme_bw()
}

# Function of point estimates
PE <- function(data) {
  apply(data[,seq(1,45, by=5)], 2, mean)
}

# Function of average RBias
avg.rbias <- function(true.h, simdata.c, simdata.h) {
  
  # --- Under constant treatment effect
  col.names.c <- c("IPW.c", "IPW.5.c", "IPW.10.c", "IPW.15.c",
                   "ATO.c", "ATM.c", "ATEN.c", "ATC.c", "ATT.c")
  simdata.c <- simdata.c[, col.names.c]
  true.c <- 4
  meanEst.c <- apply(simdata.c, 2, mean)
  Rbias.avg.c <- 100*abs(1-meanEst.c/true.c)
  
  # --- Under heterogeneous treatment effect
  col.names.h <- c("IPW.h", "IPW.5.h", "IPW.10.h", "IPW.15.h",
                   "ATO.h", "ATM.h", "ATEN.h", "ATC.h", "ATT.h")
  simdata.h <- simdata.h[, col.names.h]
  meanEst.h <- apply(simdata.h, 2, mean)
  Rbias.avg.h <- 100*abs(1-meanEst.h/true.h)
  
  list(Rbias.avg.c = Rbias.avg.c, Rbias.avg.h = Rbias.avg.h)
}

# Function of RMSE
RMSE <- function(true.h, simdata.c, simdata.h) {
  
  # --- Under constant treatment effect
  true.c <- rep(4, 9)
  col.names.c <- c("IPW.c", "IPW.5.c", "IPW.10.c", "IPW.15.c",
                   "ATO.c", "ATM.c", "ATEN.c", "ATC.c", "ATT.c")
  simdata.c <- simdata.c[, col.names.c]
  trans.c <- t(matrix(as.numeric(rep(true.c, nrow(simdata.c))), nrow = ncol(simdata.c)))
  colnames(trans.c) <- colnames(simdata.c)
  sqdiff.c <- (as.matrix(simdata.c) - trans.c)^2
  
  avg.sqdiff.c <- apply(sqdiff.c, 2, mean)
  rmse.c <- sqrt(avg.sqdiff.c)
  
  # --- Under heterogeneous treatment effect
  col.names.h <- c("IPW.h", "IPW.5.h", "IPW.10.h", "IPW.15.h",
                   "ATO.h", "ATM.h", "ATEN.h", "ATC.h", "ATT.h")
  simdata.h <- simdata.h[, col.names.h]
  trans.h <- t(matrix(as.numeric(rep(true.h, nrow(simdata.h))), nrow = ncol(simdata.h)))
  colnames(trans.h) <- colnames(simdata.h)
  sqdiff.h <- (as.matrix(simdata.h) - trans.h)^2
  
  avg.sqdiff.h <- apply(sqdiff.h, 2, mean)
  rmse.h <- sqrt(avg.sqdiff.h)
  
  list(RMSE.c = rmse.c, RMSE.h = rmse.h)
}

# Function of RE (relative efficiency)
RE <- function(simdata.c, simdata.h) {
  
  # --- Under constant treatment effect 
  col.Var.c <- c("IPW.c.var", "IPW.5.c.var", "IPW.10.c.var", "IPW.15.c.var",
                   "ATO.c.var", "ATM.c.var", "ATEN.c.var", "ATC.c.var", "ATT.c.var")
  col.PE.c <- c("IPW.c", "IPW.5.c", "IPW.10.c", "IPW.15.c",
                   "ATO.c", "ATM.c", "ATEN.c", "ATC.c", "ATT.c")
  
  PE.c <- simdata.c[, col.PE.c]
  Var.c <- simdata.c[, col.Var.c]
  
  ep.var.c <- apply(PE.c, 2, var)
  avg.sand.c <- apply(Var.c, 2, mean)
  
  RE.c <- ep.var.c/avg.sand.c
  
  # --- Under heterogeneous treatment effect
  col.Var.h <- c("IPW.h.var", "IPW.5.h.var", "IPW.10.h.var", "IPW.15.h.var",
                   "ATO.h.var", "ATM.h.var", "ATEN.h.var", "ATC.h.var", "ATT.h.var")
  col.PE.h <- c("IPW.h", "IPW.5.h", "IPW.10.h", "IPW.15.h",
                   "ATO.h", "ATM.h", "ATEN.h", "ATC.h", "ATT.h")
  
  PE.h <- simdata.h[, col.PE.h]
  Var.h <- simdata.h[, col.Var.h]
  
  ep.var.h <- apply(PE.h, 2, var)
  avg.sand.h <- apply(Var.h, 2, mean)
  
  RE.h <- ep.var.h/avg.sand.h

  list(RE.c = RE.c, RE.h = RE.h)
}

# Function of CP (coverage probability)
CP <- function(simdata.c, simdata.h) {
  
  # --- Under constant treatment effect 
  col.names.c <- c("IPW.c.ifci", "IPW.5.c.ifci", "IPW.10.c.ifci", "IPW.15.c.ifci",
                   "ATO.c.ifci", "ATM.c.ifci", "ATEN.c.ifci", "ATC.c.ifci", "ATT.c.ifci")
  simdata.c <- simdata.c[, col.names.c]
  CP.c <- apply(simdata.c, 2, mean)
  
  # --- Under heterogeneous treatment effect
  col.names.h <- c("IPW.h.ifci", "IPW.5.h.ifci", "IPW.10.h.ifci", "IPW.15.h.ifci",
                   "ATO.h.ifci", "ATM.h.ifci", "ATEN.h.ifci", "ATC.h.ifci", "ATT.h.ifci")
  simdata.h <- simdata.h[, col.names.h]
  CP.h <- apply(simdata.h, 2, mean)

  list(CP.c = CP.c, CP.h = CP.h)
}
