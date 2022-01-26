### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Variance estimations ATE ATT ATC ~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

####### Summary statistics (Tables for paper)

### by Yi Liu
### Create date: Nov 20, 2021

library(xtable)
source("SumStat_func.R")

load("truth.RData")
names(md1.true_eff) <- c("ATE", "ATT", "ATC")
names(md2.true_eff) <- c("ATE", "ATT", "ATC")
names(md3.true_eff) <- c("ATE", "ATT", "ATC")
truth.c <- data.frame(ATC=4, ATT=4)

load("md1_sims_WBoot.RData")
load("md2_sims_WBoot.RData")
load("md3_sims_WBoot.RData")
load("md1_sims_WBoot_cons.RData")
load("md2_sims_WBoot_cons.RData")
load("md3_sims_WBoot_cons.RData")

newnames <- c("ATE", "ATE.wbvar", "ATE.lwr", "ATE.upr", "ATE.ifci",
              "ATC", "ATC.wbvar", "ATC.lwr", "ATC.upr", "ATC.ifci",
              "ATT", "ATT.wbvar", "ATT.lwr", "ATT.upr", "ATT.ifci",
              "Row.Num")


### Model 1
load("md1_sims.RData")
colnames(c.sim.aug.new.cc) <- newnames
colnames(c.sim.aug.new.cm) <- newnames
colnames(c.sim.aug.new.mc) <- newnames
colnames(c.sim.aug.new.mm) <- newnames
colnames(sim.aug.new.cc) <- newnames
colnames(sim.aug.new.cm) <- newnames
colnames(sim.aug.new.mc) <- newnames
colnames(sim.aug.new.mm) <- newnames


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(4, 6),
                    #Est = rep(c(mean(wb1.rad.c.cc$ATC), mean(wb1.rad.c.cc$ATT)), 3),
                    ARBias = rep(bias(truth.c, wb1.rad.c.cc), 3),
                    MRBias = rep(bias.med(truth.c, wb1.rad.c.cc), 3),
                    RMSE = rep(RMSE(truth.c, wb1.rad.c.cc), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(c.sim.aug.new.cc$ATC.wbvar), mean(c.sim.aug.new.cc$ATT.wbvar),
                                mean(wb1.rad.c.cc$ATC.wbvar), mean(wb1.rad.c.cc$ATT.wbvar),
                                mean(wb1.exp.c.cc$ATC.wbvar), mean(wb1.exp.c.cc$ATT.wbvar))),
                    SD.med = sqrt(c(median(c.sim.aug.new.cc$ATC.wbvar), median(c.sim.aug.new.cc$ATT.wbvar),
                                    median(wb1.rad.c.cc$ATC.wbvar), median(wb1.rad.c.cc$ATT.wbvar),
                                    median(wb1.exp.c.cc$ATC.wbvar), median(wb1.exp.c.cc$ATT.wbvar))),
                    ESD = rep(c(sd(wb1.rad.c.cc$ATC), sd(wb1.rad.c.cc$ATT)), 3),
                    RE = c(RE(c.sim.aug.new.cc), RE(wb1.rad.c.cc), RE(wb1.exp.c.cc)),
                    RE.med = c(RE.med(c.sim.aug.new.cc), RE.med(wb1.rad.c.cc), RE.med(wb1.exp.c.cc)),
                    CP = c(CP(truth.c, c.sim.aug.new.cc),
                           CP(truth.c, wb1.rad.c.cc), 
                           CP(truth.c, wb1.exp.c.cc)))
print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(as.numeric(md1.true_eff[, c("ATC", "ATT")]), 3),
                    #Est = rep(c(mean(wb1.rad.h.cc$ATC), mean(wb1.rad.h.cc$ATT)), 3),
                    ARBias = rep(bias(md1.true_eff, wb1.rad.h.cc), 3),
                    MRBias = rep(bias.med(md1.true_eff, wb1.rad.h.cc), 3),
                    RMSE = rep(RMSE(md1.true_eff, wb1.rad.h.cc), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(sim.aug.new.cc$ATC.wbvar), mean(sim.aug.new.cc$ATT.wbvar),
                                mean(wb1.rad.h.cc$ATC.wbvar), mean(wb1.rad.h.cc$ATT.wbvar),
                                mean(wb1.exp.h.cc$ATC.wbvar), mean(wb1.exp.h.cc$ATT.wbvar))),
                    SD.med = sqrt(c(median(sim.aug.new.cc$ATC.wbvar), median(sim.aug.new.cc$ATT.wbvar),
                                    median(wb1.rad.h.cc$ATC.wbvar), median(wb1.rad.h.cc$ATT.wbvar),
                                    median(wb1.exp.h.cc$ATC.wbvar), median(wb1.exp.h.cc$ATT.wbvar))),
                    ESD = rep(c(sd(wb1.rad.h.cc$ATC), sd(wb1.rad.h.cc$ATT)), 3),
                    RE = c(RE(sim.aug.new.cc), RE(wb1.rad.h.cc), RE(wb1.exp.h.cc)),
                    RE.med = c(RE.med(sim.aug.new.cc), RE.med(wb1.rad.h.cc), RE.med(wb1.exp.h.cc)),
                    CP = c(CP(md1.true_eff, sim.aug.new.cc),
                           CP(md1.true_eff, wb1.rad.h.cc), 
                           CP(md1.true_eff, wb1.exp.h.cc)))

print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(4, 6),
                    #Est = rep(c(mean(wb1.rad.c.cm$ATC), mean(wb1.rad.c.cm$ATT)), 3),
                    ARBias = rep(bias(truth.c, wb1.rad.c.cm), 3),
                    MRBias = rep(bias.med(truth.c, wb1.rad.c.cm), 3),
                    RMSE = rep(RMSE(truth.c, wb1.rad.c.cm), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(c.sim.aug.new.cm$ATC.wbvar), mean(c.sim.aug.new.cm$ATT.wbvar),
                                mean(wb1.rad.c.cm$ATC.wbvar), mean(wb1.rad.c.cm$ATT.wbvar),
                                mean(wb1.exp.c.cm$ATC.wbvar), mean(wb1.exp.c.cm$ATT.wbvar))),
                    SD.med = sqrt(c(median(c.sim.aug.new.cm$ATC.wbvar), median(c.sim.aug.new.cm$ATT.wbvar),
                                    median(wb1.rad.c.cm$ATC.wbvar), median(wb1.rad.c.cm$ATT.wbvar),
                                    median(wb1.exp.c.cm$ATC.wbvar), median(wb1.exp.c.cm$ATT.wbvar))),
                    ESD = rep(c(sd(wb1.rad.c.cm$ATC), sd(wb1.rad.c.cm$ATT)), 3),
                    RE = c(RE(c.sim.aug.new.cm), RE(wb1.rad.c.cm), RE(wb1.exp.c.cm)),
                    RE.med = c(RE.med(c.sim.aug.new.cm), RE.med(wb1.rad.c.cm), RE.med(wb1.exp.c.cm)),
                    CP = c(CP(truth.c, c.sim.aug.new.cm),
                           CP(truth.c, wb1.rad.c.cm), 
                           CP(truth.c, wb1.exp.c.cm)))
print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(as.numeric(md1.true_eff[, c("ATC", "ATT")]), 3),
                    #Est = rep(c(mean(wb1.rad.h.cm$ATC), mean(wb1.rad.h.cm$ATT)), 3),
                    ARBias = rep(bias(md1.true_eff, wb1.rad.h.cm), 3),
                    MRBias = rep(bias.med(md1.true_eff, wb1.rad.h.cm), 3),
                    RMSE = rep(RMSE(md1.true_eff, wb1.rad.h.cm), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(sim.aug.new.cm$ATC.wbvar), mean(sim.aug.new.cm$ATT.wbvar),
                                mean(wb1.rad.h.cm$ATC.wbvar), mean(wb1.rad.h.cm$ATT.wbvar),
                                mean(wb1.exp.h.cm$ATC.wbvar), mean(wb1.exp.h.cm$ATT.wbvar))),
                    SD.med = sqrt(c(median(sim.aug.new.cm$ATC.wbvar), median(sim.aug.new.cm$ATT.wbvar),
                                    median(wb1.rad.h.cm$ATC.wbvar), median(wb1.rad.h.cm$ATT.wbvar),
                                    median(wb1.exp.h.cm$ATC.wbvar), median(wb1.exp.h.cm$ATT.wbvar))),
                    ESD = rep(c(sd(wb1.rad.h.cm$ATC), sd(wb1.rad.h.cm$ATT)), 3),
                    RE = c(RE(sim.aug.new.cm), RE(wb1.rad.h.cm), RE(wb1.exp.h.cm)),
                    RE.med = c(RE.med(sim.aug.new.cm), RE.med(wb1.rad.h.cm), RE.med(wb1.exp.h.cm)),
                    CP = c(CP(md1.true_eff, sim.aug.new.cm),
                           CP(md1.true_eff, wb1.rad.h.cm), 
                           CP(md1.true_eff, wb1.exp.h.cm)))
print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(4, 6),
                    #Est = rep(c(mean(wb1.rad.c.mc$ATC), mean(wb1.rad.c.mc$ATT)), 3),
                    ARBias = rep(bias(truth.c, wb1.rad.c.mc), 3),
                    MRBias = rep(bias(truth.c, wb1.rad.c.mc), 3),
                    RMSE = rep(RMSE(truth.c, wb1.rad.c.mc), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(c.sim.aug.new.mc$ATC.wbvar), mean(c.sim.aug.new.mc$ATT.wbvar),
                                mean(wb1.rad.c.mc$ATC.wbvar), mean(wb1.rad.c.mc$ATT.wbvar),
                                mean(wb1.exp.c.mc$ATC.wbvar), mean(wb1.exp.c.mc$ATT.wbvar))),
                    SD.med = sqrt(c(median(c.sim.aug.new.mc$ATC.wbvar), median(c.sim.aug.new.mc$ATT.wbvar),
                                    median(wb1.rad.c.mc$ATC.wbvar), median(wb1.rad.c.mc$ATT.wbvar),
                                    median(wb1.exp.c.mc$ATC.wbvar), median(wb1.exp.c.mc$ATT.wbvar))),
                    ESD = rep(c(sd(wb1.rad.c.mc$ATC), sd(wb1.rad.c.mc$ATT)), 3),
                    RE = c(RE(c.sim.aug.new.mc), RE(wb1.rad.c.mc), RE(wb1.exp.c.mc)),
                    RE.med = c(RE.med(c.sim.aug.new.mc), RE.med(wb1.rad.c.mc), RE.med(wb1.exp.c.mc)),
                    CP = c(CP(truth.c, c.sim.aug.new.mc),
                           CP(truth.c, wb1.rad.c.mc), 
                           CP(truth.c, wb1.exp.c.mc)))
print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(as.numeric(md1.true_eff[, c("ATC", "ATT")]), 3),
                    #Est = rep(c(mean(wb1.rad.h.mc$ATC), mean(wb1.rad.h.mc$ATT)), 3),
                    ARBias = rep(bias(md1.true_eff, wb1.rad.h.mc), 3),
                    MRBias = rep(bias.med(md1.true_eff, wb1.rad.h.mc), 3),
                    RMSE = rep(RMSE(md1.true_eff, wb1.rad.h.mc), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(sim.aug.new.mc$ATC.wbvar), mean(sim.aug.new.mc$ATT.wbvar),
                                mean(wb1.rad.h.mc$ATC.wbvar), mean(wb1.rad.h.mc$ATT.wbvar),
                                mean(wb1.exp.h.mc$ATC.wbvar), mean(wb1.exp.h.mc$ATT.wbvar))),
                    SD.med = sqrt(c(median(sim.aug.new.mc$ATC.wbvar), median(sim.aug.new.mc$ATT.wbvar),
                                    median(wb1.rad.h.mc$ATC.wbvar), median(wb1.rad.h.mc$ATT.wbvar),
                                    median(wb1.exp.h.mc$ATC.wbvar), median(wb1.exp.h.mc$ATT.wbvar))),
                    ESD = rep(c(sd(wb1.rad.h.mc$ATC), sd(wb1.rad.h.mc$ATT)), 3),
                    RE = c(RE(sim.aug.new.mc), RE(wb1.rad.h.mc), RE(wb1.exp.h.mc)),
                    RE.med = c(RE.med(sim.aug.new.mc), RE.med(wb1.rad.h.mc), RE.med(wb1.exp.h.mc)),
                    CP = c(CP(md1.true_eff, sim.aug.new.mc),
                           CP(md1.true_eff, wb1.rad.h.mc), 
                           CP(md1.true_eff, wb1.exp.h.mc)))
print(xtable(table), include.rownames=FALSE)

table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(4, 6),
                    #Est = rep(c(mean(wb1.rad.c.mm$ATC), mean(wb1.rad.c.mm$ATT)), 3),
                    ARBias = rep(bias(truth.c, wb1.rad.c.mm), 3),
                    MRBias = rep(bias.med(truth.c, wb1.rad.c.mm), 3),
                    RMSE = rep(RMSE(truth.c, wb1.rad.c.mm), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(c.sim.aug.new.mm$ATC.wbvar), mean(c.sim.aug.new.mm$ATT.wbvar),
                                mean(wb1.rad.c.mm$ATC.wbvar), mean(wb1.rad.c.mm$ATT.wbvar),
                                mean(wb1.exp.c.mm$ATC.wbvar), mean(wb1.exp.c.mm$ATT.wbvar))),
                    SD.med = sqrt(c(median(c.sim.aug.new.mm$ATC.wbvar), median(c.sim.aug.new.mm$ATT.wbvar),
                                    median(wb1.rad.c.mm$ATC.wbvar), median(wb1.rad.c.mm$ATT.wbvar),
                                    median(wb1.exp.c.mm$ATC.wbvar), median(wb1.exp.c.mm$ATT.wbvar))),
                    ESD = rep(c(sd(wb1.rad.c.mm$ATC), sd(wb1.rad.c.mm$ATT)), 3),
                    RE = c(RE(c.sim.aug.new.mm), RE(wb1.rad.c.mm), RE(wb1.exp.c.mm)),
                    RE.med = c(RE.med(c.sim.aug.new.mm), RE.med(wb1.rad.c.mm), RE.med(wb1.exp.c.mm)),
                    CP = c(CP(truth.c, c.sim.aug.new.mm),
                           CP(truth.c, wb1.rad.c.mm), 
                           CP(truth.c, wb1.exp.c.mm)))
print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(as.numeric(md1.true_eff[, c("ATC", "ATT")]), 3),
                    #Est = rep(c(mean(wb1.rad.h.mm$ATC), mean(wb1.rad.h.mm$ATT)), 3),
                    ARBias = rep(bias(md1.true_eff, wb1.rad.h.mm), 3),
                    MRBias = rep(bias.med(md1.true_eff, wb1.rad.h.mm), 3),
                    RMSE = rep(RMSE(md1.true_eff, wb1.rad.h.mm), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(sim.aug.new.mm$ATC.wbvar), mean(sim.aug.new.mm$ATT.wbvar),
                                mean(wb1.rad.h.mm$ATC.wbvar), mean(wb1.rad.h.mm$ATT.wbvar),
                                mean(wb1.exp.h.mm$ATC.wbvar), mean(wb1.exp.h.mm$ATT.wbvar))),
                    SD.med = sqrt(c(median(sim.aug.new.mm$ATC.wbvar), median(sim.aug.new.mm$ATT.wbvar),
                                    median(wb1.rad.h.mm$ATC.wbvar), median(wb1.rad.h.mm$ATT.wbvar),
                                    median(wb1.exp.h.mm$ATC.wbvar), median(wb1.exp.h.mm$ATT.wbvar))),
                    ESD = rep(c(sd(wb1.rad.h.mm$ATC), sd(wb1.rad.h.mm$ATT)), 3),
                    RE = c(RE(sim.aug.new.mm), RE(wb1.rad.h.mm), RE(wb1.exp.h.mm)),
                    RE.med = c(RE.med(sim.aug.new.mm), RE.med(wb1.rad.h.mm), RE.med(wb1.exp.h.mm)),
                    CP = c(CP(md1.true_eff, sim.aug.new.mm),
                           CP(md1.true_eff, wb1.rad.h.mm), 
                           CP(md1.true_eff, wb1.exp.h.mm)))
print(xtable(table), include.rownames=FALSE)


### Model 2
load("md2_sims.RData")
colnames(c.sim.aug.new.cc) <- newnames
colnames(c.sim.aug.new.cm) <- newnames
colnames(c.sim.aug.new.mc) <- newnames
colnames(c.sim.aug.new.mm) <- newnames
colnames(sim.aug.new.cc) <- newnames
colnames(sim.aug.new.cm) <- newnames
colnames(sim.aug.new.mc) <- newnames
colnames(sim.aug.new.mm) <- newnames


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(4, 6),
                    #Est = rep(c(mean(wb2.rad.c.cc$ATC), mean(wb2.rad.c.cc$ATT)), 3),
                    ARBias = rep(bias(truth.c, wb2.rad.c.cc), 3),
                    MRBias = rep(bias.med(truth.c, wb2.rad.c.cc), 3),
                    RMSE = rep(RMSE(truth.c, wb2.rad.c.cc), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(c.sim.aug.new.cc$ATC.wbvar), mean(c.sim.aug.new.cc$ATT.wbvar),
                                mean(wb2.rad.c.cc$ATC.wbvar), mean(wb2.rad.c.cc$ATT.wbvar),
                                mean(wb2.exp.c.cc$ATC.wbvar), mean(wb2.exp.c.cc$ATT.wbvar))),
                    SD.med = sqrt(c(median(c.sim.aug.new.cc$ATC.wbvar), median(c.sim.aug.new.cc$ATT.wbvar),
                                    median(wb2.rad.c.cc$ATC.wbvar), median(wb2.rad.c.cc$ATT.wbvar),
                                    median(wb2.exp.c.cc$ATC.wbvar), median(wb2.exp.c.cc$ATT.wbvar))),
                    ESD = rep(c(sd(wb2.rad.c.cc$ATC), sd(wb2.rad.c.cc$ATT)), 3),
                    RE = c(RE(c.sim.aug.new.cc), RE(wb2.rad.c.cc), RE(wb2.exp.c.cc)),
                    RE.med = c(RE.med(c.sim.aug.new.cc), RE.med(wb2.rad.c.cc), RE.med(wb2.exp.c.cc)),
                    CP = c(CP(truth.c, c.sim.aug.new.cc),
                           CP(truth.c, wb2.rad.c.cc), 
                           CP(truth.c, wb2.exp.c.cc)))
print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(as.numeric(md2.true_eff[, c("ATC", "ATT")]), 3),
                    #Est = rep(c(mean(wb2.rad.h.cc$ATC), mean(wb2.rad.h.cc$ATT)), 3),
                    ARBias = rep(bias(md2.true_eff, wb2.rad.h.cc), 3),
                    MRBias = rep(bias.med(md2.true_eff, wb2.rad.h.cc), 3),
                    RMSE = rep(RMSE(md2.true_eff, wb2.rad.h.cc), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(sim.aug.new.cc$ATC.wbvar), mean(sim.aug.new.cc$ATT.wbvar),
                                mean(wb2.rad.h.cc$ATC.wbvar), mean(wb2.rad.h.cc$ATT.wbvar),
                                mean(wb2.exp.h.cc$ATC.wbvar), mean(wb2.exp.h.cc$ATT.wbvar))),
                    SD.med = sqrt(c(median(sim.aug.new.cc$ATC.wbvar), median(sim.aug.new.cc$ATT.wbvar),
                                    median(wb2.rad.h.cc$ATC.wbvar), median(wb2.rad.h.cc$ATT.wbvar),
                                    median(wb2.exp.h.cc$ATC.wbvar), median(wb2.exp.h.cc$ATT.wbvar))),
                    ESD = rep(c(sd(wb2.rad.h.cc$ATC), sd(wb2.rad.h.cc$ATT)), 3),
                    RE = c(RE(sim.aug.new.cc), RE(wb2.rad.h.cc), RE(wb2.exp.h.cc)),
                    RE.med = c(RE.med(sim.aug.new.cc), RE.med(wb2.rad.h.cc), RE.med(wb2.exp.h.cc)),
                    CP = c(CP(md2.true_eff, sim.aug.new.cc),
                           CP(md2.true_eff, wb2.rad.h.cc), 
                           CP(md2.true_eff, wb2.exp.h.cc)))

print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(4, 6),
                    #Est = rep(c(mean(wb2.rad.c.cm$ATC), mean(wb2.rad.c.cm$ATT)), 3),
                    ARBias = rep(bias(truth.c, wb2.rad.c.cm), 3),
                    MRBias = rep(bias.med(truth.c, wb2.rad.c.cm), 3),
                    RMSE = rep(RMSE(truth.c, wb2.rad.c.cm), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(c.sim.aug.new.cm$ATC.wbvar), mean(c.sim.aug.new.cm$ATT.wbvar),
                                mean(wb2.rad.c.cm$ATC.wbvar), mean(wb2.rad.c.cm$ATT.wbvar),
                                mean(wb2.exp.c.cm$ATC.wbvar), mean(wb2.exp.c.cm$ATT.wbvar))),
                    SD.med = sqrt(c(median(c.sim.aug.new.cm$ATC.wbvar), median(c.sim.aug.new.cm$ATT.wbvar),
                                    median(wb2.rad.c.cm$ATC.wbvar), median(wb2.rad.c.cm$ATT.wbvar),
                                    median(wb2.exp.c.cm$ATC.wbvar), median(wb2.exp.c.cm$ATT.wbvar))),
                    ESD = rep(c(sd(wb2.rad.c.cm$ATC), sd(wb2.rad.c.cm$ATT)), 3),
                    RE = c(RE(c.sim.aug.new.cm), RE(wb2.rad.c.cm), RE(wb2.exp.c.cm)),
                    RE.med = c(RE.med(c.sim.aug.new.cm), RE.med(wb2.rad.c.cm), RE.med(wb2.exp.c.cm)),
                    CP = c(CP(truth.c, c.sim.aug.new.cm),
                           CP(truth.c, wb2.rad.c.cm), 
                           CP(truth.c, wb2.exp.c.cm)))
print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(as.numeric(md2.true_eff[, c("ATC", "ATT")]), 3),
                    #Est = rep(c(mean(wb2.rad.h.cm$ATC), mean(wb2.rad.h.cm$ATT)), 3),
                    ARBias = rep(bias(md2.true_eff, wb2.rad.h.cm), 3),
                    MRBias = rep(bias.med(md2.true_eff, wb2.rad.h.cm), 3),
                    RMSE = rep(RMSE(md2.true_eff, wb2.rad.h.cm), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(sim.aug.new.cm$ATC.wbvar), mean(sim.aug.new.cm$ATT.wbvar),
                                mean(wb2.rad.h.cm$ATC.wbvar), mean(wb2.rad.h.cm$ATT.wbvar),
                                mean(wb2.exp.h.cm$ATC.wbvar), mean(wb2.exp.h.cm$ATT.wbvar))),
                    SD.med = sqrt(c(median(sim.aug.new.cm$ATC.wbvar), median(sim.aug.new.cm$ATT.wbvar),
                                    median(wb2.rad.h.cm$ATC.wbvar), median(wb2.rad.h.cm$ATT.wbvar),
                                    median(wb2.exp.h.cm$ATC.wbvar), median(wb2.exp.h.cm$ATT.wbvar))),
                    ESD = rep(c(sd(wb2.rad.h.cm$ATC), sd(wb2.rad.h.cm$ATT)), 3),
                    RE = c(RE(c.sim.aug.new.cm), RE(wb2.rad.h.cm), RE(wb2.exp.h.cm)),
                    RE.med = c(RE.med(c.sim.aug.new.cm), RE.med(wb2.rad.h.cm), RE.med(wb2.exp.h.cm)),
                    CP = c(CP(md2.true_eff, sim.aug.new.cm),
                           CP(md2.true_eff, wb2.rad.h.cm), 
                           CP(md2.true_eff, wb2.exp.h.cm)))
print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(4, 6),
                    #Est = rep(c(mean(wb2.rad.c.mc$ATC), mean(wb2.rad.c.mc$ATT)), 3),
                    ARBias = rep(bias(truth.c, wb2.rad.c.mc), 3),
                    MRBias = rep(bias.med(truth.c, wb2.rad.c.mc), 3),
                    RMSE = rep(RMSE(truth.c, wb2.rad.c.mc), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(c.sim.aug.new.mc$ATC.wbvar), mean(c.sim.aug.new.mc$ATT.wbvar),
                                mean(wb2.rad.c.mc$ATC.wbvar), mean(wb2.rad.c.mc$ATT.wbvar),
                                mean(wb2.exp.c.mc$ATC.wbvar), mean(wb2.exp.c.mc$ATT.wbvar))),
                    SD.med = sqrt(c(median(c.sim.aug.new.mc$ATC.wbvar), median(c.sim.aug.new.mc$ATT.wbvar),
                                    median(wb2.rad.c.mc$ATC.wbvar), median(wb2.rad.c.mc$ATT.wbvar),
                                    median(wb2.exp.c.mc$ATC.wbvar), median(wb2.exp.c.mc$ATT.wbvar))),
                    ESD = rep(c(sd(wb2.rad.c.mc$ATC), sd(wb2.rad.c.mc$ATT)), 3),
                    RE = c(RE(c.sim.aug.new.mc), RE(wb2.rad.c.mc), RE(wb2.exp.c.mc)),
                    RE.med = c(RE.med(c.sim.aug.new.mc), RE.med(wb2.rad.c.mc), RE.med(wb2.exp.c.mc)),
                    CP = c(CP(truth.c, c.sim.aug.new.mc),
                           CP(truth.c, wb2.rad.c.mc), 
                           CP(truth.c, wb2.exp.c.mc)))
print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(as.numeric(md2.true_eff[, c("ATC", "ATT")]), 3),
                    #Est = rep(c(mean(wb2.rad.h.mc$ATC), mean(wb2.rad.h.mc$ATT)), 3),
                    ARBias = rep(bias(md2.true_eff, wb2.rad.h.mc), 3),
                    MRBias = rep(bias.med(md2.true_eff, wb2.rad.h.mc), 3),
                    RMSE = rep(RMSE(md2.true_eff, wb2.rad.h.mc), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(sim.aug.new.mc$ATC.wbvar), mean(sim.aug.new.mc$ATT.wbvar),
                                mean(wb2.rad.h.mc$ATC.wbvar), mean(wb2.rad.h.mc$ATT.wbvar),
                                mean(wb2.exp.h.mc$ATC.wbvar), mean(wb2.exp.h.mc$ATT.wbvar))),
                    SD.med = sqrt(c(median(sim.aug.new.mc$ATC.wbvar), median(sim.aug.new.mc$ATT.wbvar),
                                    median(wb2.rad.h.mc$ATC.wbvar), median(wb2.rad.h.mc$ATT.wbvar),
                                    median(wb2.exp.h.mc$ATC.wbvar), median(wb2.exp.h.mc$ATT.wbvar))),
                    ESD = rep(c(sd(wb2.rad.h.mc$ATC), sd(wb2.rad.h.mc$ATT)), 3),
                    RE = c(RE(sim.aug.new.mc), RE(wb2.rad.h.mc), RE(wb2.exp.h.mc)),
                    RE.med = c(RE.med(sim.aug.new.mc), RE.med(wb2.rad.h.mc), RE.med(wb2.exp.h.mc)),
                    CP = c(CP(md2.true_eff, sim.aug.new.mc),
                           CP(md2.true_eff, wb2.rad.h.mc), 
                           CP(md2.true_eff, wb2.exp.h.mc)))
print(xtable(table), include.rownames=FALSE)

table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(4, 6),
                    #Est = rep(c(mean(wb2.rad.c.mm$ATC), mean(wb2.rad.c.mm$ATT)), 3),
                    ARBias = rep(bias(truth.c, wb2.rad.c.mm), 3),
                    MRBias = rep(bias.med(truth.c, wb2.rad.c.mm), 3),
                    RMSE = rep(RMSE(truth.c, wb2.rad.c.mm), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(c.sim.aug.new.mm$ATC.wbvar), mean(c.sim.aug.new.mm$ATT.wbvar),
                                mean(wb2.rad.c.mm$ATC.wbvar), mean(wb2.rad.c.mm$ATT.wbvar),
                                mean(wb2.exp.c.mm$ATC.wbvar), mean(wb2.exp.c.mm$ATT.wbvar))),
                    SD.med = sqrt(c(median(c.sim.aug.new.mm$ATC.wbvar), median(c.sim.aug.new.mm$ATT.wbvar),
                                    median(wb2.rad.c.mm$ATC.wbvar), median(wb2.rad.c.mm$ATT.wbvar),
                                    median(wb2.exp.c.mm$ATC.wbvar), median(wb2.exp.c.mm$ATT.wbvar))),
                    ESD = rep(c(sd(wb2.rad.c.mm$ATC), sd(wb2.rad.c.mm$ATT)), 3),
                    RE = c(RE(c.sim.aug.new.mm), RE(wb2.rad.c.mm), RE(wb2.exp.c.mm)),
                    RE.med = c(RE.med(c.sim.aug.new.mm), RE.med(wb2.rad.c.mm), RE.med(wb2.exp.c.mm)),
                    CP = c(CP(truth.c, c.sim.aug.new.mm),
                           CP(truth.c, wb2.rad.c.mm), 
                           CP(truth.c, wb2.exp.c.mm)))
print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(as.numeric(md2.true_eff[, c("ATC", "ATT")]), 3),
                    #Est = rep(c(mean(wb2.rad.h.mm$ATC), mean(wb2.rad.h.mm$ATT)), 3),
                    ARBias = rep(bias(md2.true_eff, wb2.rad.h.mm), 3),
                    MRBias = rep(bias.med(md2.true_eff, wb2.rad.h.mm), 3),
                    RMSE = rep(RMSE(md2.true_eff, wb2.rad.h.mm), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(sim.aug.new.mm$ATC.wbvar), mean(sim.aug.new.mm$ATT.wbvar),
                                mean(wb2.rad.h.mm$ATC.wbvar), mean(wb2.rad.h.mm$ATT.wbvar),
                                mean(wb2.exp.h.mm$ATC.wbvar), mean(wb2.exp.h.mm$ATT.wbvar))),
                    SD.med = sqrt(c(median(sim.aug.new.mm$ATC.wbvar), median(sim.aug.new.mm$ATT.wbvar),
                                    median(wb2.rad.h.mm$ATC.wbvar), median(wb2.rad.h.mm$ATT.wbvar),
                                    median(wb2.exp.h.mm$ATC.wbvar), median(wb2.exp.h.mm$ATT.wbvar))),
                    ESD = rep(c(sd(wb2.rad.h.mm$ATC), sd(wb2.rad.h.mm$ATT)), 3),
                    RE = c(RE(sim.aug.new.mm), RE(wb2.rad.h.mm), RE(wb2.exp.h.mm)),
                    RE.med = c(RE.med(sim.aug.new.mm), RE.med(wb2.rad.h.mm), RE.med(wb2.exp.h.mm)),
                    CP = c(CP(md2.true_eff, sim.aug.new.mm),
                           CP(md2.true_eff, wb2.rad.h.mm), 
                           CP(md2.true_eff, wb2.exp.h.mm)))
print(xtable(table), include.rownames=FALSE)

### Model 3
load("md3_sims.RData")
colnames(c.sim.aug.new.cc) <- newnames
colnames(c.sim.aug.new.cm) <- newnames
colnames(c.sim.aug.new.mc) <- newnames
colnames(c.sim.aug.new.mm) <- newnames
colnames(sim.aug.new.cc) <- newnames
colnames(sim.aug.new.cm) <- newnames
colnames(sim.aug.new.mc) <- newnames
colnames(sim.aug.new.mm) <- newnames


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(4, 6),
                    #Est = rep(c(mean(wb3.rad.c.cc$ATC), mean(wb3.rad.c.cc$ATT)), 3),
                    ARBias = rep(bias(truth.c, wb3.rad.c.cc), 3),
                    MRBias = rep(bias(truth.c, wb3.rad.c.cc), 3),
                    RMSE = rep(RMSE(truth.c, wb3.rad.c.cc), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(c.sim.aug.new.cc$ATC.wbvar), mean(c.sim.aug.new.cc$ATT.wbvar),
                                mean(wb3.rad.c.cc$ATC.wbvar), mean(wb3.rad.c.cc$ATT.wbvar),
                                mean(wb3.exp.c.cc$ATC.wbvar), mean(wb3.exp.c.cc$ATT.wbvar))),
                    SD.med = sqrt(c(median(c.sim.aug.new.cc$ATC.wbvar), median(c.sim.aug.new.cc$ATT.wbvar),
                                    median(wb3.rad.c.cc$ATC.wbvar), median(wb3.rad.c.cc$ATT.wbvar),
                                    median(wb3.exp.c.cc$ATC.wbvar), median(wb3.exp.c.cc$ATT.wbvar))),
                    ESD = rep(c(sd(wb3.rad.c.cc$ATC), sd(wb3.rad.c.cc$ATT)), 3),
                    RE = c(RE(c.sim.aug.new.cc), RE(wb3.rad.c.cc), RE(wb3.exp.c.cc)),
                    RE.med = c(RE.med(c.sim.aug.new.cc), RE.med(wb3.rad.c.cc), RE.med(wb3.exp.c.cc)),
                    CP = c(CP(truth.c, c.sim.aug.new.cc),
                           CP(truth.c, wb3.rad.c.cc), 
                           CP(truth.c, wb3.exp.c.cc)))
print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(as.numeric(md3.true_eff[, c("ATC", "ATT")]), 3),
                    #Est = rep(c(mean(wb3.rad.h.cc$ATC), mean(wb3.rad.h.cc$ATT)), 3),
                    ARBias = rep(bias(md3.true_eff, wb3.rad.h.cc), 3),
                    MRBias = rep(bias.med(md3.true_eff, wb3.rad.h.cc), 3),
                    RMSE = rep(RMSE(md3.true_eff, wb3.rad.h.cc), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(sim.aug.new.cc$ATC.wbvar), mean(sim.aug.new.cc$ATT.wbvar),
                                mean(wb3.rad.h.cc$ATC.wbvar), mean(wb3.rad.h.cc$ATT.wbvar),
                                mean(wb3.exp.h.cc$ATC.wbvar), mean(wb3.exp.h.cc$ATT.wbvar))),
                    SD.med = sqrt(c(median(sim.aug.new.cc$ATC.wbvar), median(sim.aug.new.cc$ATT.wbvar),
                                    median(wb3.rad.h.cc$ATC.wbvar), median(wb3.rad.h.cc$ATT.wbvar),
                                    median(wb3.exp.h.cc$ATC.wbvar), median(wb3.exp.h.cc$ATT.wbvar))),
                    ESD = rep(c(sd(wb3.rad.h.cc$ATC), sd(wb3.rad.h.cc$ATT)), 3),
                    RE = c(RE(sim.aug.new.cc), RE(wb3.rad.h.cc), RE(wb3.exp.h.cc)),
                    RE.med = c(RE.med(sim.aug.new.cc), RE.med(wb3.rad.h.cc), RE.med(wb3.exp.h.cc)),
                    CP = c(CP(md3.true_eff, sim.aug.new.cc),
                           CP(md3.true_eff, wb3.rad.h.cc), 
                           CP(md3.true_eff, wb3.exp.h.cc)))
print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(4, 6),
                    #Est = rep(c(mean(wb3.rad.c.cm$ATC), mean(wb3.rad.c.cm$ATT)), 3),
                    ARBias = rep(bias(truth.c, wb3.rad.c.cm), 3),
                    MRBias = rep(bias(truth.c, wb3.rad.c.cm), 3),
                    RMSE = rep(RMSE(truth.c, wb3.rad.c.cm), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(c.sim.aug.new.cm$ATC.wbvar), mean(c.sim.aug.new.cm$ATT.wbvar),
                                mean(wb3.rad.c.cm$ATC.wbvar), mean(wb3.rad.c.cm$ATT.wbvar),
                                mean(wb3.exp.c.cm$ATC.wbvar), mean(wb3.exp.c.cm$ATT.wbvar))),
                    SD.med = sqrt(c(median(c.sim.aug.new.cm$ATC.wbvar), median(c.sim.aug.new.cm$ATT.wbvar),
                                    median(wb3.rad.c.cm$ATC.wbvar), median(wb3.rad.c.cm$ATT.wbvar),
                                    median(wb3.exp.c.cm$ATC.wbvar), median(wb3.exp.c.cm$ATT.wbvar))),
                    ESD = rep(c(sd(wb3.rad.c.cm$ATC), sd(wb3.rad.c.cm$ATT)), 3),
                    RE = c(RE(c.sim.aug.new.cm), RE(wb3.rad.c.cm), RE(wb3.exp.c.cm)),
                    RE.med = c(RE.med(c.sim.aug.new.cm), RE.med(wb3.rad.c.cm), RE.med(wb3.exp.c.cm)),
                    CP = c(CP(truth.c, c.sim.aug.new.cm),
                           CP(truth.c, wb3.rad.c.cm), 
                           CP(truth.c, wb3.exp.c.cm)))
print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(as.numeric(md3.true_eff[, c("ATC", "ATT")]), 3),
                    #Est = rep(c(mean(wb3.rad.h.cm$ATC), mean(wb3.rad.h.cm$ATT)), 3),
                    ARBias = rep(bias(md3.true_eff, wb3.rad.h.cm), 3),
                    MRBias = rep(bias.med(md3.true_eff, wb3.rad.h.cm), 3),
                    RMSE = rep(RMSE(md3.true_eff, wb3.rad.h.cm), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(sim.aug.new.cm$ATC.wbvar), mean(sim.aug.new.cm$ATT.wbvar),
                                mean(wb3.rad.h.cm$ATC.wbvar), mean(wb3.rad.h.cm$ATT.wbvar),
                                mean(wb3.exp.h.cm$ATC.wbvar), mean(wb3.exp.h.cm$ATT.wbvar))),
                    SD.med = sqrt(c(median(sim.aug.new.cm$ATC.wbvar), median(sim.aug.new.cm$ATT.wbvar),
                                    median(wb3.rad.h.cm$ATC.wbvar), median(wb3.rad.h.cm$ATT.wbvar),
                                    median(wb3.exp.h.cm$ATC.wbvar), median(wb3.exp.h.cm$ATT.wbvar))),
                    ESD = rep(c(sd(wb3.rad.h.cm$ATC), sd(wb3.rad.h.cm$ATT)), 3),
                    RE = c(RE(sim.aug.new.cm), RE(wb3.rad.h.cm), RE(wb3.exp.h.cm)),
                    RE.med = c(RE.med(sim.aug.new.cm), RE.med(wb3.rad.h.cm), RE.med(wb3.exp.h.cm)),
                    CP = c(CP(md3.true_eff, sim.aug.new.cm),
                           CP(md3.true_eff, wb3.rad.h.cm), 
                           CP(md3.true_eff, wb3.exp.h.cm)))
print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(4, 6),
                    #Est = rep(c(mean(wb3.rad.c.mc$ATC), mean(wb3.rad.c.mc$ATT)), 3),
                    ARBias = rep(bias(truth.c, wb3.rad.c.mc), 3),
                    MRBias = rep(bias.med(truth.c, wb3.rad.c.mc), 3),
                    RMSE = rep(RMSE(truth.c, wb3.rad.c.mc), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(c.sim.aug.new.mc$ATC.wbvar), mean(c.sim.aug.new.mc$ATT.wbvar),
                                mean(wb3.rad.c.mc$ATC.wbvar), mean(wb3.rad.c.mc$ATT.wbvar),
                                mean(wb3.exp.c.mc$ATC.wbvar), mean(wb3.exp.c.mc$ATT.wbvar))),
                    SD.med = sqrt(c(median(c.sim.aug.new.mc$ATC.wbvar), median(c.sim.aug.new.mc$ATT.wbvar),
                                    median(wb3.rad.c.mc$ATC.wbvar), median(wb3.rad.c.mc$ATT.wbvar),
                                    median(wb3.exp.c.mc$ATC.wbvar), median(wb3.exp.c.mc$ATT.wbvar))),
                    ESD = rep(c(sd(wb3.rad.c.mc$ATC), sd(wb3.rad.c.mc$ATT)), 3),
                    RE = c(RE(c.sim.aug.new.mc), RE(wb3.rad.c.mc), RE(wb3.exp.c.mc)),
                    RE.med = c(RE.med(c.sim.aug.new.mc), RE.med(wb3.rad.c.mc), RE.med(wb3.exp.c.mc)),
                    CP = c(CP(truth.c, c.sim.aug.new.mc),
                           CP(truth.c, wb3.rad.c.mc), 
                           CP(truth.c, wb3.exp.c.mc)))
print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(as.numeric(md3.true_eff[, c("ATC", "ATT")]), 3),
                    #Est = rep(c(mean(wb3.rad.h.mc$ATC), mean(wb3.rad.h.mc$ATT)), 3),
                    ARBias = rep(bias(md3.true_eff, wb3.rad.h.mc), 3),
                    MRBias = rep(bias.med(md3.true_eff, wb3.rad.h.mc), 3),
                    RMSE = rep(RMSE(md3.true_eff, wb3.rad.h.mc), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(sim.aug.new.mc$ATC.wbvar), mean(sim.aug.new.mc$ATT.wbvar),
                                mean(wb3.rad.h.mc$ATC.wbvar), mean(wb3.rad.h.mc$ATT.wbvar),
                                mean(wb3.exp.h.mc$ATC.wbvar), mean(wb3.exp.h.mc$ATT.wbvar))),
                    SD.med = sqrt(c(median(sim.aug.new.mc$ATC.wbvar), median(sim.aug.new.mc$ATT.wbvar),
                                    median(wb3.rad.h.mc$ATC.wbvar), median(wb3.rad.h.mc$ATT.wbvar),
                                    median(wb3.exp.h.mc$ATC.wbvar), median(wb3.exp.h.mc$ATT.wbvar))),
                    ESD = rep(c(sd(wb3.rad.h.mc$ATC), sd(wb3.rad.h.mc$ATT)), 3),
                    RE = c(RE(sim.aug.new.mc), RE(wb3.rad.h.mc), RE(wb3.exp.h.mc)),
                    RE.med = c(RE.med(sim.aug.new.mc), RE.med(wb3.rad.h.mc), RE.med(wb3.exp.h.mc)),
                    CP = c(CP(md3.true_eff, sim.aug.new.mc),
                           CP(md3.true_eff, wb3.rad.h.mc), 
                           CP(md3.true_eff, wb3.exp.h.mc)))
print(xtable(table), include.rownames=FALSE)

table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(4, 6),
                    #Est = rep(c(mean(wb3.rad.c.mm$ATC), mean(wb3.rad.c.mm$ATT)), 3),
                    ARBias = rep(bias(truth.c, wb3.rad.c.mm), 3),
                    MRBias = rep(bias.med(truth.c, wb3.rad.c.mm), 3),
                    RMSE = rep(RMSE(truth.c, wb3.rad.c.mm), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(c.sim.aug.new.mm$ATC.wbvar), mean(c.sim.aug.new.mm$ATT.wbvar),
                                mean(wb3.rad.c.mm$ATC.wbvar), mean(wb3.rad.c.mm$ATT.wbvar),
                                mean(wb3.exp.c.mm$ATC.wbvar), mean(wb3.exp.c.mm$ATT.wbvar))),
                    SD.med = sqrt(c(median(c.sim.aug.new.mm$ATC.wbvar), median(c.sim.aug.new.mm$ATT.wbvar),
                                    median(wb3.rad.c.mm$ATC.wbvar), median(wb3.rad.c.mm$ATT.wbvar),
                                    median(wb3.exp.c.mm$ATC.wbvar), median(wb3.exp.c.mm$ATT.wbvar))),
                    ESD = rep(c(sd(wb3.rad.c.mm$ATC), sd(wb3.rad.c.mm$ATT)), 3),
                    RE = c(RE(c.sim.aug.new.mm), RE(wb3.rad.c.mm), RE(wb3.exp.c.mm)),
                    RE.med = c(RE.med(c.sim.aug.new.mm), RE.med(wb3.rad.c.mm), RE.med(wb3.exp.c.mm)),
                    CP = c(CP(truth.c, c.sim.aug.new.mm),
                           CP(truth.c, wb3.rad.c.mm), 
                           CP(truth.c, wb3.exp.c.mm)))
print(xtable(table), include.rownames=FALSE)


table <- data.frame(Est. = rep(c("ATC", "ATT"), 3),
                    #Truth = rep(as.numeric(md3.true_eff[, c("ATC", "ATT")]), 3),
                    #Est = rep(c(mean(wb3.rad.h.mm$ATC), mean(wb3.rad.h.mm$ATT)), 3),
                    ARBias = rep(bias(md3.true_eff, wb3.rad.h.mm), 3),
                    MRBias = rep(bias.med(md3.true_eff, wb3.rad.h.mm), 3),
                    RMSE = rep(RMSE(md3.true_eff, wb3.rad.h.mm), 3),
                    Method = c("Sandwich", "Sandwich", "WB (Rad)", "WB (Rad)", 
                               "WB (Exp)", "WB (Exp)"),
                    SD = sqrt(c(mean(sim.aug.new.mm$ATC.wbvar), mean(sim.aug.new.mm$ATT.wbvar),
                                mean(wb3.rad.h.mm$ATC.wbvar), mean(wb3.rad.h.mm$ATT.wbvar),
                                mean(wb3.exp.h.mm$ATC.wbvar), mean(wb3.exp.h.mm$ATT.wbvar))),
                    SD.med = sqrt(c(median(sim.aug.new.mm$ATC.wbvar), median(sim.aug.new.mm$ATT.wbvar),
                                    median(wb3.rad.h.mm$ATC.wbvar), median(wb3.rad.h.mm$ATT.wbvar),
                                    median(wb3.exp.h.mm$ATC.wbvar), median(wb3.exp.h.mm$ATT.wbvar))),
                    ESD = rep(c(sd(wb3.rad.h.mm$ATC), sd(wb3.rad.h.mm$ATT)), 3),
                    RE = c(RE(sim.aug.new.mm), RE(wb3.rad.h.mm), RE(wb3.exp.h.mm)),
                    RE.med = c(RE.med(sim.aug.new.mm), RE.med(wb3.rad.h.mm), RE.med(wb3.exp.h.mm)),
                    CP = c(CP(md3.true_eff, sim.aug.new.mm),
                           CP(md3.true_eff, wb3.rad.h.mm), 
                           CP(md3.true_eff, wb3.exp.h.mm)))
print(xtable(table), include.rownames=FALSE)







