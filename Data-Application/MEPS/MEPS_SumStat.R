### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~ PS Methods Expenditure Data Application ~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### Summary statistics

### by Yi Liu

library(PSweight)
library(dplyr)
library(xtable)

load("meps.RData")
meps$healthExp <- as.numeric(meps$healthExp)
load("MEPS_Causal_Effects.RData")

cov.names = c("pcs", "mcs", "age","bmi","sinceCheckUp","male","married", "h.excellent", 
              "h.verygood", "h.good", "h.fair", "h.poor", "mh.excellent", "mh.verygood", 
              "mh.good", "mh.fair", "mh.poor","anylim", "exercise", "hibp", "chd", "angina",
              "mi", "stroke", "emphysema", "cholesterol", "cancer", "diabetes", "arthritis", 
              "asthma", "smoke")
find.var.num <- function(varname, dat) for (i in 1:dim(dat)[2]) if (names(dat)[i]==varname) return(i)

dat1 <- meps[meps$hispanic == 1 | meps$white == 1, 
             unlist(sapply(c("white", "healthExp", cov.names), find.var.num, dat = meps))]
dat2 <- meps[meps$black == 1 | meps$white == 1, 
             unlist(sapply(c("white", "healthExp", cov.names), find.var.num, dat = meps))]
dat3 <- meps[meps$asian == 1 | meps$white == 1, 
             unlist(sapply(c("white", "healthExp", cov.names), find.var.num, dat = meps))]

ps.mult <- white ~ pcs + mcs + age + bmi + sinceCheckUp + male + married + h.excellent + h.verygood + 
  h.good + h.fair + h.poor + mh.excellent + mh.verygood + mh.good + mh.fair + mh.poor + anylim +
  exercise + hibp + chd + angina + mi + stroke + emphysema + cholesterol + cancer + diabetes + 
  arthritis + asthma + smoke

out.form <- healthExp ~ pcs + mcs + age + bmi + sinceCheckUp + male + married + h.excellent + h.verygood + 
  h.good + h.fair + h.poor + mh.excellent + mh.verygood + mh.good + mh.fair + mh.poor + anylim +
  exercise + hibp + chd + angina + mi + stroke + emphysema + cholesterol + cancer + diabetes + 
  arthritis + asthma + smoke

### Propensity score plots
bal.mult <- SumStat(ps.formula = ps.mult, data = dat1, weight = "IPW")
png("ps_meps_hispanic.png", res=72*2, width = 1200, height = 800)
plot(bal.mult, type = "hist")
dev.off()

bal.mult <- SumStat(ps.formula = ps.mult, data = dat2, weight = "IPW")
png("ps_meps_black.png", res=72*2, width = 1200, height = 800)
plot(bal.mult, type = "hist")
dev.off()

bal.mult <- SumStat(ps.formula = ps.mult, data = dat3, weight = "IPW")
png("ps_meps_asian.png", res=72*2, width = 1200, height = 800)
plot(bal.mult, type = "hist")
dev.off()

## Love plots
dat1$nonwhite <- as.numeric(dat1$white==0)
dat2$nonwhite <- as.numeric(dat2$white==0)
dat3$nonwhite <- as.numeric(dat3$white==0)
atc.mult <- nonwhite ~ pcs + mcs + age + bmi + sinceCheckUp + male + married + h.excellent + h.verygood + 
  h.good + h.fair + h.poor + mh.excellent + mh.verygood + mh.good + mh.fair + mh.poor + anylim +
  exercise + hibp + chd + angina + mi + stroke + emphysema + cholesterol + cancer + diabetes + 
  arthritis + asthma + smoke

### White-Hispanic
trim1 <- PStrim(data = dat1, ps.formula = ps.mult, delta = 0.05)
trim2 <- PStrim(data = dat1, ps.formula = ps.mult, delta = 0.1)
trim3 <- PStrim(data = dat1, ps.formula = ps.mult, delta = 0.15)
covar <- names(SumStat(ps.formula = ps.mult, weight = "IPW", data = dat1)$IPW.sumstat[, "ASD weighted var 1-2"])

bal.IPW <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "IPW", data = dat1)$IPW.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATE", covar = covar)
colnames(bal.IPW) <- c("ASD", "Method", "covar")
bal.IPW.5 <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "IPW", data = trim1$data)$IPW.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATE (0.05)", covar = covar)
colnames(bal.IPW.5) <- c("ASD", "Method", "covar")
bal.IPW.10 <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "IPW", data = trim2$data)$IPW.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATE (0.1)", covar = covar)
colnames(bal.IPW.10) <- c("ASD", "Method", "covar")
bal.IPW.15 <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "IPW", data = trim3$data)$IPW.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATE (0.15)", covar = covar)
colnames(bal.IPW.15) <- c("ASD", "Method", "covar")
bal.ATO <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "overlap", data = dat1)$overlap.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATO", covar = covar)
colnames(bal.ATO) <- c("ASD", "Method", "covar")
bal.ATM <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "matching", data = dat1)$matching.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATM", covar = covar)
colnames(bal.ATM) <- c("ASD", "Method", "covar")
bal.ATEN <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "entropy", data = dat1)$entropy.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATEN", covar = covar)
colnames(bal.ATEN) <- c("ASD", "Method", "covar")
bal.ATT <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "treated", data = dat1)$treated.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATT", covar = covar)
colnames(bal.ATT) <- c("ASD", "Method", "covar")
bal.ATC <- as.data.frame(SumStat(ps.formula = atc.mult, weight = "treated", data = dat1)$treated.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATC", covar = covar)
colnames(bal.ATC) <- c("ASD", "Method", "covar")

df <- rbind(bal.IPW, bal.IPW.5, bal.IPW.10, bal.IPW.15, bal.ATO, bal.ATM, bal.ATEN, bal.ATT, bal.ATC)
df$Method <- factor(df$Method, levels = c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", "ATO", "ATM", "ATEN", "ATT", "ATC"))

png("asd_hispanic.png", res = 200, width = 1000, height = 800)
ggplot(df, aes(x = ASD, y = covar, col = Method)) + geom_point(alpha = 0.85) + 
  labs(x = "Standarized Mean Difference", y = "Covariates") + 
  scale_color_manual(values = c("royalblue", "deepskyblue", "darkseagreen1", 
                                "forestgreen", "sienna", "indianred2", "plum2", 
                                "goldenrod2", "dimgrey")) + 
  geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed", color = "black", size = 0.2) + 
  theme_bw()
dev.off()

n <- nrow(dat1)
n1 <- nrow(trim1$data) # sample size after trimming
n2 <- nrow(trim2$data)
n3 <- nrow(trim3$data)
n1/n # percentage of sample size
n2/n
n3/n

### White-Black
trim1 <- PStrim(data = dat2, ps.formula = ps.mult, delta = 0.05)
trim2 <- PStrim(data = dat2, ps.formula = ps.mult, delta = 0.1)
trim3 <- PStrim(data = dat2, ps.formula = ps.mult, delta = 0.15)
covar <- names(SumStat(ps.formula = ps.mult, weight = "IPW", data = dat2)$IPW.sumstat[, "ASD weighted var 1-2"])

bal.IPW <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "IPW", data = dat2)$IPW.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATE", covar = covar)
colnames(bal.IPW) <- c("ASD", "Method", "covar")
bal.IPW.5 <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "IPW", data = trim1$data)$IPW.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATE (0.05)", covar = covar)
colnames(bal.IPW.5) <- c("ASD", "Method", "covar")
bal.IPW.10 <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "IPW", data = trim2$data)$IPW.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATE (0.1)", covar = covar)
colnames(bal.IPW.10) <- c("ASD", "Method", "covar")
bal.IPW.15 <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "IPW", data = trim3$data)$IPW.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATE (0.15)", covar = covar)
colnames(bal.IPW.15) <- c("ASD", "Method", "covar")
bal.ATO <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "overlap", data = dat2)$overlap.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATO", covar = covar)
colnames(bal.ATO) <- c("ASD", "Method", "covar")
bal.ATM <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "matching", data = dat2)$matching.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATM", covar = covar)
colnames(bal.ATM) <- c("ASD", "Method", "covar")
bal.ATEN <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "entropy", data = dat2)$entropy.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATEN", covar = covar)
colnames(bal.ATEN) <- c("ASD", "Method", "covar")
bal.ATT <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "treated", data = dat2)$treated.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATT", covar = covar)
colnames(bal.ATT) <- c("ASD", "Method", "covar")
bal.ATC <- as.data.frame(SumStat(ps.formula = atc.mult, weight = "treated", data = dat2)$treated.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATC", covar = covar)
colnames(bal.ATC) <- c("ASD", "Method", "covar")

df <- rbind(bal.IPW, bal.IPW.5, bal.IPW.10, bal.IPW.15, bal.ATO, bal.ATM, bal.ATEN, bal.ATT, bal.ATC)
df$Method <- factor(df$Method, levels = c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", "ATO", "ATM", "ATEN", "ATT", "ATC"))

png("asd_black.png", res = 200, width = 1000, height = 800)
ggplot(df, aes(x = ASD, y = covar, col = Method)) + geom_point(alpha = 0.85) + 
  labs(x = "Standarized Mean Difference", y = "Covariates") + 
  scale_color_manual(values = c("royalblue", "deepskyblue", "darkseagreen1", 
                                "forestgreen", "sienna", "indianred2", "plum2", 
                                "goldenrod2", "dimgrey")) + 
  geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed", color = "black", size = 0.2) + 
  theme_bw()
dev.off()

n <- nrow(dat2)
n1 <- nrow(trim1$data) # sample size after trimming
n2 <- nrow(trim2$data)
n3 <- nrow(trim3$data)
n1/n # percentage of sample size
n2/n
n3/n

### White-Asian
trim1 <- PStrim(data = dat3, ps.formula = ps.mult, delta = 0.05)
trim2 <- PStrim(data = dat3, ps.formula = ps.mult, delta = 0.1)
trim3 <- PStrim(data = dat3, ps.formula = ps.mult, delta = 0.15)
covar <- names(SumStat(ps.formula = ps.mult, weight = "IPW", data = dat3)$IPW.sumstat[, "ASD weighted var 1-2"])

bal.IPW <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "IPW", data = dat3)$IPW.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATE", covar = covar)
colnames(bal.IPW) <- c("ASD", "Method", "covar")
bal.IPW.5 <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "IPW", data = trim1$data)$IPW.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATE (0.05)", covar = covar)
colnames(bal.IPW.5) <- c("ASD", "Method", "covar")
bal.IPW.10 <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "IPW", data = trim2$data)$IPW.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATE (0.1)", covar = covar)
colnames(bal.IPW.10) <- c("ASD", "Method", "covar")
bal.IPW.15 <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "IPW", data = trim3$data)$IPW.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATE (0.15)", covar = covar)
colnames(bal.IPW.15) <- c("ASD", "Method", "covar")
bal.ATO <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "overlap", data = dat3)$overlap.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATO", covar = covar)
colnames(bal.ATO) <- c("ASD", "Method", "covar")
bal.ATM <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "matching", data = dat3)$matching.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATM", covar = covar)
colnames(bal.ATM) <- c("ASD", "Method", "covar")
bal.ATEN <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "entropy", data = dat3)$entropy.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATEN", covar = covar)
colnames(bal.ATEN) <- c("ASD", "Method", "covar")
bal.ATT <- as.data.frame(SumStat(ps.formula = ps.mult, weight = "treated", data = dat3)$treated.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATT", covar = covar)
colnames(bal.ATT) <- c("ASD", "Method", "covar")
bal.ATC <- as.data.frame(SumStat(ps.formula = atc.mult, weight = "treated", data = dat3)$treated.sumstat[, "ASD weighted var 1-2"]) %>% mutate(type = "ATC", covar = covar)
colnames(bal.ATC) <- c("ASD", "Method", "covar")

df <- rbind(bal.IPW, bal.IPW.5, bal.IPW.10, bal.IPW.15, bal.ATO, bal.ATM, bal.ATEN, bal.ATT, bal.ATC)
df$Method <- factor(df$Method, levels = c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", "ATO", "ATM", "ATEN", "ATT", "ATC"))

png("asd_asian.png", res = 200, width = 1000, height = 800)
ggplot(df, aes(x = ASD, y = covar, col = Method)) + geom_point(alpha = 0.85) + 
  labs(x = "Standarized Mean Difference", y = "Covariates") + 
  scale_color_manual(values = c("royalblue", "deepskyblue", "darkseagreen1", 
                                "forestgreen", "sienna", "indianred2", "plum2", 
                                "goldenrod2", "dimgrey")) + 
  geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed", color = "black", size = 0.2) + 
  theme_bw()
dev.off()

n <- nrow(dat3)
n1 <- nrow(trim1$data) # sample size after trimming
n2 <- nrow(trim2$data)
n3 <- nrow(trim3$data)
n1/n # percentage of sample size
n2/n
n3/n

## Ratio of propensity score variances in two groups
### White-Hispanic
sum(dat1$white)/nrow(dat1) # p
fit <- glm(ps.mult, family = binomial(link = "logit"), data = dat1)
e.h <- as.numeric(fit$fitted.values)
var(e.h[dat1$white==1])/var(e.h[dat1$white==0]) # ratio

### White-Black
sum(dat1$white)/nrow(dat2) # p
fit <- glm(ps.mult, family = binomial(link = "logit"), data = dat2)
e.h <- as.numeric(fit$fitted.values)
var(e.h[dat2$white==1])/var(e.h[dat2$white==0]) # ratio

### White-Asian
sum(dat1$white)/nrow(dat3) # p
fit <- glm(ps.mult, family = binomial(link = "logit"), data = dat3)
e.h <- as.numeric(fit$fitted.values)
var(e.h[dat3$white==1])/var(e.h[dat3$white==0]) # ratio

## Summarizing the causal effects table
rows <- c("ATE", "ATE (0.05)", "ATE (0.1)", "ATE (0.15)", "ATO", "ATM", "ATEN", "ATT", "ATC")
est1.wt  <- c(ATE.dat1$tau, ATE.dat1.05$tau, ATE.dat1.10$tau, ATE.dat1.15$tau, 
              ATO.dat1$tau, ATM.dat1$tau, ATEN.dat1$tau, ATT.dat1$tau, ATC.dat1$tau)
est1.aug <- c(ATE.dat1.aug$tau, ATE.dat1.05.aug$tau, ATE.dat1.10.aug$tau, ATE.dat1.15.aug$tau, 
              ATO.dat1.aug$tau, ATM.dat1.aug$tau, ATEN.dat1.aug$tau, ATT.dat1.aug$tau, ATC.dat1.aug$tau)
var1.wt  <- c(ATE.dat1$se, ATE.dat1.05$se, ATE.dat1.10$se, ATE.dat1.15$se, 
              ATO.dat1$se, ATM.dat1$se, ATEN.dat1$se, ATT.dat1$se, ATC.dat1$se)
var1.aug <- c(ATE.dat1.aug$se, ATE.dat1.05.aug$se, ATE.dat1.10.aug$se, ATE.dat1.15.aug$se, 
              ATO.dat1.aug$se, ATM.dat1.aug$se, ATEN.dat1.aug$se, ATT.dat1.aug$se, ATC.dat1.aug$se)

est2.wt  <- c(ATE.dat2$tau, ATE.dat2.05$tau, ATE.dat2.10$tau, ATE.dat2.15$tau, 
              ATO.dat2$tau, ATM.dat2$tau, ATEN.dat2$tau, ATT.dat2$tau, ATC.dat2$tau)
est2.aug <- c(ATE.dat2.aug$tau, ATE.dat2.05.aug$tau, ATE.dat2.10.aug$tau, ATE.dat2.15.aug$tau, 
              ATO.dat2.aug$tau, ATM.dat2.aug$tau, ATEN.dat2.aug$tau, ATT.dat2.aug$tau, ATC.dat2.aug$tau)
var2.wt  <- c(ATE.dat2$se, ATE.dat2.05$se, ATE.dat2.10$se, ATE.dat2.15$se, 
              ATO.dat2$se, ATM.dat2$se, ATEN.dat2$se, ATT.dat2$se, ATC.dat2$se)
var2.aug <- c(ATE.dat2.aug$se, ATE.dat2.05.aug$se, ATE.dat2.10.aug$se, ATE.dat2.15.aug$se, 
              ATO.dat2.aug$se, ATM.dat2.aug$se, ATEN.dat2.aug$se, ATT.dat2.aug$se, ATC.dat2.aug$se)

est3.wt  <- c(ATE.dat3$tau, ATE.dat3.05$tau, ATE.dat3.10$tau, ATE.dat3.15$tau, 
              ATO.dat3$tau, ATM.dat3$tau, ATEN.dat3$tau, ATT.dat3$tau, ATC.dat3$tau)
est3.aug <- c(ATE.dat3.aug$tau, ATE.dat3.05.aug$tau, ATE.dat3.10.aug$tau, ATE.dat3.15.aug$tau, 
              ATO.dat3.aug$tau, ATM.dat3.aug$tau, ATEN.dat3.aug$tau, ATT.dat3.aug$tau, ATC.dat3.aug$tau)
var3.wt  <- c(ATE.dat3$se, ATE.dat3.05$se, ATE.dat3.10$se, ATE.dat3.15$se, 
              ATO.dat3$se, ATM.dat3$se, ATEN.dat3$se, ATT.dat3$se, ATC.dat3$se)
var3.aug <- c(ATE.dat3.aug$se, ATE.dat3.05.aug$se, ATE.dat3.10.aug$se, ATE.dat3.15.aug$se, 
              ATO.dat3.aug$se, ATM.dat3.aug$se, ATEN.dat3.aug$se, ATT.dat3.aug$se, ATC.dat3.aug$se)

p.value  <- function(x) ifelse(x<0.001, "<0.001", as.character(x))
SumTab1  <- data.frame(Estimand = rows, Est.wt = est1.wt, SE.wt = var1.wt, p.wt = p.value(1-2*pnorm(est1.wt/var1.wt)),
                       Est.aug = est1.aug, SE.aug = var1.aug, p.aug = p.value(1-2*pnorm(est1.aug/var1.aug)))
SumTab2  <- data.frame(Estimand = rows, Est.wt = est2.wt, SE.wt = var2.wt, p.wt = p.value(1-2*pnorm(est2.wt/var2.wt)),
                       Est.aug = est2.aug, SE.aug = var2.aug, p.aug = p.value(1-2*pnorm(est2.aug/var2.aug)))
SumTab3  <- data.frame(Estimand = rows, Est.wt = est3.wt, SE.wt = var3.wt, p.wt = p.value(1-2*pnorm(est3.wt/var3.wt)),
                       Est.aug = est3.aug, SE.aug = var3.aug, p.aug = p.value(1-2*pnorm(est3.aug/var3.aug)))

print(xtable(SumTab1, digits = 2), include.rownames=FALSE)
print(xtable(SumTab2, digits = 2), include.rownames=FALSE)
print(xtable(SumTab3, digits = 2), include.rownames=FALSE)
