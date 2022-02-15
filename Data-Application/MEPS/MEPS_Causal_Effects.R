### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~ PS Methods Expenditure Data Application ~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### Data cleaning and causal estimates

### by Yi Liu

library(PSweight)
library(dplyr)

load("meps.Rdata")
meps$healthExp <- as.numeric(meps$healthExp)
cov.names = c("pcs", "mcs", "age","bmi","sinceCheckUp","male","married", "h.excellent", 
              "h.verygood", "h.good", "h.fair", "h.poor", "mh.excellent", "mh.verygood", 
              "mh.good", "mh.fair", "mh.poor","anylim", "exercise", "hibp", "chd", "angina",
              "mi", "stroke", "emphysema", "cholesterol", "cancer", "diabetes", "arthritis", 
              "asthma", "smoke")
find.var.num <- function(varname, dat) for (i in 1:dim(dat)[2]) if (names(dat)[i]==varname) return(i)

### Three data with different p
dat1 <- meps[meps$hispanic == 1 | meps$white == 1, 
             unlist(sapply(c("white", "healthExp", cov.names), find.var.num, dat = meps))]
table(dat1$white)

dat2 <- meps[meps$black == 1 | meps$white == 1, 
             unlist(sapply(c("white", "healthExp", cov.names), find.var.num, dat = meps))]
table(dat2$white)

dat3 <- meps[meps$asian == 1 | meps$white == 1, 
             unlist(sapply(c("white", "healthExp", cov.names), find.var.num, dat = meps))]
table(dat3$white)

# Design module
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

# Analysis module
source("newSand_func.R")

## Hajek-like estimator
### White-Hispanic
y <- dat1$healthExp
z <- dat1$white
X <- dat1[,cov.names] %>% as.matrix()

ATE.dat1  <- ATE(y=y, z=z, X=X)
ATO.dat1  <- ATO(y=y, z=z, X=X)
ATM.dat1  <- ATM(y=y, z=z, X=X)
ATEN.dat1 <- ATEN(y=y, z=z, X=X)
ATC.dat1  <- ATC(y=y, z=z, X=X)
ATT.dat1  <- ATT(y=y, z=z, X=X)

trim1 <- PStrim(data = dat1, ps.formula = ps.mult, delta = 0.05)$data
y <- trim1$healthExp
z <- trim1$white
X <- trim1[,cov.names] %>% as.matrix()
ATE.dat1.05 <- ATE(y=y, z=z, X=X)

trim2 <- PStrim(data = dat1, ps.formula = ps.mult, delta = 0.10)$data
y <- trim2$healthExp
z <- trim2$white
X <- trim2[,cov.names] %>% as.matrix()
ATE.dat1.10 <- ATE(y=y, z=z, X=X)

trim3 <- PStrim(data = dat1, ps.formula = ps.mult, delta = 0.15)$data
y <- trim3$healthExp
z <- trim3$white
X <- trim3[,cov.names] %>% as.matrix()
ATE.dat1.15 <- ATE(y=y, z=z, X=X)

### White-Black
y <- dat2$healthExp
z <- dat2$white
X <- dat2[,cov.names] %>% as.matrix()

ATE.dat2  <- ATE(y=y, z=z, X=X)
ATO.dat2  <- ATO(y=y, z=z, X=X)
ATM.dat2  <- ATM(y=y, z=z, X=X)
ATEN.dat2 <- ATEN(y=y, z=z, X=X)
ATC.dat2  <- ATC(y=y, z=z, X=X)
ATT.dat2  <- ATT(y=y, z=z, X=X)

trim1 <- PStrim(data = dat2, ps.formula = ps.mult, delta = 0.05)$data
y <- trim1$healthExp
z <- trim1$white
X <- trim1[,cov.names] %>% as.matrix()
ATE.dat2.05 <- ATE(y=y, z=z, X=X)

trim2 <- PStrim(data = dat2, ps.formula = ps.mult, delta = 0.10)$data
y <- trim2$healthExp
z <- trim2$white
X <- trim2[,cov.names] %>% as.matrix()
ATE.dat2.10 <- ATE(y=y, z=z, X=X)

trim3 <- PStrim(data = dat2, ps.formula = ps.mult, delta = 0.15)$data
y <- trim3$healthExp
z <- trim3$white
X <- trim3[,cov.names] %>% as.matrix()
ATE.dat2.15 <- ATE(y=y, z=z, X=X)

### White-Asian
y <- dat3$healthExp
z <- dat3$white
X <- dat3[,cov.names] %>% as.matrix()

ATE.dat3  <- ATE(y=y, z=z, X=X)
ATO.dat3  <- ATO(y=y, z=z, X=X)
ATM.dat3  <- ATM(y=y, z=z, X=X)
ATEN.dat3 <- ATEN(y=y, z=z, X=X)
ATC.dat3  <- ATC(y=y, z=z, X=X)
ATT.dat3  <- ATT(y=y, z=z, X=X)

trim1 <- PStrim(data = dat3, ps.formula = ps.mult, delta = 0.05)$data
y <- trim1$healthExp
z <- trim1$white
X <- trim1[,cov.names] %>% as.matrix()
ATE.dat3.05 <- ATE(y=y, z=z, X=X)

trim2 <- PStrim(data = dat3, ps.formula = ps.mult, delta = 0.10)$data
y <- trim2$healthExp
z <- trim2$white
X <- trim2[,cov.names] %>% as.matrix()
ATE.dat3.10 <- ATE(y=y, z=z, X=X)

trim3 <- PStrim(data = dat3, ps.formula = ps.mult, delta = 0.15)$data
y <- trim3$healthExp
z <- trim3$white
X <- trim3[,cov.names] %>% as.matrix()
ATE.dat3.15 <- ATE(y=y, z=z, X=X)

## Augmented estimator
### White-Hispanic
y <- dat1$healthExp
z <- dat1$white
X <- dat1[,cov.names] %>% as.matrix()

ATE.dat1.aug  <- ATE(y=y, z=z, X=X, DR=TRUE, X.out=X)
ATO.dat1.aug  <- ATO(y=y, z=z, X=X, DR=TRUE, X.out=X)
ATM.dat1.aug  <- ATM(y=y, z=z, X=X, DR=TRUE, X.out=X)
ATEN.dat1.aug <- ATEN(y=y, z=z, X=X, DR=TRUE, X.out=X)
ATC.dat1.aug  <- ATC(y=y, z=z, X=X, DR=TRUE, X.out=X)
ATT.dat1.aug  <- ATT(y=y, z=z, X=X, DR=TRUE, X.out=X)

trim1 <- PStrim(data = dat1, ps.formula = ps.mult, delta = 0.05)$data
y <- trim1$healthExp
z <- trim1$white
X <- trim1[,cov.names] %>% as.matrix()
ATE.dat1.05.aug <- ATE(y=y, z=z, X=X, DR=TRUE, X.out=X)

trim2 <- PStrim(data = dat1, ps.formula = ps.mult, delta = 0.10)$data
y <- trim2$healthExp
z <- trim2$white
X <- trim2[,cov.names] %>% as.matrix()
ATE.dat1.10.aug <- ATE(y=y, z=z, X=X, DR=TRUE, X.out=X)

trim3 <- PStrim(data = dat1, ps.formula = ps.mult, delta = 0.15)$data
y <- trim3$healthExp
z <- trim3$white
X <- trim3[,cov.names] %>% as.matrix()
ATE.dat1.15.aug <- ATE(y=y, z=z, X=X, DR=TRUE, X.out=X)

### White-Black
y <- dat2$healthExp
z <- dat2$white
X <- dat2[,cov.names] %>% as.matrix()

ATE.dat2.aug  <- ATE(y=y, z=z, X=X, DR=TRUE, X.out=X)
ATO.dat2.aug  <- ATO(y=y, z=z, X=X, DR=TRUE, X.out=X)
ATM.dat2.aug  <- ATM(y=y, z=z, X=X, DR=TRUE, X.out=X)
ATEN.dat2.aug <- ATEN(y=y, z=z, X=X, DR=TRUE, X.out=X)
ATC.dat2.aug  <- ATC(y=y, z=z, X=X, DR=TRUE, X.out=X)
ATT.dat2.aug  <- ATT(y=y, z=z, X=X, DR=TRUE, X.out=X)

trim1 <- PStrim(data = dat2, ps.formula = ps.mult, delta = 0.05)$data
y <- trim1$healthExp
z <- trim1$white
X <- trim1[,cov.names] %>% as.matrix()
ATE.dat2.05.aug <- ATE(y=y, z=z, X=X, DR=TRUE, X.out=X)

trim2 <- PStrim(data = dat2, ps.formula = ps.mult, delta = 0.10)$data
y <- trim2$healthExp
z <- trim2$white
X <- trim2[,cov.names] %>% as.matrix()
ATE.dat2.10.aug <- ATE(y=y, z=z, X=X, DR=TRUE, X.out=X)

trim3 <- PStrim(data = dat2, ps.formula = ps.mult, delta = 0.15)$data
y <- trim3$healthExp
z <- trim3$white
X <- trim3[,cov.names] %>% as.matrix()
ATE.dat2.15.aug <- ATE(y=y, z=z, X=X, DR=TRUE, X.out=X)

### White-Asian
y <- dat3$healthExp
z <- dat3$white
X <- dat3[,cov.names] %>% as.matrix()

ATE.dat3.aug  <- ATE(y=y, z=z, X=X, DR=TRUE, X.out=X)
ATO.dat3.aug  <- ATO(y=y, z=z, X=X, DR=TRUE, X.out=X)
ATM.dat3.aug  <- ATM(y=y, z=z, X=X, DR=TRUE, X.out=X)
ATEN.dat3.aug <- ATEN(y=y, z=z, X=X, DR=TRUE, X.out=X)
ATC.dat3.aug  <- ATC(y=y, z=z, X=X, DR=TRUE, X.out=X)
ATT.dat3.aug  <- ATT(y=y, z=z, X=X, DR=TRUE, X.out=X)

trim1 <- PStrim(data = dat3, ps.formula = ps.mult, delta = 0.05)$data
y <- trim1$healthExp
z <- trim1$white
X <- trim1[,cov.names] %>% as.matrix()
ATE.dat3.05.aug <- ATE(y=y, z=z, X=X, DR=TRUE, X.out=X)

trim2 <- PStrim(data = dat3, ps.formula = ps.mult, delta = 0.10)$data
y <- trim2$healthExp
z <- trim2$white
X <- trim2[,cov.names] %>% as.matrix()
ATE.dat3.10.aug <- ATE(y=y, z=z, X=X, DR=TRUE, X.out=X)

trim3 <- PStrim(data = dat3, ps.formula = ps.mult, delta = 0.15)$data
y <- trim3$healthExp
z <- trim3$white
X <- trim3[,cov.names] %>% as.matrix()
ATE.dat3.15.aug <- ATE(y=y, z=z, X=X, DR=TRUE, X.out=X)

# Save the causal effects data
save(file = "MEPS_Causal_Effects.RData", 
     ATE.dat1,
     ATE.dat1.05,
     ATE.dat1.10,
     ATE.dat1.15,
     ATO.dat1,
     ATM.dat1,
     ATEN.dat1,
     ATC.dat1,
     ATT.dat1,
     
     ATE.dat2,
     ATE.dat2.05,
     ATE.dat2.10,
     ATE.dat2.15,
     ATO.dat2,
     ATM.dat2,
     ATEN.dat2,
     ATC.dat2,
     ATT.dat2,
     
     ATE.dat3,
     ATE.dat3.05,
     ATE.dat3.10,
     ATE.dat3.15,
     ATO.dat3,
     ATM.dat3,
     ATEN.dat3,
     ATC.dat3,
     ATT.dat3,
     
     ATE.dat1.aug,
     ATE.dat1.05.aug,
     ATE.dat1.10.aug,
     ATE.dat1.15.aug,
     ATO.dat1.aug,
     ATM.dat1.aug,
     ATEN.dat1.aug,
     ATC.dat1.aug,
     ATT.dat1.aug,
     
     ATE.dat2.aug,
     ATE.dat2.05.aug,
     ATE.dat2.10.aug,
     ATE.dat2.15.aug,
     ATO.dat2.aug,
     ATM.dat2.aug,
     ATEN.dat2.aug,
     ATC.dat2.aug,
     ATT.dat2.aug,
     
     ATE.dat3.aug,
     ATE.dat3.05.aug,
     ATE.dat3.10.aug,
     ATE.dat3.15.aug,
     ATO.dat3.aug,
     ATM.dat3.aug,
     ATEN.dat3.aug,
     ATC.dat3.aug,
     ATT.dat3.aug)
