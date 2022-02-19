### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

####### Toy simulation for introduction

### by Yi Liu

library(dplyr)
library(PSweight)
library(xtable)

### Data generating process (DGP)
n <- 10000
X <- cbind(X1=rnorm(n=n, 6, 3), X2=rbinom(n=n, 1, 0.75)) # covariates

sigmoid = function(x) 1/(1+exp(-x))
XD <- cbind(rep(1,n), X) # design matrix

### Simulation 1: large p
alpha1 <- c(-3.5, 0.2, 0.8) # coefficients
z1 <- rbinom(n,1,prob=sigmoid(alpha1 %*% t(XD))) # treatment assignments
p1 <- mean(z1) # p
y1 <- as.numeric(c(-1,2) %*% t(X) + z1*((X[,1]+X[,2])^2) + rnorm(n=n,0,1)) # outcome
X.out1 <- cbind(X, (X[,1]+X[,2])^2)

### Simulation 2: small p
alpha2 <- c(0, 0.2, 0.8) # coefficients
z2 <- rbinom(n,1,prob=sigmoid(alpha2 %*% t(XD))) # treatment assignments
p2 <- mean(z2) # p
y2 <- as.numeric(c(-1,2) %*% t(X) + z2*((X[,1]+X[,2])^2) + rnorm(n=n,0,1)) # outcome
X.out2 <- cbind(X, (X[,1]+X[,2])^2)

### Propensity score plots
df <- data.frame(z1=z1, z2=z2, X1=X[,1], X2=X[,2], y1=y2, y2=y2) # data frame

bal.mult1 <- SumStat(ps.formula=z1 ~ X1+X2, data=df, weight="overlap")
# plot(bal.mult1, type="hist", breaks=30) 
bal.mult2 <- SumStat(ps.formula=z2 ~ X1+X2, data=df, weight="overlap")
# plot(bal.mult2, type="hist", breaks=30) 

### Variance ratio
fit <- glm(z1 ~ X1+X2, family = binomial(link = "logit"), data = df)
e.h <- as.numeric(fit$fitted.values)
r1 <- var(e.h[df$z1==1])/var(e.h[df$z1==0]) # ratio

fit <- glm(z2 ~ X1+X2, family = binomial(link = "logit"), data = df)
e.h <- as.numeric(fit$fitted.values)
r2 <- var(e.h[df$z2==1])/var(e.h[df$z2==0]) # ratio

### Causal effects
source("newSand_func.R")

ATE(y1, z1, X)  -> ATE.md1
ATO(y1, z1, X)  -> ATO.md1
ATM(y1, z1, X)  -> ATM.md1
ATEN(y1, z1, X) -> ATEN.md1
ATT(y1, z1, X)  -> ATT.md1
ATC(y1, z1, X)  -> ATC.md1

ATE(y1, z1, X, DR=TRUE, X.out=X.out1)  -> ATE.md1.aug
ATO(y1, z1, X, DR=TRUE, X.out=X.out1)  -> ATO.md1.aug
ATM(y1, z1, X, DR=TRUE, X.out=X.out1)  -> ATM.md1.aug
ATEN(y1, z1, X, DR=TRUE, X.out=X.out1) -> ATEN.md1.aug
ATT(y1, z1, X, DR=TRUE, X.out=X.out1)  -> ATT.md1.aug
ATC(y1, z1, X, DR=TRUE, X.out=X.out1)  -> ATC.md1.aug

ATE(y2, z2, X)  -> ATE.md2
ATO(y2, z2, X)  -> ATO.md2
ATM(y2, z2, X)  -> ATM.md2
ATEN(y2, z2, X) -> ATEN.md2
ATT(y2, z2, X)  -> ATT.md2
ATC(y2, z2, X)  -> ATC.md2

ATE(y2, z2, X, DR=TRUE, X.out=X.out2)  -> ATE.md2.aug
ATO(y2, z2, X, DR=TRUE, X.out=X.out2)  -> ATO.md2.aug
ATM(y2, z2, X, DR=TRUE, X.out=X.out2)  -> ATM.md2.aug
ATEN(y2, z2, X, DR=TRUE, X.out=X.out2) -> ATEN.md2.aug
ATT(y2, z2, X, DR=TRUE, X.out=X.out2)  -> ATT.md2.aug
ATC(y2, z2, X, DR=TRUE, X.out=X.out2)  -> ATC.md2.aug

### Summarize results
rows <- c("ATE", "ATO", "ATM", "ATEN", "ATT", "ATC")
est1.wt  <- c(ATE.md1$tau, ATO.md1$tau, ATM.md1$tau, ATEN.md1$tau, ATT.md1$tau, ATC.md1$tau)
est1.aug <- c(ATE.md1.aug$tau, ATO.md1.aug$tau, ATM.md1.aug$tau, ATEN.md1.aug$tau, ATT.md1.aug$tau, ATC.md1.aug$tau)
var1.wt  <- c(ATE.md1$se, ATO.md1$se, ATM.md1$se, ATEN.md1$se, ATT.md1$se, ATC.md1$se)
var1.aug <- c(ATE.md1.aug$se, ATO.md1.aug$se, ATM.md1.aug$se, ATEN.md1.aug$se, ATT.md1.aug$se, ATC.md1.aug$se)

est2.wt  <- c(ATE.md2$tau, ATO.md2$tau, ATM.md2$tau, ATEN.md2$tau, ATT.md2$tau, ATC.md2$tau)
est2.aug <- c(ATE.md2.aug$tau, ATO.md2.aug$tau, ATM.md2.aug$tau, ATEN.md2.aug$tau, ATT.md2.aug$tau, ATC.md2.aug$tau)
var2.wt  <- c(ATE.md2$se, ATO.md2$se, ATM.md2$se, ATEN.md2$se, ATT.md2$se, ATC.md2$se)
var2.aug <- c(ATE.md2.aug$se, ATO.md2.aug$se, ATM.md2.aug$se, ATEN.md2.aug$se, ATT.md2.aug$se, ATC.md2.aug$se)

p.value  <- function(x) ifelse(x<0.001, "<0.001", as.character(x))
SumTab1  <- data.frame(Estimand = rows, Est.wt = est1.wt, SE.wt = var1.wt, p.wt = p.value(1-2*pnorm(est1.wt/var1.wt)),
                       Est.aug = est1.aug, SE.aug = var1.aug, p.aug = p.value(1-2*pnorm(est1.aug/var1.aug)))
SumTab2  <- data.frame(Estimand = rows, Est.wt = est2.wt, SE.wt = var2.wt, p.wt = p.value(1-2*pnorm(est2.wt/var2.wt)),
                       Est.aug = est2.aug, SE.aug = var2.aug, p.aug = p.value(1-2*pnorm(est2.aug/var2.aug)))

print(xtable(SumTab1, digits = 2), include.rownames=FALSE)
print(xtable(SumTab2, digits = 2), include.rownames=FALSE)

png("toy_ps1.png", res=72*2, width = 1200, height = 700)
plot(bal.mult1, type="hist", breaks=30)
dev.off()

png("toy_ps2.png", res=72*2, width = 1200, height = 700)
plot(bal.mult2, type="hist", breaks=30)
dev.off()

data.frame(p1=p1, p2=p2, r1=r1, r2=r2)
