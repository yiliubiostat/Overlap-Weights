### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ What are we weighting for ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### Toy example for illustration the impact of p=P(Z=1)
### by Yi Liu

library(dplyr)
library(ggplot2)

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

### Simulation 2: small p
alpha2 <- c(0, 0.2, 0.8) # coefficients
z2 <- rbinom(n,1,prob=sigmoid(alpha2 %*% t(XD))) # treatment assignments
p2 <- mean(z2) # p
y2 <- as.numeric(c(-1,2) %*% t(X) + z2*((X[,1]+X[,2])^2) + rnorm(n=n,0,1)) # outcome

df <- data.frame(z1=z1, z2=z2, X1=X[,1], X2=X[,2], y1=y2, y2=y2) # data frame

### Variance ratio
XD <- cbind(rep(1,n), X)
ps1 <- 1/(1+exp(-alpha1 %*% t(XD))) %>% as.numeric() # true propensity score
ps2 <- 1/(1+exp(-alpha2 %*% t(XD))) %>% as.numeric()

r1 <- var(ps1[df$z1==1])/var(ps1[df$z1==0]) # ratio
r2 <- var(ps2[df$z2==1])/var(ps2[df$z2==0]) # ratio

### Causal effects
delta <- (df$X1+df$X2)^2 # mean(delta) is the estimated ATE
ATE <- 6^2+9 + 0.75^2+0.25*0.75 + 2*6*0.75 # true ATE mathematically
# mean(delta)-ATE # see bias
ATT1 <- sum(delta*ps1)/sum(ps1) # For other estimands, we just approximate them
ATT2 <- sum(delta*ps2)/sum(ps2)

ATC1 <- sum(delta*(1-ps1))/sum(1-ps1)
ATC2 <- sum(delta*(1-ps2))/sum(1-ps2)

ATO1 <- sum(delta*ps1*(1-ps1))/sum(ps1*(1-ps1))
ATO2 <- sum(delta*ps2*(1-ps2))/sum(ps2*(1-ps2))

u1 = apply(cbind(ps1, 1-ps1), 1, min)
ATM1 <- sum(delta*u1)/sum(u1)
u2 = apply(cbind(ps2, 1-ps2), 1, min)
ATM2 <- sum(delta*u2)/sum(u2)

ksi1 = -(ps1*log(ps1) + (1-ps1)*log(1-ps1))
ATEN1 <- sum(delta*ksi1)/sum(ksi1)
ksi2 = -(ps2*log(ps2) + (1-ps2)*log(1-ps2))
ATEN2 <- sum(delta*ksi2)/sum(ksi2)

rm(u1,u2,ksi1,ksi2,delta)

### Propensity score plots
### Propensity score plots
df <- df %>% mutate(ps1=ps1, ps2=ps2)
png("toy_ps1.png", res=72*2, width = 850, height = 600)
ggplot(df, aes(x=ps1, fill=as.factor(z1))) + 
  geom_histogram(position="identity", alpha=0.5, color="black", bins=40) +
  labs(x="Propensity score", y="") + 
  scale_fill_manual(values=c("darkblue", "gray"), name = "Group") + 
  theme_classic()
dev.off()

png("toy_ps2.png", res=72*2, width = 850, height = 600)
ggplot(df, aes(x=ps2, fill=as.factor(z2))) + 
  geom_histogram(position="identity", alpha=0.5, color="black", bins=40) +
  labs(x="Propensity score", y="") + 
  scale_fill_manual(values=c("darkblue", "gray"), name = "Group") + 
  theme_classic()
dev.off()

round(c(ATE, ATO1, ATM1, ATEN1, ATT1, ATC1), 2)
round(c(ATE, ATO2, ATM2, ATEN2, ATT2, ATC2), 2)
data.frame(p1=p1, p2=p2, r1=r1, r2=r2)
