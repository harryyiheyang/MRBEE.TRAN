library(data.table)
library(dplyr)
library(MRBEEX)
library(devtools)
load_all()
options(bitmapType="cairo")
EUR_Clumped <- readRDS("data/EUR_Clumped.rds") %>% list2env(.,envir=.GlobalEnv)
BETA=Z/sqrt(N)
SE=1/sqrt(N)
bX=BETA[,-10]
bXse=SE[,-10]
by=BETA[,10]
byse=SE[,10]
NAM=c(colnames(bX),"Type 2 Diabetes")
fit_EUR=MRBEE_IPOD_SuSiE(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy[NAM,NAM],Lvec=c(1:8))

EAS_Clumped <- readRDS("data/EAS_Clumped.rds") %>% list2env(.,envir=.GlobalEnv)
BETA=Z/sqrt(N)
SE=1/sqrt(N)
bX=BETA[,-10]
bXse=SE[,-10]
by=BETA[,10]
byse=SE[,10]
NAM=c(colnames(bX),"Type 2 Diabetes")
fit_EAS=MRBEE_IPOD_SuSiE(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy[NAM,NAM],Lvec=c(1:8))

cbind(fit_EUR$theta,fit_EAS$theta)
plot(fit_EUR$theta,fit_EAS$theta)
c(cor(fit_EUR$theta,fit_EAS$theta),cor(fit_EUR$theta,fit_EAS$theta,method="kendall"),cor(fit_EUR$theta,fit_EAS$theta,method="spearman"))

theta.source=fit_EUR$theta
fit_EAS_Tran=TR_MRBEE_IPOD(by,bX,byse,bXse,Rxy[NAM,NAM],Lvec=c(1:5),theta.source=fit_EUR$theta,theta.source.cov=fit_EUR$theta.cov,transfer.coef=1)
Estimate=cbind(theta.source,fit_EAS_Tran$theta,fit_EAS$theta)
SE=cbind(fit_EUR$theta.se,fit_EAS_Tran$theta.se,fit_EAS$theta.se)
Z=Estimate/SE
colnames(Z)=colnames(Estimate)=c("EUR","EAS_Tran","EAS")
print(Z)
print(Estimate)
fit_EAS_Tran$delta
par(mfrow=c(2,1))
barplot(fit_EAS$gamma)
barplot(fit_EAS_Tran$gamma)
par(mfrow=c(1,1))
plot(fit_EAS_Tran$Bic)

