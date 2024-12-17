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
fit_EUR=MRBEE_IMRP_SuSiE(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,L=10)

EAS_Transfer <- readRDS("data/EAS_Transfer.rds") %>% list2env(.,envir=.GlobalEnv)
BETA=Z/sqrt(N)
SE=1/sqrt(N)
bX=BETA[,-c(1,11)]
bXse=SE[,-c(1,11)]
by=BETA[,11]
byse=SE[,11]
bz=BETA[,1]
bzse=SE[,1]
fit_EAS=MRBEE_IMRP_SuSiE(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,L=8)
fit_EAS$theta/fit_EAS$theta.se

cbind(fit_EUR$theta,fit_EAS$theta)
plot(fit_EUR$theta,fit_EAS$theta)
c(cor(fit_EUR$theta,fit_EAS$theta),cor(fit_EUR$theta,fit_EAS$theta,method="kendall"),cor(fit_EUR$theta,fit_EAS$theta,method="spearman"))

theta.source=fit_EUR$theta
susie.iter=500;pip.thres=0.2;max.iter=100;max.eps=1e-4;pv.thres=0.05;var.est="variance";FDR=T;adjust.method="Sidak";reliability.thres=0.8;ridge.diff=100
L=8;transfer.coef=1
Rxy=Rxyz[-1,-1]
fit_EAS_Tran=TR_MRBEE_IMRP(by,bX,byse,bXse,Rxy,L=8,theta.source=theta.source,transfer.coef=1)
cbind(sqrt(fit_EAS_Tran$delta.se^2+fit_EUR$theta.se^2),fit_EAS$theta.se)
Estimate=cbind(theta.source,fit_EAS_Tran$delta+theta.source*fit_EAS_Tran$estimate.transfer.coef,fit_EAS_Tran$delta,fit_EAS$theta)
SE=cbind(fit_EUR$theta.se,sqrt(fit_EAS_Tran$delta.se^2+fit_EUR$theta.se^2),fit_EAS_Tran$delta.se,fit_EAS$theta.se)
Z=Estimate/SE
colnames(Z)=colnames(Estimate)=c("EUR","EAS_Tran","EAS_Diff","EAS")
print(Z)
print(Estimate)
