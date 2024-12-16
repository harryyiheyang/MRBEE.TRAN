library(data.table)
library(dplyr)
library(MRBEEX)
load_all()
options(bitmapType="cairo")
EUR_Clumped <- readRDS("data/EUR_Clumped.rds") %>% list2env(.,envir=.GlobalEnv)
BETA=Z/sqrt(N)
SE=1/sqrt(N)
bX=BETA[,-10]
bXse=SE[,-10]
by=BETA[,10]
byse=SE[,10]
fit_EUR=MRBEE_IMRP(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,maxdiff=2)

EAS_Transfer <- readRDS("data/EAS_Transfer.rds") %>% list2env(.,envir=.GlobalEnv)
BETA=Z/sqrt(N)
SE=1/sqrt(N)
bX=BETA[,-c(1,11)]
bXse=SE[,-c(1,11)]
by=BETA[,11]
byse=SE[,11]
bz=BETA[,1]
bzse=SE[,1]
fit_EAS=MRBEE_IMRP(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,maxdiff=2)
fit_EAS$theta/fit_EAS$theta.se

cbind(fit_EUR$theta,fit_EAS$theta)
plot(fit_EUR$theta,fit_EAS$theta)
c(cor(fit_EUR$theta,fit_EAS$theta),cor(fit_EUR$theta,fit_EAS$theta,method="kendall"),cor(fit_EUR$theta,fit_EAS$theta,method="spearman"))

indEUR=which(abs(fit_EUR$theta/fit_EUR$theta.se)>2)
theta.source=fit_EUR$theta
theta.source[-indEUR]=0
fit_EAS_Tran=TR_MRBEE_IMRP(by,bX,byse,bXse,bz,bzse,Rxyz,L=6,theta.source=theta.source,transfer.coef=1)
cbind(sqrt(fit_EAS_Tran$delta.se[-1]^2+fit_EUR$theta.se^2),fit_EAS$theta.se)
Z=cbind(theta.source,fit_EAS_Tran$delta[-1]+theta.source,fit_EAS_Tran$delta[-1],fit_EAS$theta)
colnames(Z)=c("EUR","EAS_Tran","EAS_Diff","EAS")
print(Z)
fit_EAS_Tran$delta[1]
