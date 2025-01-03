library(data.table)
library(dplyr)
library(MRBEEX)
library(devtools)
load_all()
options(bitmapType="cairo")
EUR_Clumped <- readRDS("data/EUR_Clumped.rds") %>% list2env(.,envir=.GlobalEnv)
BETA=Z/sqrt(N)
SE=1/sqrt(N)
bx=BETA[,3]
bxse=SE[,3]
by=BETA[,10]
byse=SE[,10]
NAME=colnames(BETA)
NAM=NAME[c(1,10)]
fit_EUR=MRBEE_IMRP_UV(by=by,bx=bx,byse=byse,bxse=bxse,Rxy=Rxy[NAM,NAM])

EAS_Clumped <- readRDS("data/EAS_Clumped.rds") %>% list2env(.,envir=.GlobalEnv)
BETA=Z/sqrt(N)
SE=1/sqrt(N)
bx=BETA[,3]
bxse=SE[,3]
by=BETA[,10]
byse=SE[,10]
NAM=NAME[c(1,10)]
fit_EAS=MRBEE_IMRP_UV(by=by,bx=bx,byse=byse,bxse=bxse,Rxy=Rxy[NAM,NAM])
fit_EAS_Tran=TR_MRBEE_IMRP_UV(by=by,bx=bx,byse=byse,bxse=bxse,Rxy=Rxy[NAM,NAM],theta.source=fit_EUR$theta)
c(fit_EAS_Tran$theta,fit_EAS$theta,fit_EUR$theta)
par(mfrow=c(2,1))
barplot(fit_EAS$gamma)
barplot(fit_EAS_Tran$gamma)
