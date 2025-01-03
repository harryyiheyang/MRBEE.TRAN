#' Estimate Causal Effect with MRcML and SuSiE
#'
#' This function estimates the causal effect using a bias-correction estimating equation, considering potential pleiotropy and measurement errors, and using SuSiE to select the causal exposures.
#'
#' @param by A vector (n x 1) of the GWAS effect size of outcome.
#' @param bX A matrix (n x p) of the GWAS effect sizes of p exposures.
#' @param byse A vector (n x 1) of the GWAS effect size SE of outcome.
#' @param bXse A matrix (n x p) of the GWAS effect size SEs of p exposures.
#' @param bX A vector (n x 1) of the liner predictor of the transferred exposure's effect.
#' @param bXse A vector (n x 1) of the standard error of the liner predictor of the transferred exposure's effect.
#' @param Rxy A matrix (p+1 x p+1) of the correlation matrix of the p exposures and outcome. The first one should be the transferred linear predictor and last one should be the outcome.
#' @param L A scale of the number of single effects used in SuSiE. Defaults to \code{5}.
#' @param susie.iter A scale of the maximum number of iterations used in SuSiE. Defaults to \code{200}.
#' @param tau A scale of tuning parameter in MCP. Defaults to \code{9}.
#' @param admm.rho A scale of tuning parameter in ADMM. Defaults to \code{3}.
#' @param pip.thres A scale of PIP threshold for calibyating causality used in SuSiE.
#' @param max.iter Maximum number of iterations for causal effect estimation. Defaults to \code{10}.
#' @param max.eps Tolerance for stopping criteria. Defaults to \code{1e-4}.
#' @param ridge.diff A scale of parameter on the differences of causal effect estimate in one credible set. Defaults to \code{10}.
#' @param sampling.time A scale of number of subsampling in estimating the standard error. Defaults to \code{100}.
#' @param sampling.iter A scale of iteration in subsampling in estimating the standard error. Defaults to \code{10}.
#' @param theta.ini An initial estimate of causal effect.
#'
#'
#' @return A list containing the estimated causal effect, its covariance, and pleiotropy.
#' @importFrom susieR susie_suff_stat coef.susie susie
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply matrixListProduct matrixInverse
#' @importFrom MASS rlm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'
MRcML_IPOD_SuSiE=function(by,bX,byse,bXse,Rxy,L=min(10,ncol(bX)),tau=9,admm.rho=3,susie.iter=200,pip.thres=0.2,max.iter=10,max.eps=1e-4,ridge.diff=100,sampling.time=100,sampling.iter=10,theta.ini=F){
if(theta.ini[1]==F){
fit.ini=MRBEE_IPOD_SuSiE(by,bX,byse,bXse,Rxy)
theta.ini=fit.ini$theta
}
######### Basic Processing  ##############
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
n=length(by)
p=ncol(bX)
GList=IVweight_CML(byse,bXse,Rxy)
OmegaxyList=GList$OmegaxyList
RxyList=GList$RxyList
e=c(by-bX%*%theta.ini)
gamma.ini=e
gamma.ini[which(abs(gamma.ini)<tau)]=0
########## Initial Estimation ############
error=2
iter=0
fit.susie=NULL
theta=theta.ini
gamma=gamma.ini
gamma1=u=gamma*0
evec=c()
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
bXest=CML_Bupdate(by=by-gamma,bX=bX,theta=theta,OmegaxyList=OmegaxyList)
varomega=OmegaxyList[,p+1,p+1]
XtX=matrixMultiply(t(bXest),bXest*varomega)
XtX=t(XtX)/2+XtX/2
g2=theta*0
for(i in indvalid){
s1=bX[i,]-bXest[i,]
s2=OmegaxyList[i,1:p,p+1]
s3=outer(s1,bXest[i,])
g2=g2+matrixVectorMultiply(s3,s2)
}
Xty=matrixVectorMultiply(t(bXest),(by-gamma)*varomega)+g2
yty=sum((by-gamma)^2*varomega)
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=L,n=n,estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F)
theta=coef.susie(fit.susie)[-1]*(fit.susie$pip>=pip.thres)
indtheta=which(theta!=0)
Diff=generate_block_matrix(summary(fit.susie)$vars,n/diag(XtX),theta)
if(length(indtheta)==1){
theta=CML_thetaupdate_uv(by=by-gamma,bX=bX,bXest=bXest,OmegaxyList=OmegaxyList,indvalid=indvalid,indtheta=indtheta,glatent=g2)
}
if(length(indtheta)>1){
theta=CML_thetaupdate(by=by-gamma,bX=bX,bXest=bXest,OmegaxyList=OmegaxyList,indvalid=indvalid,indtheta=indtheta,Diff=Diff,ridge.diff=ridge.diff,glatent=g2)
}
gamma=(by-matrixVectorMultiply(bXest,theta)-u+admm.rho*gamma1)/(1+admm.rho)
gamma1=mcp(gamma+u/admm.rho,tau/admm.rho)
u=u+admm.rho*(gamma-gamma1)
iter=iter+1
if(iter>5){
error=sqrt(sum((theta-theta1)^2))
}
evec[iter]=error
}
############################### inference #########################
names(theta)=colnames(bX)
indtheta=which(theta!=0)
res=gamma1*byse1
names(res)=rownames(bX)

ThetaList=matrix(0,sampling.time,p)
colnames(ThetaList)=colnames(bX)
cat("Bootstrapping process:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
for(j in 1:sampling.time) {
setTxtProgressBar(pb, j)
indj=sample(n,0.5*n,replace=F)
nj=length(indj)
bXj=bX[indj,]
bXestj=bXest[indj,]
byj=by[indj]
bXsej=bXse[indj,]
bysej=byse[indj]
thetaj=theta*runif(length(theta),0.95,1.05)
OmegaxyListj=GList$OmegaxyList[indj,,]
RxyListj=GList$RxyList[indj,,]
gammaj=gamma[indj]
gamma1j=uj=gammaj*0
indvalidj=which(gamma1j==0)
fit.susiej=fit.susie
for(iterj in 1:sampling.iter){
bXestj=CML_Bupdate(by=byj-gammaj,bX=bXj,theta=thetaj,OmegaxyList=OmegaxyListj)
indvalidj=which(gamma1j==0)
varomegaj=OmegaxyListj[,p+1,p+1]
XtXj=matrixMultiply(t(bXestj),bXestj*varomegaj)
XtXj=t(XtXj)/2+XtXj/2
g2j=theta*0
for(i in indvalidj){
s1j=bXj[i,]-bXestj[i,]
s2j=OmegaxyListj[i,1:p,p+1]
s3j=outer(s1j,bXestj[i,])
g2j=g2j+matrixVectorMultiply(s3j,s2j)
}
Xtyj=matrixVectorMultiply(t(bXestj),(byj-gammaj)*varomegaj)+g2j
ytyj=sum((byj)^2*varomegaj)
fit.susie=susie_suff_stat(XtX=XtXj,Xty=Xtyj,yty=ytyj,L=L,n=nj,estimate_prior_method="EM",residual_variance=1,s_init=fit.susiej,standardize=F,max_iter=15,intercept=F)
thetaj=coef.susie(fit.susiej)[-1]*(fit.susiej$pip>=pip.thres)
indthetaj=which(thetaj!=0)
Diffj=generate_block_matrix(summary(fit.susiej)$vars,nj/diag(XtXj),thetaj)
if(length(indthetaj)==1){
thetaj=CML_thetaupdate_uv(by=byj-gammaj,bX=bXj,bXest=bXestj,OmegaxyList=OmegaxyListj,indvalid=indvalidj,indtheta=indthetaj,glatent=g2j)
}
if(length(indthetaj)>1){
thetaj=CML_thetaupdate(by=byj-gammaj,bX=bXj,bXest=bXestj,OmegaxyList=OmegaxyListj,indvalid=indvalidj,indtheta=indthetaj,Diff=Diffj,ridge.diff=ridge.diff,glatent=g2j)
}
gammaj=as.vector(byj-matrixVectorMultiply(bXestj,thetaj)-uj+admm.rho*gamma1j)/(1+admm.rho)
gamma1j=mcp(gammaj+uj/admm.rho,tau/admm.rho)
uj=uj+admm.rho*(gammaj-gamma1j)
}
ThetaList[j,]=thetaj
}
close(pb)
theta.se=colSD(ThetaList)*sqrt(n/length(indvalid))
covtheta=cov(ThetaList)*n/length(indvalid)
colnames(covtheta)=rownames(covtheta)=names(theta.se)=colnames(bX)

A=list()
A$theta=theta
A$gamma=res
A$theta.se=sqrt(diag(covtheta))
A$theta.cov=covtheta
A$susie.theta=fit.susie
A$bXest=bXest
A$theta.list=ThetaList
A$algorithm.error=evec
return(A)

}
