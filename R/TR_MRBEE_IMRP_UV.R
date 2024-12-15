#' Univariable MR Estimation using MRBEE Method
#'
#' This function performs univariable Mendelian Randomization (MR) analysis using
#' the MRBEE method.
#'
#' @param by Vector of GWAS effect sizes for the outcome (n x 1).
#' @param bx Vector of GWAS effect sizes for the exposure (n x 1).
#' @param byse Vector of standard errors (SE) for the GWAS effect sizes of the outcome (n x 1).
#' @param bxse Vector of SE for the GWAS effect sizes of the exposure (n x 1).
#' @param theta.source A scale of the causal effect estimate learning from the source data.
#' @param Rxy Correlation matrix (p+1 x p+1) of the exposures and outcome, with the outcome being the last.
#' @param susie.iter A scale of the maximum number of iterations used in SuSiE.
#' @param pip.thres A scale of PIP theshold for calibrating causality used in SuSiE.
#' @param max.iter Maximum number of iterations for the estimation process. Defaults to 30.
#' @param max.eps Tolerance level for convergence. Defaults to 1e-4.
#' @param pv.thres P-value threshold for pleiotropy detection. Defaults to 0.05.
#' @param var.est Method for estimating the standard error in the pleiotropy test. Can be "robust", "variance", or "ordinal".
#' @param FDR Logical indicating whether to apply False Discovery Rate (FDR) correction. Defaults to TRUE.
#' @param adjust.method Method for estimating q-values, defaults to "Sidak".
#' @param reliability.thres A scale of threshold for the minimum value of the reliability ratio. If the original reliability ratio is less than this threshold, only part of the estimation error is removed so that the working reliability ratio equals this threshold.
#' @param sampling.time A scale of number of replications in bootstrapping procedure.
#' @param sampling.iter A scale of iterations per bootstrapping procedure. Default is \code{10}.
#'
#' @return A list containing the estimated causal effect, its covariance, and pleiotropy.
#' @importFrom MASS rlm
#' @importFrom susieR susie coef.susie
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export


TR_MRBEE_IMRP_UV=function(by,bx,byse,bxse,Rxy,theta.source,susie.iter=100,pip.thres=0.3,max.iter=30,max.eps=1e-4,pv.thres=0.05,var.est="robust",FDR=T,adjust.method="Sidak",reliability.thres=0.8,sampling.time=100,sampling.iter=8){
by=by/byse
byseinv=1/byse
bx=bx*byseinv
bxse=bxse*byseinv
byse1=byse
byse=byse/byse
n=length(by)
r=reliability.adj.uv(bx,bxse,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bxse,Rxy)
########## Initial Estimation ############
eta0=bx%*%theta.source
br=by-eta0
colnames(bZ)[1]="Source Effect"
fit=rlm(br~bx-1)
theta.ini=fit$coefficient
theta=theta.ini
theta1=10000
e=c(br-bx*theta)
indvalid=which(abs(e)<=3*stats::mad(e))
indvalid=validadj(abs(e),indvalid,0.5)
## making the fraction of valid IVs must be larger than 50%
########## Iteration ###################
error=sqrt(sum((theta-theta1)^2))
iter=0
fit.susie=NULL
while(error>max.eps&iter<max.iter){
theta1=theta
e=c(br-bx*theta)
pv=imrpdetect(x=e,theta=theta[-1],RxyList=RxyList,var.est=var.est,FDR=FDR,adjust.method=adjust.method,indvalid=indvalid)
indvalid=which(pv>pv.thres)
if (length(indvalid) < length(pv) * 0.5) {
indvalid.cut = which(pv > stats::quantile(pv, 0.5))
indvalid = union(indvalid, indvalid.cut)
}
fit.susie=susie(y=br[indvalid],X=as.matrix(bx[indvalid]),L=1,estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F)
theta=coef.susie(fit.susie)[-1]*(fit.susie$pip>=pip.thres)
indtheta=which(theta!=0)
if(length(indtheta)==1){
h=sum(bx[indvalid]^2)-sum(bxse[indvalid]^2*Rxy[1,1])
g=sum(bx[indvalid]*br[indvalid])-Rxy[1,2]*sum(bxse[indvalid]*byse[indvalid])
theta=g/h
}
iter=iter+1
if(iter>5) error=sqrt(sum((theta-theta1)^2))
}
############################### inference #########################
r=c((br-bx*theta))*byse1
r[indvalid]=0
names(r)=names(bx)

ThetaList=c(1:sampling.time)*0
cat("Bootstrapping process:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
for(j in 1:sampling.time) {
setTxtProgressBar(pb, j)
indj <- sample(1:n, n, replace = T)
bxj=bx[indj]
brj=br[indj]
bxsej=bxse[indj]
bysej=byse[indj]
indvalidj=indvalid[indj]
thetaj=theta
theta1j=thetaj*0
fit.susiej=NULL
RxyListj=IVweight(byse[indj],bxse[indj],Rxy)
for(jiter in 1:sampling.iter){
ej=c(brj-bxj*thetaj)
pvj=imrpdetect(x=ej,theta=thetaj,RxyList=RxyListj,var.est=var.est,FDR=FDR,adjust.method=adjust.method,indvalid=indvalidj)
indvalidj=which(pvj>pv.thres)
if (length(indvalidj) < length(pvj) * 0.5) {
indvalid.cut = which(pv > stats::quantile(pvj, 0.5))
indvalidj = union(indvalidj, indvalid.cut)
}
fit.susiej=susie(y=brj[indvalidj],X=as.matrix(bxj[indvalidj]),L=1,estimate_prior_method="EM",residual_variance=1,s_init=fit.susiej,standardize=F,max_iter=susie.iter,intercept=F)
thetaj=coef.susie(fit.susiej)[-1]*(fit.susiej$pip>=pip.thres)
indthetaj=which(thetaj!=0)
if(length(indtheta)==1){
h=sum(bxj[indvalidj]^2)-sum(bxsej[indvalidj]^2*Rxy[1,1])
g=sum(bxj[indvalidj]*brj[indvalidj])-Rxy[1,2]*sum(bxsej[indvalidj]*bysej[indvalidj])
thetaj=g/h
}
}
ThetaList[j, ] <- thetaj
}
close(pb)
theta.se=sd(ThetaList)
theta.var=var(ThetaList)
A=list()
A$theta=theta
A$gamma=r
A$theta.se=theta.se
A$theta.z=theta/theta.se
A$theta.pip=mean(ThetaList!=0)
A$reliability.adjust=r
A$susie.theta=fit.susie
A$thetalist=ThetaList
return(A)
}
