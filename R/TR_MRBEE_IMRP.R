#' Estimate Causal Effect with MRBEE
#'
#' This function estimates the causal effect using a bias-correction estimating equation, considering potential pleiotropy and measurement errors.
#'
#' @param by A vector (n x 1) of the GWAS effect size of outcome.
#' @param bX A matrix (n x p) of the GWAS effect sizes of p exposures.
#' @param byse A vector (n x 1) of the GWAS effect size SE of outcome.
#' @param bXse A matrix (n x p) of the GWAS effect size SEs of p exposures.
#' @param bz A vector (n x 1) of the liner predictor of the transferred exposure's effect.
#' @param bzse A vector (n x 1) of the standard error of the liner predictor of the transferred exposure's effect.
#' @param transfer.coef A scale of transfer.coef of theta.source to theta.target. Defaults to 1. If \code{transfer.coef="adaptive"}, then the median regression coefficient between theta.source and theta.target.naive is used.
#' @param theta.source A vector (p x 1) of the causal effect estimate learning from the source data.
#' @param Rxyz A matrix (p+2 x p+2) of the correlation matrix of the p exposures and outcome. The first one should be the transferred linear predictor and last one should be the outcome.
#' @param L A scale of the number of single effects used in SuSiE.
#' @param susie.iter A scale of the maximum number of iterations used in SuSiE.
#' @param pip.thres A scale of PIP theshold for calibrating causality used in SuSiE.
#' @param max.iter Maximum number of iterations for causal effect estimation. Defaults to 100.
#' @param max.eps Tolerance for stopping criteria. Defaults to 1e-4.
#' @param pv.thres P-value threshold in pleiotropy detection. Defaults to 0.05.
#' @param var.est Method for estimating the variance of residual in pleiotropy test. Can be "robust", "variance", or "ordinal". Defaults is robust that estimates the variance of residual using median absolute deviation (MAD).
#' @param FDR Logical. Whether to apply the FDR to convert the p-value to q-value. Defaults to TRUE.
#' @param adjust.method Method for estimating q-value. Defaults to "Sidak".
#' @param reliability.thres A scale of threshold for the minimum value of the reliability ratio. If the original reliability ratio is less than this threshold, only part of the estimation error is removed so that the working reliability ratio equals this threshold.
#' @param ridge.diff A scale of parameter on the differences of causal effect estimate in one credible set. Defaults to \code{10}.
#' @return A list containing the estimated causal effect, its covariance, and pleiotropy.
#' @importFrom susieR susie_suff_stat coef.susie susie
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom MASS rlm
#' @importFrom MRBEEX MRBEE_IMRP
#' @export
#'
TR_MRBEE_IMRP=function(by,bX,byse,bXse,bz,bzse,Rxyz,L=min(10,ncol(bX)),transfer.coef=1,theta.source,susie.iter=500,pip.thres=0.2,max.iter=100,max.eps=1e-4,pv.thres=0.05,var.est="robust",FDR=T,adjust.method="Sidak",reliability.thres=0.8,ridge.diff=100){
######### Basic Processing  ##############
Rxy=Rxyz[-1,-1]
fit.no.tran=MRBEE_IMRP(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,maxdiff=2)
theta.naive=fit.no.tran$theta
if(transfer.coef[1]=="adaptive"){
fit.naive=rlm(theta.naive~theta.source-1)
delta.naive=coef(fit.naive)
}else{
delta.naive=1
}
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
bz=bz/byse
bzse=bzse/byse
byse=byse/byse
n=length(by)
p=ncol(bX)
bZ=cbind(bz,bX)
bZse=cbind(bzse,bXse)
colnames(bZ)[1]=colnames(bZse)[1]="Transfer Coefficient"
r=reliability.adj(bZ,bZse,thres=reliability.thres)
r=c(r,1)
Rxyz=t(t(Rxyz)*r)*r
RxyzList=IVweight(byse,bZse,Rxyz)
Rxyzall=biasterm(RxyList=RxyzList,c(1:n))
########## Initial Estimation ############
br=by-bz*delta.naive
fit=susie(y=br,X=bZ,L=5)
delta.ini=coef.susie(fit)[-1]*(fit$pip>0.3)
delta=delta.ini
delta1=10000
e=c(by-bZ%*%delta)
indvalid=which(abs(e)<=3*stats::mad(e))
indvalid=validadj(abs(e),indvalid,0.5) ## making the fraction of valid IVs must be larger than 50%
########## Iteration ###################
error=2
iter=0
fit.susie=NULL
delta=delta.ini
while(error>max.eps&iter<max.iter){
delta1=delta
e=c(by-bZ%*%delta)
pv=imrpdetect(x=e,theta=delta,RxyList=RxyzList,var.est=var.est,FDR=FDR,adjust.method=adjust.method,indvalid=indvalid)
indvalid=which(pv>pv.thres)
if (length(indvalid) < length(pv) * 0.5) {
indvalid.cut = which(pv > stats::quantile(pv, 0.5))
indvalid = union(indvalid, indvalid.cut)
}
if(length(indvalid)==n){
Rxyzsum=Rxyzall
}else{
Rxyzsum=Rxyzall-biasterm(RxyList=RxyzList,setdiff(1:n,indvalid))
}
ZtZ=t(bZ[indvalid,])%*%bZ[indvalid,]
Zty=c(t(bZ[indvalid,])%*%br[indvalid])
yty=sum((br[indvalid])^2)
fit.susie=susie_suff_stat(XtX=ZtZ,Xty=Zty,yty=yty,L=L,n=length(indvalid),estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F)
delta=coef.susie(fit.susie)[-1]*(fit.susie$pip>=pip.thres)
inddelta=which(delta!=0)
Diff=generate_block_matrix(summary(fit.susie)$vars,length(indvalid)/diag(ZtZ),delta)
if(length(inddelta)==1){
xtx=ZtZ[inddelta,inddelta]-Rxyzsum[inddelta,inddelta]
xty=Zty[inddelta]-Rxyzsum[inddelta,p+2]
delta[inddelta]=xty/xtx
}
if(length(inddelta)>1){
xtx=ZtZ[inddelta,inddelta]-Rxyzsum[inddelta,inddelta]+ridge.diff*Diff[inddelta,inddelta]
xty=Zty[inddelta]-Rxyzsum[inddelta,p+2]
delta[inddelta]=c(solve(xtx)%*%xty)
}
iter=iter+1
if(iter>5){
error=sqrt(sum((delta-delta1)^2))
}
}
############################### inference #########################
names(delta)=colnames(bZ)
inddelta=which(delta!=0)
r=c(by-bZ%*%delta)*byse1
r[indvalid]=0
names(r)=rownames(bZ)

adjf=n/(length(indvalid)-length(inddelta)-1)
if(length(inddelta)>0){
Hinv=solve(xtx)
D=bZ[indvalid,inddelta]%*%(Hinv%*%t(bZ[indvalid,inddelta]))
D=(rep(1,length(indvalid))-diag(D))
D[which(D<0.25)]=0.25
E=-bZ[indvalid,inddelta]*(e[indvalid]/D)
for(i in 1:length(indvalid)){
E[i,]=E[i,]+RxyzList[indvalid[i],p+2,inddelta]-c(RxyzList[indvalid[i],inddelta,inddelta]%*%delta[inddelta])
}
V=t(E)%*%E
covdelta=diag(delta)*0
covdelta[inddelta,inddelta]=(Hinv%*%V%*%Hinv)*adjf
}else{
covdelta=diag(delta)*0
}
colnames(covdelta)=rownames(covdelta)=colnames(bZ)
A=list()
A$delta=delta
A$gamma=r
A$delta.se=sqrt(diag(covdelta))
A$delta.cov=covdelta
A$reliability.adjust=r
A$susie.delta=fit.susie
return(A)

}
