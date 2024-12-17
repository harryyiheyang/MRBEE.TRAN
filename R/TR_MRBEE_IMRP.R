#' Estimate Non-Transferable Causal Effect with MRBEE and SuSiE
#'
#' This function estimates the non-transferable causal effect using a bias-correction estimating equation, considering potential pleiotropy and measurement errors, and using SuSiE to select the non-transferable causal effect.
#'
#' @param by A vector (n x 1) of the GWAS effect size of outcome.
#' @param bX A matrix (n x p) of the GWAS effect sizes of p exposures.
#' @param byse A vector (n x 1) of the GWAS effect size SE of outcome.
#' @param bXse A matrix (n x p) of the GWAS effect size SEs of p exposures.
#' @param transfer.coef A scale of transfer.coef of theta.source to theta.target. Defaults to 1. If \code{transfer.coef="adaptive"}, then the median regression coefficient between theta.source and theta.target.naive is used.
#' @param theta.source A vector (p x 1) of the causal effect estimate learning from the source data.
#' @param Rxy A matrix (p+1 x p+1) of the correlation matrix of the p exposures and outcome. The first one should be the transferred linear predictor and last one should be the outcome.
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
TR_MRBEE_IMRP=function(by,bX,byse,bXse,Rxy,L=min(10,ncol(bX)),transfer.coef=1,theta.source,susie.iter=500,pip.thres=0.2,max.iter=100,max.eps=1e-4,pv.thres=0.05,var.est="variance",FDR=T,adjust.method="Sidak",reliability.thres=0.8,ridge.diff=100){
######### Basic Processing  ##############
fit.no.tran=MRBEE_IMRP(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy)
# we don't consider sparse estimate because it may cause singular design in rlm(.)
theta.naive=fit.no.tran$theta
if(transfer.coef[1]=="adaptive"){
fit.naive=rlm(theta.naive~theta.source-1)
delta.naive=coef(fit.naive)
pv=pchisq((delta.naive-1)^2/summary(fit.naive)[[4]][2]^2,1,lower.tail=F)
if(pv>0.05){
delta.naive=1
}
}else{
delta.naive=transfer.coef
}
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
n=length(by)
p=ncol(bX)
r=reliability.adj(bX,bXse,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy)
Rxyall=biasterm(RxyList=RxyList,c(1:n))
########## Initial Estimation ############
br=as.vector(by-bX%*%theta.source*delta.naive)
fit=susie(y=br,X=bX,L=5)
delta.ini=coef.susie(fit)[-1]*(fit$pip>0.3)
delta=delta.ini
delta1=10000
e=c(br-bX%*%delta)
indvalid=which(abs(e)<=3*stats::mad(e))
indvalid=validadj(abs(e),indvalid,0.5) ## making the fraction of valid IVs must be larger than 50%
########## Iteration ###################
error=2
iter=0
fit.susie=NULL
delta=delta.ini
while(error>max.eps&iter<max.iter){
delta1=delta
e=c(br-bX%*%delta)
pv=imrpdetect(x=e,theta=delta,RxyList=RxyList,var.est=var.est,FDR=FDR,adjust.method=adjust.method,indvalid=indvalid)
indvalid=which(pv>pv.thres)
if (length(indvalid) < length(pv) * 0.5) {
indvalid.cut = which(pv > stats::quantile(pv, 0.5))
indvalid = union(indvalid, indvalid.cut)
}
if(length(indvalid)==n){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:n,indvalid))
}
theta.complement=center.classifying(delta,theta.source*delta.naive)
br.complement=c(br-bX%*%theta.complement)
XtX=t(bX[indvalid,])%*%bX[indvalid,]
XtX=t(XtX)/2+XtX/2
Xty=c(t(bX[indvalid,])%*%br.complement[indvalid])
yty=sum((br.complement[indvalid])^2)
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=L,n=length(indvalid),estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F)
delta.latent=coef.susie(fit.susie)[-1]*(fit.susie$pip>=pip.thres)
inddelta=which(delta.latent!=0)
Diff=generate_block_matrix(summary(fit.susie)$vars,length(indvalid)/diag(XtX),delta.latent)
if(length(inddelta)==1){
xtx=XtX[inddelta,inddelta]-Rxysum[inddelta,inddelta]
xty=Xty[inddelta]-Rxysum[inddelta,p+1]
delta.latent[inddelta]=xty/xtx
}
if(length(inddelta)>1){
xtx=XtX[inddelta,inddelta]-Rxysum[inddelta,inddelta]+ridge.diff*Diff[inddelta,inddelta]
xty=Xty[inddelta]-Rxysum[inddelta,p+1]
delta.latent[inddelta]=c(solve(xtx)%*%xty)
}
delta=delta.latent+theta.complement
iter=iter+1
if(iter>5){
error=sqrt(sum((delta-delta1)^2))
}
}
############################### inference #########################
names(delta)=colnames(bX)
inddelta=which(delta.latent!=0)
res=c(br-bX%*%delta)*byse1
res[indvalid]=0
names(res)=rownames(bX)

adjf=n/(length(indvalid)-length(inddelta)-1)
if(length(inddelta)>0){
Hinv=solve(xtx)
D=bX[indvalid,inddelta]%*%(Hinv%*%t(bX[indvalid,inddelta]))
D=(rep(1,length(indvalid))-diag(D))
D[which(D<0.25)]=0.25
E=-bX[indvalid,inddelta]*(e[indvalid]/D)
for(i in 1:length(indvalid)){
E[i,]=E[i,]+RxyList[indvalid[i],p+1,inddelta]-c(RxyList[indvalid[i],inddelta,inddelta]%*%delta[inddelta])
}
V=t(E)%*%E
covdelta=diag(delta)*0
covdelta[inddelta,inddelta]=(Hinv%*%V%*%Hinv)*adjf
}else{
covdelta=diag(delta)*0
}
colnames(covdelta)=rownames(covdelta)=colnames(bX)
A=list()
A$delta=delta
A$gamma=res
A$delta.se=sqrt(diag(covdelta))
A$delta.cov=covdelta
A$reliability.adjust=r
A$susie.delta=fit.susie
A$estimate.transfer.coef=delta.naive
A$delta.latent=delta.latent
return(A)

}
