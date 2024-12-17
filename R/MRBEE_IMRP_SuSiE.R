#' Estimate Causal Effect with MRBEE and SuSiE
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
#' @param L A scale of the number of single effects used in SuSiE.
#' @param susie.iter A scale of the maximum number of iterations used in SuSiE.
#' @param pip.thres A scale of PIP theshold for calibyating causality used in SuSiE.
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
MRBEE_IMRP_SuSiE=function(by,bX,byse,bXse,Rxy,L=min(10,ncol(bX)),susie.iter=500,pip.thres=0.2,max.iter=100,max.eps=1e-4,pv.thres=0.05,var.est="variance",FDR=T,adjust.method="Sidak",reliability.thres=0.8,ridge.diff=100){
######### Basic Processing  ##############
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
fit=susie(y=by,X=bX,L=5)
theta.ini=coef.susie(fit)[-1]*(fit$pip>0.3)
theta=theta.ini
theta1=10000
e=c(by-bX%*%theta)
indvalid=which(abs(e)<=3*stats::mad(e))
indvalid=validadj(abs(e),indvalid,0.5) ## making the fraction of valid IVs must be larger than 50%
########## Iteration ###################
error=2
iter=0
fit.susie=NULL
theta=theta.ini
while(error>max.eps&iter<max.iter){
theta1=theta
e=c(by-bX%*%theta)
pv=imrpdetect(x=e,theta=theta,RxyList=RxyList,var.est=var.est,FDR=FDR,adjust.method=adjust.method,indvalid=indvalid)
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
XtX=t(bX[indvalid,])%*%bX[indvalid,]
XtX=t(XtX)/2+XtX/2
Xty=c(t(bX[indvalid,])%*%by[indvalid])
yty=sum((by[indvalid])^2)
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=L,n=length(indvalid),estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F)
theta=coef.susie(fit.susie)[-1]*(fit.susie$pip>=pip.thres)
indtheta=which(theta!=0)
Diff=generate_block_matrix(summary(fit.susie)$vars,length(indvalid)/diag(XtX),theta)
if(length(indtheta)==1){
xtx=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
xtx=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]+ridge.diff*Diff[indtheta,indtheta]
xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=c(solve(xtx)%*%xty)
}
iter=iter+1
if(iter>5){
error=sqrt(sum((theta-theta1)^2))
}
}
############################### inference #########################
names(theta)=colnames(bX)
indtheta=which(theta!=0)
res=c(by-bX%*%theta)*byse1
res[indvalid]=0
names(res)=rownames(bX)

adjf=n/(length(indvalid)-length(indtheta)-1)
if(length(indtheta)>0){
Hinv=solve(xtx)
D=bX[indvalid,indtheta]%*%(Hinv%*%t(bX[indvalid,indtheta]))
D=(rep(1,length(indvalid))-diag(D))
D[which(D<0.25)]=0.25
E=-bX[indvalid,indtheta]*(e[indvalid]/D)
for(i in 1:length(indvalid)){
E[i,]=E[i,]+RxyList[indvalid[i],p+1,indtheta]-c(RxyList[indvalid[i],indtheta,indtheta]%*%theta[indtheta])
}
V=t(E)%*%E
covtheta=diag(theta)*0
covtheta[indtheta,indtheta]=(Hinv%*%V%*%Hinv)*adjf
}else{
covtheta=diag(theta)*0
}
colnames(covtheta)=rownames(covtheta)=colnames(bX)
A=list()
A$theta=theta
A$gamma=res
A$theta.se=sqrt(diag(covtheta))
A$theta.cov=covtheta
A$reliability.adjust=r
A$susie.theta=fit.susie
return(A)

}
