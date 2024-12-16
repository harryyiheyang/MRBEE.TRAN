#' Estimate Causal Effect with MRBEE
#'
#' This function estimates the causal effect using a bias-correction estimating equation, considering potential pleiotropy and measurement errors.
#'
#' @param by A vector (n x 1) of the GWAS effect size of outcome.
#' @param bX A matrix (n x p) of the GWAS effect sizes of p exposures.
#' @param byse A vector (n x 1) of the GWAS effect size SE of outcome.
#' @param bXse A matrix (n x p) of the GWAS effect size SEs of p exposures.
#' @param theta.source A vector (p x 1) of the causal effect estimate learning from the source data.
#' @param Rxy A matrix (p+1 x p+1) of the correlation matrix of the p exposures and outcome. The last one should be the outcome.
#' @param L A scale of the number of single effects used in SuSiE.
#' @param susie.iter A scale of the maximum number of iterations used in SuSiE.
#' @param pip.thres A scale of PIP theshold for calibrating causality used in SuSiE.
#' @param max.iter Maximum number of iterations for causal effect estimation. Defaults to 30.
#' @param max.eps Tolerance for stopping criteria. Defaults to 1e-4.
#' @param pv.thres P-value threshold in pleiotropy detection. Defaults to 0.05.
#' @param var.est Method for estimating the variance of residual in pleiotropy test. Can be "robust", "variance", or "ordinal". Defaults is robust that estimates the variance of residual using median absolute deviation (MAD).
#' @param FDR Logical. Whether to apply the FDR to convert the p-value to q-value. Defaults to TRUE.
#' @param adjust.method Method for estimating q-value. Defaults to "Sidak".
#' @param reliability.thres A scale of threshold for the minimum value of the reliability ratio. If the original reliability ratio is less than this threshold, only part of the estimation error is removed so that the working reliability ratio equals this threshold.
#' @param sampling.time A scale of number of replications in bootstrapping procedure.
#' @param sampling.iter A scale of iterations per bootstrapping procedure. Default is \code{10}.
#' @param ridge.diff A scale of parameter on the differences of causal effect estimate in one credible set. Defaults to \code{10}.
#' @return A list containing the estimated causal effect, its covariance, and pleiotropy.
#' @importFrom susieR susie_suff_stat coef.susie susie
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'
TR_MRBEE_IMRP=function(by,bX,byse,bXse,theta.source,Rxy,L=min(10,ncol(bX)),susie.iter=500,pip.thres=10,max.iter=30,max.eps=1e-4,pv.thres=0.05,var.est="robust",FDR=T,adjust.method="Sidak",reliability.thres=0.8,sampling.time=100,sampling.iter=8,ridge.diff=10){
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
eta0=bX%*%theta.source
bZ=cbind(eta0,bX)
colnames(bZ)[1]="Source Effect"
br=by-eta0
fit=susie(y=br,X=bX,L=5)
theta.ini=c(1,coef.susie(fit)[-1])
theta=theta.ini
theta1=10000
e=c(by-bZ%*%theta)
indvalid=which(abs(e)<=3*stats::mad(e))
indvalid=validadj(abs(e),indvalid,0.5) ## making the fraction of valid IVs must be larger than 50%
########## Iteration ###################
error=sqrt(sum((theta-theta1)^2))
iter=0
fit.susie=NULL
while(error>max.eps&iter<max.iter){
theta1=theta
e=c(by-bZ%*%theta)
pv=imrpdetect(x=e,theta=theta[-1],RxyList=RxyList,var.est=var.est,FDR=FDR,adjust.method=adjust.method,indvalid=indvalid)
indvalid=which(pv>pv.thres)
if (length(indvalid) < length(pv) * 0.5) {
indvalid.cut = which(pv > stats::quantile(pv, 0.5))
indvalid = union(indvalid, indvalid.cut)
}
if(length(indvalid)==n){
Rxysum=Rxyall
Rxysum=cbind(0,rbind(0,Rxysum))
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:n,indvalid))
Rxysum=cbind(0,rbind(0,Rxysum))
}
ZtZ=t(bZ[indvalid,])%*%bZ[indvalid,]
Zty=c(t(bZ[indvalid,])%*%by[indvalid])
yty=sum(by[indvalid]^2)
fit.susie=susie_suff_stat(XtX=ZtZ,Xty=Zty,yty=yty,L=L,n=length(indvalid),estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F)
theta=coef.susie(fit.susie)[-1]*(fit.susie$pip>=pip.thres)
indtheta=which(theta!=0)
Diff=generate_block_matrix(summary(fit.susie)$vars,length(indvalid)/diag(ZtZ),theta)
if(length(indtheta)==1){
xtx=ZtZ[indtheta,indtheta]-Rxysum[indtheta,indtheta]
xty=Zty[indtheta]-Rxysum[indtheta,p+2]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
XtX=ZtZ[indtheta,indtheta]-Rxysum[indtheta,indtheta]+ridge.diff*Diff[indtheta,indtheta]
Xty=Zty[indtheta]-Rxysum[indtheta,p+2]
theta[indtheta]=c(solve(XtX)%*%Xty)
}
iter=iter+1
if(iter>5) error=sqrt(sum((theta-theta1)^2))
}
############################### inference #########################
theta=theta
names(theta)=colnames(bZ)
indtheta=which(theta!=0)
r=c(by-bZ%*%theta)*byse1
r[indvalid]=0
names(r)=rownames(bX)

ThetaList=matrix(0,sampling.time,p+1)
colnames(ThetaList)=colnames(bZ)
cat("Bootstrapping process:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
for(j in 1:sampling.time) {
setTxtProgressBar(pb, j)
indj <- sample(1:n, n, replace = T)
bZj=bZ[indj,]
byj=by[indj]
indvalidj=indvalid[indj]
thetaj=theta
theta1j=thetaj*0
fit.susiej=NULL
RxyListj=IVweight(byse[indj],bXse[indj],Rxy)
Rxyallj=biasterm(RxyList=RxyListj,c(1:n))
for(jiter in 1:sampling.iter){
ej=c(byj-bZj%*%thetaj)
pvj=imrpdetect(x=ej,theta=thetaj[-1],RxyList=RxyListj,var.est=var.est,FDR=FDR,adjust.method=adjust.method,indvalid=indvalidj)
indvalidj=which(pvj>pv.thres)
if (length(indvalidj) < length(pvj) * 0.5) {
indvalid.cut = which(pv > stats::quantile(pvj, 0.5))
indvalidj = union(indvalidj, indvalid.cut)
}
if(length(indvalidj)==n){
Rxysumj=Rxyallj
Rxysumj=cbind(0,rbind(0,Rxysumj))
}else{
Rxysumj=Rxyallj-biasterm(RxyList=RxyListj,setdiff(1:n,indvalidj))
Rxysumj=cbind(0,rbind(0,Rxysumj))
}
ZtZj=t(bZj[indvalidj,])%*%bZj[indvalidj,]
Ztyj=c(t(bZj[indvalidj,])%*%byj[indvalidj])
ytyj=sum(byj[indvalidj]^2)
fit.susiej=susie_suff_stat(XtX=ZtZj,Xty=Ztyj,yty=ytyj,L=L,n=n,estimate_prior_method="EM",residual_variance=1,s_init=fit.susiej,standardize=F,max_iter=susie.iter,intercept=F)
thetaj=coef.susie(fit.susiej)[-1]*(fit.susiej$pip>=pip.thres)
indthetaj=which(thetaj!=0)
Diffj=generate_block_matrix(summary(fit.susiej)$vars,length(indvalidj)/diag(ZtZj),thetaj)
if(length(indthetaj)==1){
xtx=ZtZj[indthetaj,indthetaj]-Rxysumj[indthetaj,indthetaj]
xty=Ztyj[indthetaj]-Rxysumj[indthetaj,p+2]
thetaj[indthetaj]=xty/xtx
}
if(length(indthetaj)>1){
XtX=ZtZj[indthetaj,indthetaj]-Rxysumj[indthetaj,indthetaj]+ridge.diff*Diffj
Xty=Ztyj[indthetaj]-Rxysumj[indthetaj,p+2]
thetaj[indthetaj]=c(solve(XtX)%*%Xty)
}

}
ThetaList[j, ] <- thetaj
}
close(pb)
theta.se=colSD(ThetaList)
theta.cov=cov(ThetaList)
colnames(theta.cov)=rownames(theta.cov)=names(theta.se)=colnames(bZ)

A=list()
A$theta=theta
A$gamma=r
A$theta.se=theta.se
A$theta.cov=theta.cov
A$theta.pip=colMeans(ThetaList!=0)
A$reliability.adjust=r
A$susie.theta=fit.susie
A$thetalist=ThetaList
return(A)

}
