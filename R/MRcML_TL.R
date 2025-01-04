#' Estimate Non-Transferable Causal Effect with MRcML and SuSiE
#'
#' This function estimates the non-transferable causal effect using a conditional score method, considering potential pleiotropy and measurement errors, and using SuSiE to select the non-transferable causal effect.
#'
#' @param by A vector (n x 1) of the GWAS effect size of outcome.
#' @param bX A matrix (n x p) of the GWAS effect sizes of p exposures.
#' @param byse A vector (n x 1) of the GWAS effect size SE of outcome.
#' @param bXse A matrix (n x p) of the GWAS effect size SEs of p exposures.
#' @param bX A vector (n x 1) of the liner predictor of the transferred exposure's effect.
#' @param bXse A vector (n x 1) of the standard error of the liner predictor of the transferred exposure's effect.
#' @param Rxy A matrix (p+1 x p+1) of the correlation matrix of the p exposures and outcome. The first one should be the transferred linear predictor and last one should be the outcome.
#' @param theta.source A vector (p x 1) of the causal effect estimate learning from the source data.
#' @param theta.source.cov A matrix (p x p) of the covariance matrix of the causal effect estimate learning from the source data.
#' @param transfer.coef A scale of transfer.coef of theta.source to theta.target. Default is \code{1}.
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
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'
MRcML_TL=function(by,bX,byse,bXse,Rxy,L=min(10,ncol(bX)),tau=9,admm.rho=3,susie.iter=200,pip.thres=0.2,max.iter=10,max.eps=1e-4,ridge.diff=100,sampling.time=100,sampling.iter=10,theta.ini=F,theta.source,theta.source.cov,transfer.coef=1){
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
delta=theta.source-theta.ini
while(error>max.eps&iter<max.iter){
delta1=delta
indvalid=which(gamma1==0)
bXest=CML_Bupdate(by=by-gamma,bX=bX,theta=theta,OmegaxyList=OmegaxyList)
varomega=OmegaxyList[,p+1,p+1]
XtX=matrixMultiply(t(bXest),bXest*varomega)
XtX=t(XtX)/2+XtX/2
g2=delta*0
for(i in indvalid){
s1=bX[i,]-bXest[i,]
s2=OmegaxyList[i,1:p,p+1]
s3=outer(s1,bXest[i,])
g2=g2+matrixVectorMultiply(s3,s2)
}
fit.cluster=center.classifying(delta,-theta.source)
delta.complement=fit.cluster$complement
delta.cluster=fit.cluster$cluster
br.complement=c(by-bXest%*%(theta.source+delta.complement)-gamma)
Xty=matrixVectorMultiply(t(bXest),br.complement*varomega)+g2
yty=sum(br.complement^2*varomega)
tryCatch({
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=L,n=n,estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F,residual_variance_lowerbound=1)
},error = function(e) {
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=L,n=n,estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F,estimate_residual_variance=F)
})
delta.latent=coef.susie(fit.susie)[-1]*(fit.susie$pip>=pip.thres)
inddelta=which(delta.latent!=0)
Diff=generate_block_matrix(summary(fit.susie)$vars,n/diag(XtX),delta.latent)
if(length(inddelta)==1){
delta.latent=CML_thetaupdate_uv(by=br.complement,bX=bX,bXest=bXest,OmegaxyList=OmegaxyList,indvalid=indvalid,indtheta=inddelta,glatent=g2)
}
if(length(inddelta)>1){
delta.latent=CML_thetaupdate(by=br.complement,bX=bX,bXest=bXest,OmegaxyList=OmegaxyList,indvalid=indvalid,indtheta=inddelta,Diff=Diff,ridge.diff=ridge.diff,glatent=g2)
}
delta=delta.latent+delta.complement
theta=delta+theta.source
gamma=(by-matrixVectorMultiply(bXest,theta)-u+admm.rho*gamma1)/(1+admm.rho)
gamma1=mcp(gamma+u/admm.rho,tau/admm.rho)
u=u+admm.rho*(gamma-gamma1)
iter=iter+1
if(iter>5){
error=sqrt(sum((delta-delta1)^2))
}
evec[iter]=sqrt(sum((delta-delta1)^2))
}
############################### inference #########################
names(theta)=names(delta)=colnames(bX)
inddelta=which(delta!=0)
res=by-matrixVectorMultiply(bXest,theta)*byse1
res[indvalid]=0
names(res)=rownames(bX)

ThetaList=DeltaList=matrix(0,sampling.time,p)
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
deltaj=theta.source-thetaj
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
fit.clusterj=center.classifying(deltaj,-theta.source)
delta.complementj=fit.clusterj$complement
delta.clusterj=fit.clusterj$clusterj
br.complementj=c(byj-bXestj%*%(theta.source+delta.complementj)-gammaj)
Xtyj=matrixVectorMultiply(t(bXestj),br.complementj*varomegaj)+g2j
ytyj=sum(br.complementj^2*varomegaj)
fit.susie=susie_suff_stat(XtX=XtXj,Xty=Xtyj,yty=ytyj,L=L,n=nj,estimate_prior_method="EM",residual_variance=1,s_init=fit.susiej,standardize=F,max_iter=susie.iter,intercept=F)
delta.latentj=coef.susie(fit.susiej)[-1]*(fit.susiej$pip>=pip.thres)
inddeltaj=which(delta.latentj!=0)
Diffj=generate_block_matrix(summary(fit.susiej)$vars,nj/diag(XtXj),delta.latentj)
if(length(inddeltaj)==1){
delta.latentj=CML_thetaupdate_uv(by=br.complementj,bX=bXj,bXest=bXestj,OmegaxyList=OmegaxyListj,indvalid=indvalidj,indtheta=inddeltaj,glatent=g2j)
}
if(length(inddeltaj)>1){
delta.latentj=CML_thetaupdate(by=br.complementj,bX=bXj,bXest=bXestj,OmegaxyList=OmegaxyListj,indvalid=indvalidj,indtheta=inddeltaj,Diff=Diffj,ridge.diff=ridge.diff,glatent=g2j)
}
deltaj=delta.latentj+delta.complementj
thetaj=deltaj+theta.source
gammaj=as.vector(byj-matrixVectorMultiply(bXestj,thetaj)-uj+admm.rho*gamma1j)/(1+admm.rho)
gamma1j=mcp(gammaj+uj/admm.rho,tau/admm.rho)
uj=uj+admm.rho*(gammaj-gamma1j)
}
ThetaList[j,]=thetaj
DeltaList[j,]=deltaj
}
close(pb)
theta.cov=cov(ThetaList)*n/length(indvalid)
theta.cov[which(delta==0),which(delta==0)]=theta.cov[which(delta==0),which(delta==0)]+theta.source.cov[which(delta==0),which(delta==0)]
theta.se=sqrt(diag(theta.cov))
colnames(theta.cov)=rownames(theta.cov)=names(theta.se)=colnames(bX)
delta.se=colSD(DeltaList)*sqrt(n/length(indvalid))
delta.cov=cov(DeltaList)*n/length(indvalid)
colnames(delta.cov)=rownames(delta.cov)=names(delta.se)=colnames(bX)

A=list()
A$delta=delta
A$theta=theta.source+delta
A$gamma=gamma1*byse1
A$delta.se=delta.se
A$theta.se=theta.se
A$delta.cov=delta.cov
A$theta.cov=theta.cov
A$susie.delta=fit.susie
A$delta.latent=delta.latent
A$delta.list=DeltaList
A$theta.list=ThetaList
A$algorithm.error=evec
return(A)

}
