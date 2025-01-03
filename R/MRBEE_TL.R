#' Estimate Non-Transferable Causal Effect with MRBEE and SuSiE
#'
#' This function estimates the non-transferable causal effect using a bias-correction estimating equation, considering potential pleiotropy and measurement errors, and using SuSiE to select the non-transferable causal effect.
#'
#' @param by A vector (n x 1) of the GWAS effect size of outcome.
#' @param bX A matrix (n x p) of the GWAS effect sizes of p exposures.
#' @param byse A vector (n x 1) of the GWAS effect size SE of outcome.
#' @param bXse A matrix (n x p) of the GWAS effect size SEs of p exposures.
#' @param theta.source A vector (p x 1) of the causal effect estimate learning from the source data.
#' @param theta.source.cov A matrix (p x p) of the covariance matrix of the causal effect estimate learning from the source data.
#' @param Rxy A matrix (p+1 x p+1) of the correlation matrix of the p exposures and outcome. The first one should be the transferred linear predictor and last one should be the outcome.
#' @param transfer.coef A scale of transfer.coef of theta.source to theta.target. Default is \code{1}.
#' @param tauvec The candidate vector of tuning parameters for the MCP penalty function. Default is \code{seq(3, 30, by=3)}.
#' @param admm.rho The tuning parameter in the nested ADMM algorithm. Default is \code{3}.
#' @param Lvec A vector of the number of single effects used in SuSiE. Default is \code{c(1:6)}.
#' @param susie.iter A scale of the maximum number of iterations used in SuSiE. Default is \code{200}.
#' @param pip.thres A scale of PIP theshold for calibyating causality used in SuSiE. Default is \code{0.3}.
#' @param ebic.delta A scale of tuning parameter of causal effect estimate in extended BIC. Default is \code{1}.
#' @param ebic.gamma A scale of tuning parameter of horizontal pleiotropy in extended BIC. Default is \code{2}.
#' @param max.iter Maximum number of iterations for causal effect estimation. Default is \code{50}.
#' @param max.eps Tolerance for stopping criteria. Default is \code{1e-4}.
#' @param reliability.thres A scale of threshold for the minimum value of the reliability ratio. If the original reliability ratio is less than this threshold, only part of the estimation error is removed so that the working reliability ratio equals this threshold. Default is \code{0.8}.
#' @param ridge.diff A scale of parameter on the differences of causal effect estimate in one credible set. Defaults to \code{10}.
#' @param sampling.time A scale of number of subsampling in estimating the standard error. Default is \code{100}.
#' @param sampling.iter A scale of iteration in subsampling in estimating the standard error. Default is \code{10}.
#'
#' @return A list containing the estimated causal effect, its covariance, and pleiotropy.
#' @importFrom susieR susie_suff_stat coef.susie susie
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply matrixListProduct
#' @importFrom MASS rlm
#' @importFrom MRBEEX MRBEE_IMRP
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'
MRBEE_TL=function(by,bX,byse,bXse,Rxy,theta.source,theta.source.cov,tauvec=seq(3,30,3),admm.rho=3,Lvec=c(1:6),ebic.delta=1,ebic.gamma=2,transfer.coef=1,susie.iter=200,pip.thres=0.3,max.iter=50,max.eps=1e-4,reliability.thres=0.8,ridge.diff=100,sampling.time=100,sampling.iter=10){
######### Basic Processing  ##############
fit.no.tran=MRBEE_IMRP(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy)
theta.source=transfer.coef*theta.source
theta.source.cov=transfer.coef^2*theta.source.cov
theta.ini=fit.no.tran$theta
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
br=as.vector(by-bX%*%theta.source)
BtB=t(bX)%*%bX
########## Initial Estimation ############
e=by-matrixVectorMultiply(bX,theta.ini)
gamma.ini=e
gamma.ini[which(abs(gamma.ini)<5)]=0
########## Iteration ###################
Bic=matrix(0,length(Lvec),length(tauvec))
Btheta=array(0,c(length(Lvec),length(tauvec),p))
Bgamma=array(0,c(length(Lvec),length(tauvec),n))
for(i in 1:length(Lvec)){
fit.susie=NULL
delta=theta.source-theta.ini
for(v in length(tauvec):1){
error=2
iter=0
gamma=gamma.ini
gamma1=u=gamma*0
while(error>max.eps&iter<max.iter){
delta1=delta
indvalid=which(gamma1==0)
if(length(indvalid)==n){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:n,indvalid))
}
fit.cluster=center.classifying(delta,-theta.source)
delta.complement=fit.cluster$complement
delta.cluster=fit.cluster$cluster
br.complement=c(br-bX%*%delta.complement-gamma)
addbias=addbiasterm(RxyList,theta.source+delta.complement,indvalid)
XtX=BtB-Rxysum[1:p,1:p]
XtX=t(XtX)/2+XtX/2
Xty=matrixVectorMultiply(t(bX),br.complement)-Rxysum[1+p,1:p]+addbias
yty=sum((br.complement)^2)
tryCatch({
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[i],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F,residual_variance_lowerbound=1)
},error = function(e) {
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[i],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F,estimate_residual_variance=F)
})
delta.latent=coef.susie(fit.susie)[-1]*(fit.susie$pip>=pip.thres)
inddelta=which(delta.latent!=0)
Diff=generate_block_matrix(summary(fit.susie)$vars,n/diag(BtB),delta.latent)
delta=delta*0
if(length(inddelta)==1){
xtx=XtX[inddelta,inddelta]
xty=Xty[inddelta]
delta.latent[inddelta]=xty/xtx
}
if(length(inddelta)>1){
xtx=XtX[inddelta,inddelta]+ridge.diff*Diff[inddelta,inddelta]
xty=Xty[inddelta]
delta.latent[inddelta]=c(solve(xtx)%*%xty)
}
delta=delta.latent+delta.complement
theta=delta+theta.source
gamma=(by-matrixVectorMultiply(bX,theta)-u+admm.rho*gamma1)/(1+admm.rho)
gamma1=mcp(gamma+u/admm.rho,tauvec[v]/admm.rho)
u=u+admm.rho*(gamma-gamma1)
iter=iter+1
if(iter>5){
error=sqrt(sum((delta-delta1)^2))
}
}
e=by-matrixVectorMultiply(bX,theta)-gamma
vare=sum(e^2)/(length(indvalid)-length(inddelta))
Bic[i,v]=log(vare)+(log(n)+log(p)*ebic.delta)/n*length(inddelta)+(1+ebic.gamma)*log(n)*(n-length(indvalid))/n
Btheta[i,v,]=theta
Bgamma[i,v,]=gamma1
}
}
############################### final estimate ##########################
istar=bimin(Bic)[1]
vstar=bimin(Bic)[2]
theta=Btheta[istar,vstar,]
gamma=Bgamma[istar,vstar,]
delta=theta.source-theta
fit.susie=NULL
error=2
iter=0
gamma1=u=0*by
while(error>max.eps&iter<max.iter){
delta1=delta
indvalid=which(gamma1==0)
if(length(indvalid)==n){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:n,indvalid))
}
fit.cluster=center.classifying(delta,-theta.source)
delta.complement=fit.cluster$complement
delta.cluster=fit.cluster$cluster
br.complement=c(br-bX%*%delta.complement-gamma)
addbias=addbiasterm(RxyList,theta.source+delta.complement,indvalid)
XtX=BtB-Rxysum[1:p,1:p]
XtX=t(XtX)/2+XtX/2
Xty=matrixVectorMultiply(t(bX),br.complement)-Rxysum[1+p,1:p]+addbias
yty=sum((br.complement)^2)
tryCatch({
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[istar],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F,residual_variance_lowerbound=1)
},error = function(e) {
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[istar],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F,estimate_residual_variance=F)
})
delta.latent=coef.susie(fit.susie)[-1]*(fit.susie$pip>=pip.thres)
inddelta=which(delta.latent!=0)
Diff=generate_block_matrix(summary(fit.susie)$vars,n/diag(BtB),delta.latent)
delta=delta*0
if(length(inddelta)==1){
xtx=XtX[inddelta,inddelta]
xty=Xty[inddelta]
delta.latent[inddelta]=xty/xtx
}
if(length(inddelta)>1){
xtx=XtX[inddelta,inddelta]+ridge.diff*Diff[inddelta,inddelta]
xty=Xty[inddelta]
delta.latent[inddelta]=c(solve(xtx)%*%xty)
}
delta=delta.latent+delta.complement
theta=delta+theta.source
gamma=(by-matrixVectorMultiply(bX,theta)-u+admm.rho*gamma1)/(1+admm.rho)
gamma1=mcp(gamma+u/admm.rho,tauvec[vstar]/admm.rho)
u=u+admm.rho*(gamma-gamma1)
iter=iter+1
if(iter>5){
error=sqrt(sum((delta-delta1)^2))
}
}
############################### inference #########################
names(delta)=colnames(bX)
theta=theta.source+delta
res=gamma1*byse1
names(res)=rownames(bX)
ThetaList=DeltaList=matrix(0,sampling.time,p)
colnames(ThetaList)=colnames(DeltaList)=colnames(bX)
cat("Bootstrapping process:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
j=1
while(j<=sampling.time) {
setTxtProgressBar(pb, j)
indicator <- FALSE
tryCatch({
indj=sample(n,0.5*n,replace=F)
nj=length(indj)
bXj=bX[indj,]
byj=by[indj]
bXsej=bXse[indj,]
bysej=byse[indj]
brj=br[indj]
thetaj=theta*runif(length(theta),0.95,1.05)
RxyListj=RxyList[indj,,]
Rxyallj=biasterm(RxyList=RxyListj,c(1:nj))
gammaj=gamma[indj]
gamma1j=gamma1[indj]
uj=gammaj-gamma1j
indvalidj=which(gamma1j==0)
fit.susiej=fit.susie
deltaj=theta.source-thetaj
BtBj=matrixMultiply(t(bXj),bXj)
for(iterj in 1:sampling.iter){
indvalidj=which(gamma1j==0)
if(length(indvalidj)==nj){
Rxysumj=Rxyallj
}else{
Rxysumj=Rxyallj-biasterm(RxyList=RxyListj,setdiff(1:nj,indvalidj))
}
fit.clusterj=center.classifying(deltaj,-theta.source)
delta.complementj=fit.clusterj$complement
delta.clusterj=fit.clusterj$clusterj
br.complementj=c(brj-bXj%*%delta.complementj-gammaj)
addbiasj=addbiasterm(RxyListj,theta.source+delta.complementj,indvalidj)
XtXj=BtBj-Rxysumj[1:p,1:p]
XtXj=t(XtXj)/2+XtXj/2
Xtyj=matrixVectorMultiply(t(bXj),br.complementj)-Rxysumj[1+p,1:p]+addbiasj
ytyj=sum((br.complementj)^2)
fit.susiej=susie_suff_stat(XtX=XtXj,Xty=Xtyj,yty=ytyj,L=Lvec[istar],n=length(indvalidj),estimate_prior_method="EM",residual_variance=1,s_init=fit.susiej,standardize=F,max_iter=15,intercept=F,residual_variance_lowerbound=1)
delta.latentj=coef.susie(fit.susiej)[-1]*(fit.susiej$pip>=pip.thres)
inddeltaj=which(delta.latentj!=0)
Diffj=generate_block_matrix(summary(fit.susiej)$vars,nj/diag(BtBj),delta.latentj)
deltaj=deltaj*0
if(length(inddeltaj)==1){
xtxj=XtXj[inddeltaj,inddeltaj]
xtyj=Xtyj[inddeltaj]
delta.latentj[inddeltaj]=xtyj/xtxj
}
if(length(inddeltaj)>1){
xtxj=XtXj[inddeltaj,inddeltaj]+ridge.diff*Diffj[inddeltaj,inddeltaj]
xtyj=Xtyj[inddeltaj]
delta.latentj[inddeltaj]=c(solve(xtxj)%*%xtyj)
}
deltaj=delta.latentj+delta.complementj
thetaj=deltaj+theta.source
gammaj=(byj-matrixVectorMultiply(bXj,thetaj)-uj+admm.rho*gamma1j)/(1+admm.rho)
gamma1j=mcp(gammaj+uj/admm.rho,tauvec[vstar]/admm.rho)
uj=uj+admm.rho*(gammaj-gamma1j)
}
ThetaList[j,]=theta.source+deltaj
DeltaList[j,]=deltaj
j=j+1
}, error = function(e) {
# Error handling block
cat("Error occurred: ", e$message, "\n")
indicator <<- TRUE  # Set indicator to TRUE if an error occurs
j <<- j - 1  # Decrement the iteration counter to retry
})
if (indicator) {
next  # Retry the current iteration
}
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
A$gamma=res
A$delta.se=delta.se
A$theta.se=theta.se
A$delta.cov=delta.cov
A$theta.cov=theta.cov
A$reliability.adjust=r
A$susie.delta=fit.susie
A$delta.latent=delta.latent
A$delta.list=DeltaList
A$theta.list=ThetaList
A$Bic=Bic
A$L.optimal=Lvec[istar]
A$tau.optimal=tauvec[vstar]
return(A)

}
