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

#' @return A list containing the estimated causal effect, its covariance, and pleiotropy.
#' @importFrom susieR susie_suff_stat coef.susie susie
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply matrixListProduct
#' @importFrom MASS rlm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'
MRBEE_IPOD_SuSiE=function(by,bX,byse,bXse,Rxy,tauvec=seq(3,30,3),admm.rho=3,Lvec=c(1:6),ebic.theta=1,ebic.gamma=2,susie.iter=200,pip.thres=0.3,max.iter=50,max.eps=1e-4,reliability.thres=0.8,ridge.diff=100,sampling.time=100,sampling.iter=10){
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
XtX=matrixMultiply(t(bX),bX)
XtX=t(XtX)/2+XtX/2
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy)
Rxyall=biasterm(RxyList=RxyList,c(1:n))
########## Initial Estimation ############
fit=susie(y=by,X=bX,L=5,max_iter=1000)
theta.ini=coef.susie(fit)[-1]*(fit$pip>0.3)
theta=theta.ini
theta1=10000
e=c(by-bX%*%theta)
gamma.ini=e
gamma.ini[which(abs(gamma.ini)<5)]=0
########## Iteration ###################
Bic=matrix(0,length(Lvec),length(tauvec))
Btheta=array(0,c(length(Lvec),length(tauvec),p))
Bgamma=array(0,c(length(Lvec),length(tauvec),n))
tauvec=sort(tauvec)
for(i in 1:length(Lvec)){
fit.susie=NULL
theta=theta.ini
for(v in length(tauvec):1){
error=2
iter=0
gamma=gamma.ini
gamma1=u=gamma*0
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==n){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:n,indvalid))
}
e=by-gamma
Xty=matrixVectorMultiply(t(bX),e)
yty=sum(e^2)
tryCatch({
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[i],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F,residual_variance_lowerbound=1)
},error = function(e) {
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[i],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F,estimate_residual_variance=F)
})
theta=coef.susie(fit.susie)[-1]*(fit.susie$pip>=pip.thres)
indtheta=which(theta!=0)
Diff=generate_block_matrix(summary(fit.susie)$vars,length(indvalid)/diag(XtX),theta)
if(length(indtheta)==1){
xtx=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
xtx=XtX[indtheta,indtheta]+ridge.diff*Diff[indtheta,indtheta]-Rxysum[indtheta,indtheta]
xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=c(solve(xtx)%*%xty)
}
gamma=(by-matrixVectorMultiply(bX,theta)-u+admm.rho*gamma1)/(1+admm.rho)
gamma1=mcp(gamma+u/admm.rho,tauvec[v]/admm.rho)
u=u+admm.rho*(gamma-gamma1)
iter=iter+1
if(iter>5){
error=sqrt(sum((theta-theta1)^2))
}
}
e=e-matrixVectorMultiply(bX,theta)
vare=sum(e^2)/(length(indvalid)-length(indtheta))
Bic[i,v]=log(vare)+(log(n)+log(p)*ebic.theta)/n*length(indtheta)+(1+ebic.gamma)*log(n)*(n-length(indvalid))/n
Btheta[i,v,]=theta
Bgamma[i,v,]=gamma1
}
}
############################### final estimate ##########################
istar=bimin(Bic)[1]
vstar=bimin(Bic)[2]
theta=Btheta[istar,vstar,]
error=2
iter=0
gamma=Bgamma[istar,vstar,]
gamma1=u=gamma*0
fit.susie=NULL
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==n){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:n,indvalid))
}
e=by-gamma
Xty=matrixVectorMultiply(t(bX),e)
yty=sum(e^2)
tryCatch({
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[istar],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F,residual_variance_lowerbound=1)
},error = function(e) {
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[istar],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F,estimate_residual_variance=F)
})
theta=coef.susie(fit.susie)[-1]*(fit.susie$pip>=pip.thres)
indtheta=which(theta!=0)
Diff=generate_block_matrix(summary(fit.susie)$vars,n/diag(XtX),theta)
if(length(indtheta)==1){
xtx=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
xtx=XtX[indtheta,indtheta]+ridge.diff*Diff[indtheta,indtheta]-Rxysum[indtheta,indtheta]
xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=c(solve(xtx)%*%xty)
}
gamma=(by-matrixVectorMultiply(bX,theta)-u+admm.rho*gamma1)/(1+admm.rho)
gamma1=mcp(gamma+u/admm.rho,tauvec[vstar]/admm.rho)
u=u+admm.rho*(gamma-gamma1)
iter=iter+1
if(iter>5){
error=sqrt(sum((theta-theta1)^2))
}
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
XtXj=matrixMultiply(t(bXj),bXj)
XtXj=t(XtXj)/2+XtXj/2
thetaj=theta*runif(length(theta),0.95,1.05)
RxyListj=RxyList[indj,,]
Rxyallj=biasterm(RxyList=RxyListj,c(1:nj))
gammaj=gamma[indj]
gamma1j=uj=gammaj*0
indvalidj=which(gamma1j==0)
fit.susiej=fit.susie
for(iterj in 1:sampling.iter){
indvalidj=which(gamma1j==0)
if(length(indvalidj)==nj){
Rxysumj=Rxyallj
}else{
Rxysumj=Rxyallj-biasterm(RxyList=RxyListj,setdiff(1:nj,indvalidj))
}
ej=byj-gammaj
Xtyj=matrixVectorMultiply(t(bXj),ej)
ytyj=sum(ej^2)
fit.susiej=susie_suff_stat(XtX=XtXj,Xty=Xtyj,yty=ytyj,L=Lvec[istar],n=length(indvalidj),estimate_prior_method="EM",residual_variance=1,s_init=fit.susiej,standardize=F,max_iter=15,intercept=F,residual_variance_lowerbound=1)
thetaj=coef.susie(fit.susiej)[-1]*(fit.susiej$pip>=pip.thres)
indthetaj=which(thetaj!=0)
Diffj=generate_block_matrix(summary(fit.susiej)$vars,nj/diag(XtXj),thetaj)
if(length(indthetaj)==1){
xtxj=XtXj[indthetaj,indthetaj]-Rxysumj[indthetaj,indthetaj]
xtyj=Xtyj[indthetaj]-Rxysumj[indthetaj,p+1]
thetaj[indthetaj]=xtyj/xtxj
}
if(length(indthetaj)>1){
xtxj=XtXj[indthetaj,indthetaj]+ridge.diff*Diffj[indthetaj,indthetaj]-Rxysumj[indthetaj,indthetaj]
xtyj=Xtyj[indthetaj]-Rxysumj[indthetaj,p+1]
thetaj[indthetaj]=c(solve(xtxj)%*%xtyj)
}
gammaj=as.vector(byj-matrixVectorMultiply(bXj,thetaj)-uj+admm.rho*gamma1j)/(1+admm.rho)
gamma1j=mcp(gammaj+uj/admm.rho,tauvec[vstar]/admm.rho)
uj=uj+admm.rho*(gammaj-gamma1j)
}
ThetaList[j,]=thetaj
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
theta.se=colSD(ThetaList)*sqrt(n/length(indvalid))
theta.cov=cov(ThetaList)*n/length(indvalid)
colnames(theta.cov)=rownames(theta.cov)=names(theta.se)=colnames(bX)
A=list()
A$theta=theta
A$gamma=res
A$theta.se=theta.se
A$theta.cov=theta.cov
A$reliability.adjust=r
A$susie.theta=fit.susie
A$theta.list=ThetaList
A$Bic=Bic
A$tau.optimal=tauvec[vstar]
A$L.optimal=Lvec[istar]
return(A)

}
