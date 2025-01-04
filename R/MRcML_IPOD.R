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
#' @param tau A scale of tuning parameter in MCP. Defaults to \code{9}.
#' @param admm.rho A scale of tuning parameter in ADMM. Defaults to \code{3}.
#' @param max.iter Maximum number of iterations for causal effect estimation. Defaults to 10.
#' @param max.eps Tolerance for stopping criteria. Defaults to 1e-4.
#' @param adjust.method Method for estimating q-value. Defaults to "Sidak".
#' @param sampling.time A scale of number of subsampling in estimating the standard error.
#' @param sampling.iter A scale of iteration in subsampling in estimating the standard error.
#' @param theta.ini An initial estimate of causal effect.
#'
#' @return A list containing the estimated causal effect, its covariance, and pleiotropy.
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply matrixListProduct matrixInverse
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'
MRcML_IPOD=function(by,bX,byse,bXse,Rxy,tau=9,admm.rho=3,max.iter=10,max.eps=1e-4,sampling.time=100,sampling.iter=10,theta.ini){
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
gamma.ini=as.vector(by-bX%*%theta.ini)
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
g2=theta*0
for(i in indvalid){
s1=bX[i,]-bXest[i,]
s2=OmegaxyList[i,1:p,p+1]
s3=outer(s1,bXest[i,])
g2=g2+matrixVectorMultiply(s3,s2)
}
theta=CML_thetaupdate(by=by-gamma,bX=bX,bXest=bXest,OmegaxyList=OmegaxyList,indvalid=indvalid,indtheta=1:p,Diff=0*diag(p),ridge.diff=0,glatent=g2)
gamma=(by-matrixVectorMultiply(bXest,theta)-u+admm.rho*gamma1)/(1+admm.rho)
gamma1=mcp(gamma+u/admm.rho,tau/admm.rho)
u=u+admm.rho*(gamma-gamma1)
iter=iter+1
if(iter>5){
error=sqrt(sum((theta-theta1)^2))
}
evec[iter]=sqrt(sum((theta-theta1)^2))
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
for(iterj in 1:sampling.iter){
indvalidj=which(gamma1j==0)
bXestj=CML_Bupdate(by=byj-gammaj,bX=bXj,theta=thetaj,OmegaxyList=OmegaxyListj)
g2j=theta*0
for(i in indvalidj){
s1j=bXj[i,]-bXestj[i,]
s2j=OmegaxyListj[i,1:p,p+1]
s3j=outer(s1j,bXestj[i,])
g2j=g2j+matrixVectorMultiply(s3j,s2j)
}
thetaj=CML_thetaupdate(by=byj-gammaj,bX=bXj,bXest=bXestj,OmegaxyList=OmegaxyListj,indvalid=indvalidj,indtheta=1:p,Diff=0*diag(p),ridge.diff=0,glatent=g2j)
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
A$bXest=bXest
A$theta.list=ThetaList
A$algorithm.error=evec
return(A)

}
