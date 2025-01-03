IVweight_CML=function(byse,bXse,Rxy){
Omegaxy=solve(Rxy)
bZse=cbind(bXse,byse)
p=dim(bZse)[2]
n=dim(bZse)[1]
OmegaxyList=RxyList=array(0,c(n,p,p))
for(i in 1:n){
s=bZse[i,]
RxyList[i,,]=t(t(Rxy)*s)*s
s1=1/s
OmegaxyList[i,,]=t(t(Omegaxy)*s1)*s1
}
return(list(RxyList=RxyList,OmegaxyList=OmegaxyList))
}

IVweight=function(byse,bXse,Rxy){
bZse=cbind(bXse,byse)
p=dim(bZse)[2]
n=dim(bZse)[1]
RxyList=array(0,c(n,p,p))
for(i in 1:n){
s=bZse[i,]
RxyList[i,,]=t(t(Rxy)*s)*s
}
return(RxyList)
}

imrpdetect=function(x,theta,RxyList,indvalid,var.est="robust",FDR=T,adjust.method="Sidak"){
p=length(theta)
if(var.est=="robust"){
varx=stats::mad(x[indvalid])^2
}
if(var.est=="variance"){varx=stats::var(x[indvalid])}
if(var.est=="ordinal"){
varx=x*0
for(i in 1:length(x)){
varx[i]=c(RxyList[i,p+1,p+1]+t(theta)%*%RxyList[i,1:p,1:p]%*%theta-2*sum(theta*RxyList[i,p+1,1:p]))
}
}
pv=stats::pchisq(x^2/varx,1,lower.tail=F)
if(FDR==T){
pv=FDRestimation::p.fdr(pvalues=pv,adjust.method=adjust.method)$fdrs
}
return(as.vector(pv))
}

validadj <- function(vector1, vector2, tau) {
diff <- length(vector2) / length(vector1)
if (diff < tau) {
missing_indices <- setdiff(1:length(vector1), vector2)
sorted_missing_indices <- missing_indices[order(vector1[missing_indices])]
num_to_add <- ceiling(tau * length(vector1)) - length(vector2)
vector2 <- c(vector2, sorted_missing_indices[1:num_to_add])
}
return(vector2)
}

biasterm=function(RxyList,indvalid){
X=RxyList[1,,]*0
for(i in indvalid){
X=X+RxyList[i,,]
}
return(X)
}

reliability.adj.uv=function(bx,bxse,Theta="identity",thres=0.7){
if(Theta[1]=="identity"){
total.var=mean(bx^2)
error.var=mean(bxse^2)
reliability=(total.var-error.var)/total.var
r=1
if(reliability<thres){
r=total.var/error.var*(1-thres)
}
r=sqrt(r)
}else{
r=1
Theta=as.matrix(Theta)
total.var=mean(bx*(Theta%*%bx))
error.var=mean(bxse^2)
reliability=(total.var-error.var)/total.var
if(reliability<thres){
r=total.var/error.var*(1-thres)
}
r=sqrt(r)
}
return(r)
}

reliability.adj=function(bX,bXse,Theta="identity",thres=0.7){
if(Theta[1]=="identity"){
p=ncol(bX)
r=rep(1,p)
total.var=colMeans(bX^2)
error.var=colMeans(bXse^2)
reliability=(total.var-error.var)/total.var
ind=which(reliability<thres)
if(length(ind)>0){
r[ind]=total.var[ind]/error.var[ind]*(1-thres)
}
r=sqrt(r)
}else{
p=ncol(bX)
r=rep(1,p)
m=length(bX[,1])
Theta=as.matrix(Theta)
total.var=as.vector(diag(t(bX)%*%Theta%*%bX)/m)
error.var=colMeans(bXse^2)
reliability=(total.var-error.var)/total.var
ind=which(reliability<thres)
if(length(ind)>0){
r[ind]=total.var[ind]/error.var[ind]*(1-thres)
}
r=sqrt(r)
}
return(r)
}

generate_D_matrix <- function(s, sign_vec) {
p <- length(s)
if (length(sign_vec) != p) {
stop("Length of sign_vec must match length of s.")
}

if (p == 1) {
D <- 0
} else {
num_pairs <- p*(p-1)/2
D_all <- matrix(0, nrow = num_pairs, ncol = p)

row_idx <- 1
for (i in 1:(p-1)) {
for (j in (i+1):p) {
D_all[row_idx, i] <-  sign_vec[i] / s[i]
D_all[row_idx, j] <- -sign_vec[j] / s[j]
row_idx <- row_idx + 1
}
}

D <- t(D_all)%*%D_all
}

return(D)
}

generate_block_matrix <- function(vars_df, s, theta) {
ind=which(theta==0)
vars_df$cs[which(vars_df$variable%in%ind)]=-1
concerned_vars <- vars_df[vars_df$cs != -1, ]
cs_values <- unique(concerned_vars$cs)
max_var_index <- max(vars_df$variable)
final_matrix <- matrix(0, nrow = max_var_index, ncol = max_var_index)
for (cs_val in cs_values) {
group_vars <- concerned_vars$variable[concerned_vars$cs == cs_val]
if (length(group_vars) > 1) {
group_s <- s[group_vars]
D <- generate_D_matrix(group_s,sign(theta[group_vars]))
final_matrix[group_vars, group_vars]=D
}
}
return(final_matrix)
}

group.pip.filter=function(pip.summary,xQTL.cred.thres=0.95,xQTL.pip.thres=0.1){
ind=which(pip.summary$cs>0)
if(length(ind)>0){
J=max(pip.summary$cs[ind])
pip.summary$cs.pip=pip.summary$variable_prob
for(i in 1:J){
indi=which(pip.summary$cs==i)
summaryi=pip.summary[indi,]
pip.cred=sum(summaryi$variable_prob)
pip.summary$cs.pip[indi]=pip.cred
}
ind.keep=which(pip.summary$cs.pip>=xQTL.cred.thres&pip.summary$variable_prob>=xQTL.pip.thres)
cs=pip.summary$cs
cs.pip=pip.summary$cs.pip
cs->cs[pip.summary$variable]
cs.pip->cs.pip[pip.summary$variable]
cs[which(cs==-1)]=0
}else{
ind.keep=NULL
cs=pip.summary$cs.pip*0
cs.pip=pip.summary$cs.pip*0
}
return(list(ind.keep=pip.summary$variable[ind.keep],cs=cs,cs.pip=cs.pip,result=pip.summary))
}

colSD=function(A){
a=A[1,]
for(i in 1:ncol(A)){
a[i]=sd(A[,i])
}
return(a)
}

sloepse=function(theta.source,theta.source.cov,bX,by,bXse,byse,Rxy){
n=length(byse)
p=ncol(bXse)
bzse=byse
Rxx=Rxy[1:p,1:p]
for(i in 1:n){
bxi=as.vector(bX[i,])
bxsei=as.vector(bXse[i,])
Covi=t(t(Rxx)*bxsei)*bxsei
bzse[i]=sum(theta.source*(Covi%*%theta.source))+sum(bxi*(theta.source.cov%*%bxi))
}
return(bzse)
}

center.classifying=function(theta,theta.source){
p=length(theta)
complement=cluster=theta*0
for(i in 1:p){
s=which.min(c(abs(theta[i]),abs(theta[i]-theta.source[i])))
complement[i]=ifelse(s==1,0,theta.source[i])
cluster[i]=ifelse(s==1,1,2)
}
return(list(complement=complement,cluster=cluster))
}

addbiasterm=function(Rxylist,theta.complement,indvalid){
g=theta.complement*0
p=length(theta.complement)
for(i in indvalid){
g=g+as.vector(Rxylist[i,1:p,1:p]%*%theta.complement)
}
return(g)
}

mcp=function(x,lam,a=3){
b=abs(x)
z=soft(x,lam)/(1-1/a)
z[which(b>(a*lam))]=x[which(b>(a*lam))]
return(z)
}

soft=function(a,b){
c=abs(a)-b
c[c<0]=0
c=c*sign(a)
return(c)
}

bimin=function(mat){
min_element=min(mat)
min_indices=which(mat == min_element, arr.ind = TRUE)
if (nrow(min_indices) > 1) {
min_indices=min_indices[nrow(min_indices), ]
}
return(min_indices)
}


CML_Bupdate=function(by,bX,theta,OmegaxyList){
n=length(by)
bXest=bX*0
bZ=cbind(bX,by)
p=dim(bX)[2]
for(i in 1:n){
z=bZ[i,]
G=cbind(diag(p),theta)
S1=matrixMultiply(G,OmegaxyList[i,,])
S2=matrixMultiply(S1,t(G))
S2=matrixInverse(S2)
hatz=matrixVectorMultiply(S2,matrixVectorMultiply(S1,z))
bXest[i,]=hatz
}
return(bXest)
}

CML_Bupdate_uv=function(by,bX,theta,OmegaxyList){
n=length(by)
bXest=bX*0
bZ=cbind(bX,by)
p=1
for(i in 1:n){
z=bZ[i,]
G=c(1,theta)
S1=c(OmegaxyList[i,,]%*%G)
S2=sum(S1*G)
hatz=sum(S1*z)/S2
bXest[i]=hatz
}
return(bXest)
}

CML_thetaupdate=function(by,bX,bXest,OmegaxyList,indvalid,indtheta,glatent,Diff,ridge.diff){
p=dim(bX)[2]
varomega=OmegaxyList[indvalid,p+1,p+1]
g1=matrixVectorMultiply(t(bXest[indvalid,indtheta]),by[indvalid]*varomega)
H=matrixMultiply(t(bXest[indvalid,indtheta]),bXest[indvalid,indtheta]*varomega)
theta=bX[1,]*0
theta[indtheta]=solve(H+ridge.diff*Diff[indtheta,indtheta])%*%(glatent[indtheta]+g1)
return(theta)
}

CML_thetaupdate_uv=function(by,bX,bXest,OmegaxyList,indvalid,indtheta,glatent){
p=dim(bX)[2]
varomega=OmegaxyList[indvalid,p+1,p+1]
g1=sum(bXest[indvalid,indtheta]*by[indvalid]*varomega)
H=sum(bXest[indvalid,indtheta]^2*varomega)
theta=bX[1,]*0
theta[indtheta]=(glatent[indtheta]+g1)/H
return(theta)
}
