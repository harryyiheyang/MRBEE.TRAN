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

center.classifying=function(delta,theta.source){
p=length(delta)
complement=delta*0
for(i in 1:p){
s=which.min(c(abs(delta[i]),abs(delta[i]+theta.source[i])))
complement[i]=ifelse(s==1,0,-theta.source[i])
}
return(complement)
}
