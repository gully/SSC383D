# Authors: Gully and Amanda
# Date: Feb. 26, 2014
# Desc: Cheese example

library(bayesm)

#Taken from the examples in bayesm

#Outline:
#1) Get the data
#2) Inspect the dimensions, summarize
#3a) Wrangle the data
#3b) Example of processing for use with rhierLinearModel
#4) Run each individual regression and store results


#------------------------------
#1) Get the data
#------------------------------
data(cheese)

#------------------------------
#2) Inspect the dimensions, summarize
#------------------------------
#Find the domain and range
cat(" Quantiles of the Variables ",fill=TRUE)
mat=apply(as.matrix(cheese[,2:4]),2,quantile)
print(mat)

#------------------------------
#3) Wrangle the data
#------------------------------

#Find the unique store entries, count them
retailer=levels(cheese$RETAILER)
nreg=length(retailer) 
nvar=3

#Arrange the data into a list of the 88 stores,
#	with a list of ( y [1 x n_entries], X [n_entries x 3] for each store.
#	each store has a unique number of entries,
#	typically 61, 68, or 52.  This is like number of weeks.
#Also, wrangle the data into ln() form:
#	y = ln(Q)
#	X = [1, disp, ln(P)] week0
#		...
#	X = [1, disp, ln(P)] week61
regdata=NULL
for (reg in 1:nreg) {
y=log(cheese$VOLUME[cheese$RETAILER==retailer[reg]])
iota=c(rep(1,length(y)))
X=cbind(iota,cheese$DISP[cheese$RETAILER==retailer[reg]],
log(cheese$PRICE[cheese$RETAILER==retailer[reg]]))
regdata[[reg]]=list(y=y,X=X)
}
Z=matrix(c(rep(1,nreg)),ncol=1)
nz=ncol(Z)


#------------------------------
#4) Run each individual regression and store results
#------------------------------

#Make an 88 x 3 matrix for the least squares coefficients
lscoef=matrix(double(nreg*nvar),ncol=nvar)

#Loop over each store:
#	Fit the problem: y

reg=1
plot(regdata[[reg]]$X[,3],regdata[[reg]]$y, xlab='ln(P)', ylab='ln(Q)', xlim=c(0.5, 1.5), ylim=c(6.0, 10.0), main='Price Elasticity of demand')

for (reg in 1:nreg) {
lsfitList=lsfit(regdata[[reg]]$X,regdata[[reg]]$y,intercept=FALSE)
coef=lsfitList$coef
residuals=lsfitList$residuals

points(regdata[[reg]]$X[,3],regdata[[reg]]$y, col=reg)

if (var(regdata[[reg]]$X[,2])==0) { lscoef[reg,1]=coef[1]; lscoef[reg,3]=coef[2]}
else {lscoef[reg,]=coef }
}

R=2000
Data=list(regdata=regdata,Z=Z)
Mcmc=list(R=R,keep=1)
set.seed(66)
out=rhierLinearModel(Data=Data,Mcmc=Mcmc)
cat("Summary of Delta Draws",fill=TRUE)
summary(out$Deltadraw)
cat("Summary of Vbeta Draws",fill=TRUE)
summary(out$Vbetadraw)
if(0){
#
# plot hier coefs
plot(out$betadraw)
}
