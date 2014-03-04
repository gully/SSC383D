# Authors: Gully and Amanda
# Date: Feb. 26, 2014
# Desc: Cheese example

library(bayesm)
library(TeachingDemos)
#Taken from the examples in bayesm

#Outline:
#1) Get the data
#2) Inspect the dimensions, summarize
#3a) Wrangle the data
#3b) Example of processing for use with rhierLinearModel
#4) Run each individual regression and store results

#Desired outputs:


#------------------------------
#1) Get the data
#------------------------------

#First from the package, which is nearly identical
data(cheese)

#Second from the SSC383D GitHub .csv file
cheese_fn = '/Users/gully/astroML/SSC383D/SSC383D/Chapter05/examples/cheese.csv'
dat1= read.csv(cheese_fn)
str(dat1)

##Find where the first cheese data has finite value for its display
# str(cheese)
# ids= cheese$RETAILER==dat1$store
cheeseDispBool= cheese$DISP>0.0
dat1DispBool= dat1$disp == 1
#compDispBool= cheeseDispBool == dat1DispBool
##Yes, the class data is merely a binary 'yes' or 'no' for signs

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
for (reg in 1:nreg) {

#Ignoring the displays
lsfitList0=lsfit(regdata[[reg]]$X[,c(1,3)],regdata[[reg]]$y,intercept=FALSE)
coef0=lsfitList0$coef
residuals=lsfitList0$residuals

#Including the displays
lsfitList=lsfit(regdata[[reg]]$X,regdata[[reg]]$y,intercept=FALSE)
coef=lsfitList$coef

#Handle the two cases with NO displays
noDisplays= (var(regdata[[reg]]$X[,2])==0)
Dispi = regdata[[reg]]$X[,2] > 0.0
allDisplays=length(Dispi)==sum(Dispi)
if (noDisplays) {lscoef[reg,1]=coef[1]; lscoef[reg,3]=coef[2]}
else {lscoef[reg,]=coef }

old.par <- par(mfrow=c(2, 1))

#Plot ln(Q) vs. ln(P)
plot(regdata[[reg]]$X[,3],regdata[[reg]]$y, xlab='ln(P)', ylab='ln(Q)', xlim=c(0.5, 1.5), ylim=c(6.0, 10.0), main='Price Elasticity of demand')
abline(coef0[1], coef0[2], col="blue") # regression line (y~x)
abline(coef[1], coef[3], col="red") # regression line (y~x)

#Highlight the ones with a display:
points(regdata[[reg]]$X[Dispi,3],regdata[[reg]]$y[Dispi], col='green')

# plot(regdata[[reg]]$X[,3], residuals, ylim=c(-0.6, 0.6),  xlab='ln(P)', ylab='resid ln(Q)', xlim=c(0.5, 1.5), main='Residuals')
# abline(0,0, col='red')
# points(regdata[[reg]]$X[Dispi,3], residuals[Dispi], col='green')

#Histogram:
if(!noDisplays && !allDisplays){
  ymax=40
  hist(residuals[!Dispi], col=rgb(1,0,0,1/4), xlim=c(-1,1), ylim=c(0,ymax))
  hist(residuals[Dispi], col=rgb(0,1,0,1/4), add=TRUE)
  }
else{
  hist(c(-0.5, 0, 0, 0, 0.5), col=rgb(1,0,0,1/4), xlim=c(-1,1), ylim=c(0,ymax))
  }
  
par(old.par)

Sys.sleep(0.50)
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
#
# plot hier coefs
str(out)
plot(out$betadraw)
plot(out$Deltadraw)
plot(out$Vbetadraw)
plot(out$taudraw)
