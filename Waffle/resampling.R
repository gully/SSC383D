library(mosaic)

# Nonparametric bootstrap
heights = read.csv("heights.csv")

plot(SHGT ~ MHGT, data=heights)

lmM = lm(SHGT ~ MHGT, data=heights)
lmF = lm(SHGT ~ FHGT, data=heights)

N = nrow(heights)

NMC = 500
betasave = matrix(0, NMC, 2)
for(i in 1:NMC) {
	keep = sample(1:N, N, replace=TRUE)
	lmstar = lm(SHGT ~ MHGT, data=heights[keep,])
	#plot(SHGT ~MHGT, data=heights[keep,],xlim=c(50,75),ylim=c(60,80))
	#abline(lmstar)
	betasave[i,] = coef(lmstar)
}

hist(betasave[,1])
hist(betasave[,2])


########
# Residual-resampling bootstrap
########

chym = read.csv("chymotrypsin.csv")
plot(Rate ~ Conc, data=chym)

# Michaelis-Menten kinetics

###impute the curve given points and parameters
mmpredict = function(x, vmax, km) {
	vmax*x/{km+x}
}

#x=concentrations
#y=observed rates at each concentration
target = function(theta, x, y) {
	vmax = exp(theta[1])
	km = exp(theta[2])
	ypred = mmpredict(x, vmax, km) 
	sum(0.5*{y-ypred}^2)#non-linear least squares; optimization is better w/o boxed constraints so use exponent 
}

##should actually find gradients if care about this
mymax = optim(c(0,0), target, x = chym$Conc, y = chym$Rate) #optim defaults to Nelder-Mead if not given gradients 
theta = mymax$par #mymax is on log scale

rgrid = seq(0,0.5, length=100) #construct discreet grid along x-axis
plot(Rate ~ Conc, data=chym)
lines(rgrid, mmpredict(rgrid, exp(theta[1]), exp(theta[2]))) #theta1 and theta2 are optimized parameters
##but we don't know anything about uncertainty! what are the error bars going to be? 

### Now bootstrap
yhat = mmpredict(chym$Conc, exp(theta[1]), exp(theta[2])) #vector of fitted values corresponding to each concentration point
eps = chym$Rate - yhat #residual vector
N = nrow(chym) #sample size 

NMC = 250 #number of bootstrap samples
thetasave = matrix(0, NMC, 2)

for(i in 1:NMC) {
	estar = sample(eps, N, replace=TRUE) #take sample size N with replacement of residuals
	ystar = yhat + estar #residual fitted values + resampled residuals= fake dataset
	mymax = optim(c(0,0), target, x = chym$Conc, y = ystar)
	thetastar = mymax$par
	plot(chym$Conc, ystar,ylim=c(1,2.1))
	lines(rgrid, mmpredict(rgrid, exp(thetastar[1]), exp(thetastar[2]))) #fit curve to fake data
	thetasave[i,] = thetastar
}

# Inspect the sampling distributions
hist(exp(thetasave[,1]))
hist(exp(thetasave[,2]))
plot(Rate ~ Conc, data=chym)
lines(rgrid, mmpredict(rgrid, exp(theta[1]), exp(theta[2])))




## Parametric bootstrap
brca = read.csv('brca.csv', header=TRUE)

N = nrow(brca)

glm1 = glm(recall ~ factor(radiologist) + age5059 + age6069 + age70plus + familyhistory + biopsurg + symptoms + premeno + postmenohormone + postmenounknown + previousmammogram + density1 + density3 + density4, data=brca, family=binomial)

###are there any obvious systematic covariates we're missing??
##what covariates influence the 2 x 2 table? 

#marginal: recall | covariates
## cancer | recall, covariates
##How does R calculate the standard errors? 
		#Beta-hat converges to N(Beta*, I^(-1)(theta))
		# I^(-1) is inverse fisher information matrix
		## low curvature=more precision on mode
		##high curvature=less precision on mode


NMC = 500
betasave = matrix(0, NMC, length(coef(glm1)))
wpred = fitted(glm1) #model's predicted probablity of recall
##given fitted values, let's simulate fake data from the parametric model I assume generated my data in the first place
for(i in 1:NMC) {
	ystar = rbinom(N, 1, wpred) #fake data; 0's and 1's
	glm.boot = glm(ystar ~ factor(radiologist) + age5059 + age6069 + age70plus + familyhistory + biopsurg
		+ symptoms + premeno + postmenohormone + postmenounknown + previousmammogram + density1 + density3
		+ density4, data=brca, family=binomial)
	betasave[i,] = coef(glm.boot)
}
colnames(betasave) = names(coef(glm1))
hist(betasave[,1]) #intercept
hist(betasave[,2]) 

boot.se = apply(betasave, 2, sd) #bootstrap standard errors
asymp.se = sqrt(diag(vcov(glm1)))

cbind(boot.se, asymp.se) #compare standard errors 

hist(betasave[,16],50) #density classification 3; giant bootstrap SE -> log odds of -15 -> categories are perfectly separable-->Gaussian assumption is probably wrong



##Simple resampling technique=permutation testing
####

gardasil = read.csv('gardasil.csv')

head(gardasil)

#glm1 = glm(Completed ~ AgeGroup + InsuranceType, data=gardasil, family='binomial')
#glm2 = glm(Completed ~ AgeGroup + InsuranceType + Location, data=gardasil, family='binomial')
lm1 = lm(Completed ~ AgeGroup + InsuranceType, data=gardasil)
lm2 = lm(Completed ~ AgeGroup + InsuranceType + Location, data=gardasil)

rsq1= 1-var(resid(lm1))/var(gardasil$Completed)

NMC = 1000
rsq=rep(0,NMC)
for (i in 1:NMC){
	lm2=lm(Completed ~ AgeGroup + InsuranceType + shuffle(Location), data=gardasil)
	rsq[i] = 1 - var(resid(lm2))/var(gardasil$Completed) 
}

hist(rsq) #what is the bump in r^2 by shuffling each predictor? 
##is there an increase in r^2 b/c of junk predictors? '
##similar to F-test ( but F-test is dependent on parametric assumptions)

