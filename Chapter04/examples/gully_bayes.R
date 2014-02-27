#Author: Michael Gully-Santiago
#Date: Feb 26, 2014
#Desc: Data set for SSC383D

#Outline of the program
#1) Read in the data
#2) Wrangle the data
#3) Inspect the data and make some plots


#The design matrix looks like this:

#              age         sex         bmi          map           tc         ldl         hdl          tch          ltg         glu
#[1,]  0.038075906  0.05068012  0.06169621  0.021872355 -0.044223498 -0.03482076 -0.04340085 -0.002592262  0.019908421 -0.01764613
#[2,] -0.001882017 -0.04464164 -0.05147406 -0.026327835 -0.008448724 -0.01916334  0.07441156 -0.039493383 -0.068329744 -0.09220405
#[3,]  0.085298906  0.05068012  0.04445121 -0.005670611 -0.045599451 -0.03419447 -0.03235593 -0.002592262  0.002863771 -0.02593034
#...
#[442,]

#-----------------------------------------------------
#1) Read in the data
data(diabetes, package="BayesBridge");
#-----------------------------------------------------

#-----------------------------------------------------
#2) Wrangle the data
cov.name = colnames(diabetes$x);
y = diabetes$y;
X = diabetes$x;
#-----------------------------------------------------


#-----------------------------------------------------
#3) Inspect the data
#-----------------------------------------------------
plot(X[,1], X[,10], xlab='Age', ylab='Glucose')

plot(X[,1], X[,3], xlab='Age', ylab='BMI')

plot(X[,1], X[,4], xlab='Age', ylab='Blood Pressure')



#Which single column is the biggest predictor of effect?
plot(X[,10], y, ylab='Effect')
#BMI is stronger than age
#7 is anticorrelated with effect
#9 is decent


#-----------------------------------------------------
#4) More wrangling
#-----------------------------------------------------

# Center the data.
y = y - mean(y);
mX = colMeans(X);
for(i in 1:442){
X[i,] = X[i,] - mX;
}

#-----------------------------------------------------
#5) Advanced stuff
#-----------------------------------------------------
# Expectation maximization.
bridge.EM(y, X, 0.5, 1.0, 1e8, 1e-9, 30, use.cg=TRUE);

out1=trace.beta(y, X, alpha=0.5, ratio.grid=exp(seq(-20,20,0.1)),
tol=1e-9, max.iter=20, use.cg=FALSE, plot.it=FALSE);

gb = bridge.reg(y, X, nsamp=10000, alpha=0.5,
sig2.shape=0.0, sig2.scale=0.0, nu.shape=2.0, nu.rate=2.0);

# plot(gb$tau)
# plot(gb$alpha)
# plot(gb$u)
# plot(gb$w)
# plot(gb$shape)
# plot(gb$sig2)