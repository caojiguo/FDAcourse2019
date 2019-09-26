# R codes for Estimating Constraint functions

library('fda')
# Example when estiamting a function with the positive constraint.


# Vancouver precipitation data
VancPrec = CanadianWeather$dailyAv[,'Vancouver','Precipitation.mm']

length(VancPrec)
day = 1:365
# weather data 

fbasis = create.fourier.basis(c(1,365),365)

# define the harmonic acceleration differential operator
# L = D^3 f(t) + 0 * D^2 f(t) + w^2 * D^f(t) + 0 * f(t)
# It will not penalize a+bt, cos(wt) and sin(wt)
# w = 2pi/period of the process
# c(0,(2*pi/365)^2,0) is the vector of coefficients to 
# the lower derivatives
harmLfd = vec2Lfd(c(0,(2*pi/365)^2,0), c(1, 365))

# Tell 3 information
# 1:fbasis - Type of basis functions 
# 2: harmLfd - Define the roughness penalty
# 3: lambda - value of the smoothing parameter
harmpar = fdPar(fbasis,harmLfd,lambda=1e4)

# smoothing without considering the positive contraint;
precfd1 = smooth.basis(day,VancPrec,harmpar)

# smoothing when considering the positive contraint;
posprecfd1 = smooth.pos(day,VancPrec,harmpar)

# plot the fitte curves
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(day,VancPrec,col=2,xlab='day',ylab='precipitation',cex.lab=1.5,cex.axis=1.5)
# plot the estimated function without considering the positive constraints
lines(precfd1$fd,col=4,lwd=2)

# plot the estimated function when considering the positive constraints
posvals = eval.posfd(day,posprecfd1$Wfdobj)
lines(day,posvals,col="black",lwd=2)

# fit data with the smoothing parameter lambda = 10^2
posharmpar = fdPar(fbasis,harmLfd,lambda=1e2)
posprecfd2 = smooth.pos(day,VancPrec,posharmpar)

posvals = eval.posfd(day,posprecfd2$Wfdobj)
lines(day,posvals,col=6,lwd=2)

legend(0,8.5,legend=c('fd lambda=1e4','pos.fd lambda=1e4','pos.fd lambda=1e2'),col=c(4,"black",6),lwd=2)

# Example when estiamting a function with the monotone constraint.

# Now the growth data
names(growth)
growth

hgtm = growth$hgtm # heights of boys
hgtf = growth$hgtf # heights of girls
age = growth$age   # age of measurements


knots  <- growth$age
norder <- 5  # because the transformation is equivalent to estiamting the derivative. 
# So here we increase the order of B-splines by 1. 
nbasis <- length(knots) + norder - 2
rng <- c(min(growth$age),max(growth$age))
hgtbasis <- create.bspline.basis(rng, nbasis, norder, knots)


# define the roughness penalty
Lfdobj <- 2 # use the second derivative to define the roughness penalty
lambda <- 1e-2
growfdPar <- fdPar(hgtbasis, Lfdobj, lambda)

# estimate the function without considering the monotone contraint
# Smoothing for all 39 male subjects
# Each male subject has his data stored in columns of "hgtm"
hgtmfd <- smooth.basis(age, hgtm, growfdPar)$fd
dim(hgtm)
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(age,hgtm[,25],xlab='age',ylab='height',col=2,cex.lab=1.5,cex.axis=1.5)
lines(hgtmfd[25],col=4,lwd=2) # estimated function without considering the monotone contraint

# estimate the function when considering the monotone contraint
monhtmfd25 = smooth.monotone(age,hgtm[,25],growfdPar)
Wfdobj = monhtmfd25$Wfdobj
beta = monhtmfd25$beta

tfine = seq(1,18,len=101)
hgt25vals = eval.monfd(tfine,Wfdobj)
# plot the function when considering the monotone contraint
lines(tfine,beta[1] + beta[2]*hgt25vals,col=6,lwd=2)

# plot the derivative of the estimated function
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(hgtmfd[25],Lfdobj=1,xlab='age',ylab='Dheight',col=4,cex.lab=1.5,cex.axis=1.5,lwd=2)

Dhgt25vals = eval.monfd(tfine,Wfdobj,Lfdobj=1)
lines(tfine,beta[2]*Dhgt25vals,col=6,lwd=2)


# Example to estimate Probability density Function

# simulate some data based on the normal distribution
#  set up range for density
rangeval <- c(-3,3)
#  set up some standard normal data
x <- rnorm(50)
#  make sure values within the range
x[x < -3] <- -2.99
x[x >  3] <-  2.99

# plot the histogram of x
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
hist(x)

#  set up basis for W(x)
# We choose 11 basis functions
basisobj <- create.bspline.basis(rangeval, 11,norder=6)

# define the roughness penalty on W(x)
#  set up initial value for Wfdobj
# W(x) = 0
Wfd0 <- fd(matrix(0,11,1), basisobj)
Wlambda     = 1e-2
WfdParobj <- fdPar(Wfd0, 3, Wlambda)
#  estimate density
denslist <- density.fd(x, WfdParobj)
#  plot density

xval <- seq(-3,3,.2)
wval <- eval.fd(xval, denslist$Wfdobj)
pval <- exp(wval)/denslist$C

quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(xval, pval, type="l", ylim=c(0,1),lwd=3)
points(x,rep(0,50))

