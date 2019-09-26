install.packages("fda")
library(fda)
#provide directory of the RData
# setwd("/Users/cao/Dropbox/Teaching/FDA/SummerCourse2017/R")

# Use the scalar-on-function functional linear model to find the effect of 
# the daily tempreture curve on the annual precipitations.
# x_i(t): the daily tempreture curve for the i-th city
# y_i   : the annual precipitation for the i-th city;

# obtain the annual precipitation for 35 cities
annualprec   = log10(apply(daily$precav,2,sum))
ny = length(annualprec)
(ny = length(annualprec))
# Define the 65 Fourier basis functions
tempbasis65  = create.fourier.basis(c(0,365),65)

# Smooth the daily temperature data to obtain a smooth temperature curve
tempSmooth65 = smooth.basis(day.5, daily$tempav, tempbasis65)
tempfd65     = tempSmooth65$fd


templist      = vector("list",2)
templist[[1]] = rep(1,35) # the first functional covariate
templist[[2]] = tempfd65  # the second functional covariate



# create a constant basis for the intercept
conbasis   = create.constant.basis(c(0,365))
quartz()
plot(conbasis)
# Define the small number of basis functions for beta(t)
# We do not use roughness penalty here
nbasis = 11
betabasis5 = create.fourier.basis(c(0,365),nbasis)
betalist1  = vector("list",2)
betalist1[[1]] = conbasis
betalist1[[2]] = betabasis5

# fit the functional linear model 
fRegressList1 = fRegress(annualprec,templist,betalist1)
names(fRegressList1)
betaestlist1  = fRegressList1$betaestlist
length(betaestlist1)
# betaestlist1 has two elements. The first element is the intercept
# The second element is the slope beta(t)

# obtain beta(t)
tempbetafd1   = betaestlist1[[2]]$fd

quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempbetafd1, xlab="Day", ylab="Beta for temperature")

# obtain the intercept alpha
coef(betaestlist1[[1]])

# alpha is not equal to the mean of annualprec
mean(annualprec)

# fitted value yhat
annualprechat1 = fRegressList1$yhatfdobj
# fitted residual
annualprecres1 = annualprec - annualprechat1
# sum squared residuals
SSE1  = sum(annualprecres1^2)

# sum squared residuals for the null model y = alpha + \epsilon
SSE0    = sum((annualprec - mean(annualprec))^2)

# F test for the overall effect of x_i(t)
# H0: y = alpha + \epsilon
# H1: y = alpha + \int [beta(t)x(t)]dt + epsilon
(Fratio = ((SSE0-SSE1)/(nbasis-1)/(SSE1/(ny-nbasis))))

# 95% quantile of F(11,23)
qf(0.95,nbasis,ny-nbasis-1)

# Fratio >>qf(0.95,nbasis-1,ny-nbasis)
# indicating that x_i(t) has a significant effect on y_i

# calculate the p-value
1-pf(Fratio,nbasis-1,ny-nbasis)
# p-value is 1.7*10^(-8) indicating that x_i(t) has a significant effect on y_i



# Penalized Estimation

#Using the harmonic acceleration differential operator 
# to define roughness penalty on beta(t)
Lcoef = c(0,(2*pi/365)^2,0)
harmaccelLfd = vec2Lfd(Lcoef, c(0,365))

# We use 35 Fourier basis functions to represent beta(t)
betabasis35 = create.fourier.basis(c(0, 365), 35)

# Choosing Smoothing Parameters using cross-validation

loglam = seq(5,15,0.5)
nlam   = length(loglam)
SSE.CV = rep(NA,nlam)
for (ilam in 1:nlam) {
  print(paste("log lambda =", loglam[ilam]))
  lambda     = 10^(loglam[ilam])
  betalisti  = betalist1
  betalisti[[2]] = fdPar(betabasis35, harmaccelLfd, lambda)
  fRegi          = fRegress.CV(annualprec, templist, betalisti)
  SSE.CV[ilam]   = fRegi$SSE.CV
}
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(loglam, SSE.CV, type="b", lwd=2,
     xlab="log smoothing parameter lambda",
     ylab="Cross-validation score", cex.lab=2,cex.axis=2)

# Choose lambda which minimize SSE.CV
lambda      = 10^12.5
betafdPar.  = fdPar(betabasis35, harmaccelLfd, lambda)

betalist2      = betalist1
betalist2[[2]] = betafdPar.
# the first element of betalist2 is still the constant for the intercept

# do functional linear model
annPrecTemp    = fRegress(annualprec, templist, betalist2)

# get beta(t)
betaestlist2   = annPrecTemp$betaestlist

# get fitted value yhat
annualprechat2 = annPrecTemp$yhatfdobj

# get the effective degrees of freedom
print(annPrecTemp$df)

# do the F test 
# test the overall effect of x_i(t) on y_i
(SSE2 = sum((annualprec-annualprechat2)^2))
(Fratio2 = ((SSE0-SSE2)/(annPrecTemp$df-1))/(SSE2/(ny-annPrecTemp$df)))
c(Fratio,Fratio2)
# Fratio2 > Fratio

# 95% quantile 
qf(0.95,annPrecTemp$df-1,ny-annPrecTemp$df)

# p-value
1-pf(Fratio2,annPrecTemp$df-1,ny-annPrecTemp$df)

# plot beta(t)
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(betaestlist2[[2]]$fd, xlab="Days",ylab="beta(t)",lwd=2,cex.lab=2,cex.axis=2)

# plot yhat vs. y
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(annualprechat2, annualprec, lwd=2,cex.lab=2,cex.axis=2)
abline(lm(annualprec~annualprechat2), lty='dashed', lwd=2)
abline(0,1,lty=1, lwd=2,col="red")





# Confidence Intervals

# fitted residuals
resid   = annualprec - annualprechat2

# estimate sigma^2
SigmaE. = sum(resid^2)/(35-annPrecTemp$df)
SigmaE  = SigmaE.*diag(rep(1,35))

# for smoothing temperature chat = a matrix * y
y2cMap  = tempSmooth65$y2cMap

# obtain point-wise standard error for beta(t)
stderrList = fRegress.stderr(annPrecTemp, y2cMap, SigmaE)

betafdPar      = betaestlist2[[2]]
betafd         = betafdPar$fd
betastderrList = stderrList$betastderrlist
betastderrfd   = betastderrList[[2]]

quartz()
plot(betafd, xlab="Day", ylab="Temperature Reg. Coeff.",
     ylim=c(-6e-4,1.2e-03), lwd=2,cex.lab=2,cex.axis=2)
lines(betafd+2*betastderrfd, lty=2, lwd=2)
lines(betafd-2*betastderrfd, lty=2, lwd=2)
# The temperature has a significant positive effect on the annual precipitations from September to December. 
# There is no significant effect of the temperature on the annual precipitations in other time.
