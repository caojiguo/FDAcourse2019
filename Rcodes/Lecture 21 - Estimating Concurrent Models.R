
## Functional Responses with Functional Predictors:
##              The Concurrent Model
##

install.packages("fda")
library(fda)
#provide directory of the RData
setwd("/Users/cao/Dropbox/Teaching/FDA/SummerCourse2017/R")

# Knee Angle Predicted from Hip Angle

gaittime = seq(0.5,19.5,1)
gaitrange = c(0,20)
gaitfine = seq(0,20,len=101)
length(gaittime)
harmaccelLfd20 = vec2Lfd(c(0, (2*pi/20)^2, 0), rangeval=gaitrange)
gaitbasis = create.fourier.basis(gaitrange, nbasis=21)

gaitLoglam = seq(-4,0,0.25)
nglam   = length(gaitLoglam)

# First select smoothing for the raw data

gaitSmoothStats = array(NA, dim=c(nglam, 3),
      dimnames=list(gaitLoglam, c("log10.lambda", "df", "gcv") ) )
gaitSmoothStats[, 1] = gaitLoglam

#  loop through smoothing parameters

for (ilam in 1:nglam) {
  gaitSmooth = smooth.basisPar(gaittime, gait, gaitbasis,
                   Lfdobj=harmaccelLfd20, lambda=10^gaitLoglam[ilam])
  gaitSmoothStats[ilam, "df"]  = gaitSmooth$df
  gaitSmoothStats[ilam, "gcv"] = sum(gaitSmooth$gcv)
  # note: gcv is a matrix in this case
}

#  display and plot GCV criterion and degrees of freedom

gaitSmoothStats
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(gaitSmoothStats[, c(1, 3)], type='b',cex.lab=2,cex.axis=2)

#  set up plotting arrangements for one and two panel displays
#  allowing for larger fonts
quartz()
par(mfrow=c(1,2),mar = c(8, 8, 4, 2))
plot(gaitSmoothStats[, c(1, 3)], type="b", log="y",lwd=2,cex.lab=2,cex.axis=2)
plot(gaitSmoothStats[, 1:2], type="b", log="y",lwd=2,cex.lab=2,cex.axis=2)

#    GCV is minimized with lambda = 10^(-1.5).

gaitSmooth = smooth.basisPar(gaittime, gait,
       gaitbasis, Lfdobj=harmaccelLfd20, lambda=10^(-1.5))
gaitfd = gaitSmooth$fd

names(gaitfd$fdnames) = c("Normalized time", "Child", "Angle")
gaitfd$fdnames[[3]] = c("Hip", "Knee")
dim(gait)
y2cMap = gaitSmooth$y2cMap
dim(y2cMap)
hipfd  = gaitfd[,1]
kneefd = gaitfd[,2]

kneefdMean = mean(kneefd)

quartz()
par(mfrow=c(3,1),mar = c(8, 8, 4, 2))
plot(kneefdMean, xlab='', ylab='', ylim=c(0, 80),
     main='Mean Knee Angle', lwd=2)
abline(v=c(7.5, 14.7), lty='dashed')
plot(deriv(kneefdMean), xlab='', ylab='',
     main='Knee Angle Velocity', lwd=2)
abline(v=c(7.5, 14.7), h=0, lty='dashed')
plot(deriv(kneefdMean, 2), xlab='', ylab='',
     main='Knee Angle Acceleration', lwd=2)
abline(v=c(7.5, 14.7), h=0, lty='dashed')

quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
phaseplanePlot(gaitfine, kneefdMean,
               labels=list(evalarg=gaittime, labels=1:20),
               xlab='Knee Velocity', ylab='Knee Acceleration')

# Set up a  functional linear regression

xfdlist   = list(const=rep(1,39), hip=hipfd)
betafdPar = fdPar(gaitbasis, harmaccelLfd20)
betalist  = list(const=betafdPar, hip=betafdPar)

gaitRegress= fRegress(kneefd, xfdlist, betalist)

# Intercept

betaestlist = gaitRegress$betaestlist
kneeIntercept = predict(betaestlist$const$fd, gaitfine)

# mean knee angle

kneeMean = predict(kneefdMean, gaitfine)

# Plot intercept & mean knee angle
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
ylim1 = range(kneeIntercept, kneeMean)
plot(gaitfine, kneeIntercept, ylim=ylim1, lwd=2,
     main="Intercept and Mean Knee Angle", type='l',
     xlab='', ylab='')
lines(gaitfine, kneeMean, lty='dashed')
abline(h=0, v=c(7.5, 14.7), lty='dashed')

# Hip coefficient

hipCoef = predict(betaestlist$hip$fd, gaitfine)

kneehatfd = gaitRegress$yhatfd$fd
kneehatmat = eval.fd(gaittime, kneehatfd)
resmat. = gait[,,'Knee Angle'] - kneehatmat
SigmaE = cov(t(resmat.))

# Plot Hip Coefficient 
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
ylim2=c(0, max(hipCoef, knee.R2))
plot(gaitfine, hipCoef, lwd=2, xlab='', ylab='', ylim=ylim2, type='l',
     main='Hip Coefficient')
abline(v=c(7.5, 14.7), lty='dashed')

gaitbasismat = eval.basis(gaitfine, gaitbasis)


fRegressList1 = fRegress(kneefd, xfdlist, betalist,
                         y2cMap=y2cMap, SigmaE=SigmaE)

fRegressList2 = fRegress.stderr(fRegressList1, y2cMap, SigmaE)
betastderrlist = fRegressList2$betastderrlist

quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plotbeta(betaestlist, betastderrlist, gaitfine)

