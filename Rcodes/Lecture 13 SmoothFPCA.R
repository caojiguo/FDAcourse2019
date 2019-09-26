library(fda)

# Smoothed FPCA

# First with canadian weather data
# The usual
daybasis365 = create.fourier.basis(c(0,365),365)

harmLfd = vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))

tempfdPar = fdPar(daybasis365,harmLfd,1e4)
tempfd = smooth.basis(1:365,daily$tempav,tempfdPar)

quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempfd$fd,xlab='day',ylab='temperature',cex.lab=1.5,cex.axis=1.5,col=4)

# I tried differential values of the smoothing parameter

tempfdPar2 = fdPar(daybasis365,harmLfd,1e-2)
tempfd2 = smooth.basis(1:365,daily$tempav,tempfdPar2)

quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempfd2$fd,xlab='day',ylab='temperature',cex.lab=1.5,cex.axis=1.5,col=4)


tempfdPar3 = fdPar(daybasis365,harmLfd,1e6)
tempfd3 = smooth.basis(1:365,daily$tempav,tempfdPar3)

quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempfd3$fd,xlab='day',ylab='temperature',cex.lab=1.5,cex.axis=1.5,col=4)


# do FPCA with a roughness penalty on FPCs.
# We use the curves esitmated with the smoothing spline 
# with the smoothing parameter lambda = 10^{-2}

# Here we define the roughness penalty by the harmonic acceleration differential operator 
# with the smoothing parameter lambda = 10^6
ptemppca = pca.fd(tempfd2$fd,nharm=4,harmfdPar=tempfdPar3)

# Get FPCs
pharmfd = ptemppca$harmonics
pharmvals = eval.fd(1:365,pharmfd)

quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(1:365,pharmvals,xlab='day',ylab='PCs',
        lwd=2,lty=1,cex.lab=1.5,cex.axis=1.5,type='l')
legend(0,-0.05,c('PC1','PC2','PC3','PC4'),col=1:4,lty=1,lwd=2)
title('Temperature Principle Component Functions')



# Now let's do a bit of reconstruction for Montreal

# Do FPCA on the temperature curves without the Montreal curve
Stemppca = pca.fd(tempfd2$fd[-12],nharm=4,harmfdPar=tempfdPar3)
# Get the FPC
harms = Stemppca$harmonics
# Get the mean curve
meanfd = Stemppca$meanfd

# Montreal data
Mdat = CanadianWeather$dailyAv[,'Montreal','Temperature.C']
length(Mdat)
# evaluate FPC in the days [1:132]
Stempvals = eval.fd(day.5[1:132],harms)
# evaluate the mean curve in the days [1:132]
mtempvals = eval.fd(day.5[1:132],meanfd)

# Remove Mean curve from the Montreal data
Mdat2 = Mdat[1:132]-mtempvals

# Obtain the FPC scores of the Montreal curve
coef = lm(Mdat2~Stempvals-1)$coef
coef

# Prediction for the Montreal curve
Rfd = coef[1]*harms[1]+coef[2]*harms[2]+
  coef[3]*harms[3]+coef[4]*harms[4]+
  Stemppca$meanfd

Rvals = eval.fd(day.5,Rfd)
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(day.5,Rvals,type='l',lwd=10,col=4,xlab='day',ylab='Temperature',
     cex.lab=2.5,cex.axis=2.5)
points(day.5,Mdat,col=2,lwd=2)
points(day.5[1:132],Mdat[1:132],col=3)
lines(Stemppca$meanfd,lty=2,col=5,lwd=2)

# Multivariate FPCA

library(fda)
?gait
# Looking at the gait data

gaittime = seq(0,1,len=20)

quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(gaittime,gait[,,1],type='l',cex.lab=1.5,cex.axis=1.5,ylab='Hip')
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(gaittime,gait[,,2],type='l',cex.lab=1.5,cex.axis=1.5,ylab='Knee')

quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(gait[,,1],gait[,,2],type='l',cex.lab=1.5,cex.axis=1.5,xlab='hip',ylab='Knee')

# Setting up a general object

# define the harmonic acceleration differential operator
harmaccelLfd <- vec2Lfd(c(0, 0, (2*pi)^2, 0))
gaitbasis21 <- create.fourier.basis(nbasis=21)

# Smooth the gait data
gaitfd <- smooth.basisPar(gaittime, gait,
       gaitbasis21, Lfdobj=harmaccelLfd, lambda=1e-11)$fd

names(gaitfd$fdnames) = c("Normalized time", "Child", "Angle")
gaitfd$fdnames[[3]] = c("Hip", "Knee")

# Plot the smooth

tfine = seq(0,1,0.01)
gaitvals = eval.fd(tfine,gaitfd)
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(gaitvals[,,1],gaitvals[,,2],type='l',cex.lab=1.5,cex.axis=1.5,xlab='hip',ylab='Knee')

gait.cor = cor.fd(tfine,gaitfd)
dim(gait.cor)
library("fields")
quartz()
par(mfrow=c(1,3),mar = c(8, 8, 4, 2))
#  contour(tfine,tfine,gait.cor[,,1,3],xlab='day',ylab='day',cex.lab=1.5,cex.axis=1.5)
  image.plot(tfine,tfine,gait.cor[,,1,1],gait.cor,
             xlab='hip',ylab='hip',cex.lab=1.5,cex.axis=1.5)
  image.plot(tfine,tfine,gait.cor[,,1,2],gait.cor,
             xlab='hip',ylab='Knee',cex.lab=1.5,cex.axis=1.5)
  image.plot(tfine,tfine,gait.cor[,,1,3],gait.cor,
             xlab='Knee',ylab='Knee',cex.lab=1.5,cex.axis=1.5)      



# Now a principle components analysis

gait.pca = pca.fd(gaitfd,nharm=4)

quartz()
par(mfrow=c(2,1),mar = c(8, 8, 4, 2))
plot(gait.pca$meanfd,cex.lab=1.5,cex.axis=1.5)

# Mean cycle

meanvals = eval.fd(tfine,gait.pca$meanfd)
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(meanvals[,,1],meanvals[,,2],xlab='hip',ylab='knee',
	cex.lab=1.5,cex.axis=1.5,type='l',lwd=2,col=4)


# Plot the variance proportion
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(gait.pca$values,cex.lab=1.5,cex.axis=1.5,xlab='PC',ylab='eigenvalue',col=4,cex=2)
gait.pca$varprop
sum(gait.pca$varprop)


# plot FPCs

quartz()
par(mfrow=c(1,2),mar = c(8, 8, 4, 2))
plot(gait.pca$harmonics,lwd=2,cex.lab=1.5,cex.axis=1.5)

