library(fda)
# setwd("/Users/cao/Dropbox/Teaching/FDA/SummerCourse2018/R")


# First with canadian weather data
# The usual

daybasis365 = create.fourier.basis(c(0,365),365)

# harmonic acceleration differential operator
harmLfd = vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))
tempfdPar = fdPar(daybasis365,harmLfd,1e4)
tempfd = smooth.basis(1:365,daily$tempav,tempfdPar)

quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempfd$fd,xlab='day',ylab='temperature',cex.lab=1.5,cex.axis=1.5)

daily$place
length(daily$place)
# calculate the variance-covariance of the functional data
tempvar = var.fd(tempfd$fd)
tvvals = eval.bifd(1:365,1:365,tempvar)
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
contour(1:365,1:365,tvvals,xlab='day',ylab='day',cex.lab=1.5,cex.axis=1.5)
install.packages("fields")
library(fields)
quartz()
image.plot(1:365,1:365,tvvals,xlab='day',ylab='day',cex.lab=1.5,cex.axis=1.5)

# Correlation Coefficient
temp.cor = cor.fd(1:365,tempfd$fd)
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
contour(1:365,1:365,temp.cor,xlab='day',ylab='day',cex.lab=1.5,cex.axis=1.5)
image.plot(1:365,1:365,temp.cor,xlab='day',ylab='day',cex.lab=1.5,cex.axis=1.5)



# Do functional principal component analysis on the 35 temperature curves
# We choose 4 FPCs

temppca = pca.fd(tempfd$fd,nharm=4)
names(temppca)
temppca$varprop
#temppca$values are the eigenvalues
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(temppca$values[1:8],xlab='component',ylab='variance',col="red",
	cex.lab=1.5,cex.axis=1.5,cex=2)

# plot the cumulative percentage explained total variations
# It shows that the top 3 FPCs explains more than 99% of total variations
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(cumsum(temppca$values[1:10])/sum(temppca$values),xlab='Number of Components',
	ylab='cumulative variance explained',col=2,cex.lab=2,
	cex.axis=2,cex=2)
abline(h=0.99)

# Show the mean curves - temppca$meanfd
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempfd$fd,xlab='day',ylab='temperature',cex.lab=1.5,cex.axis=1.5,col=4)
lines(temppca$meanfd,lwd=2.5,col=2)

# functional principal components
harmfd = temppca$harmonics
harmvals = eval.fd(1:365,harmfd)
dim(harmvals) # The top 4 FPCs

# plot the second FPC
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(1:365,harmvals[,2],xlab='day',ylab='PCs',
     lwd=4,lty=1,cex.lab=2,cex.axis=2,type='l')

# plot all 4 FPCs
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(1:365,harmvals,xlab='day',ylab='PCs',
	lwd=4,lty=1,cex.lab=2.5,cex.axis=2.5,type='l')
legend(0,-0.07,c('PC1','PC2','PC3','PC4'),col=1:4,lty=1,lwd=5)
title('Temperature Principle Component Functions')

# plot the first FPC scores vs. the second FPC scores 
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(temppca$scores[,1:2],xlab='PC Score 1',ylab='PC Score 2',col=4,
	cex.lab=1.5,cex.axis=1.5,cex=1)
text(temppca$scores[,1],temppca$scores[,2],labels=daily$place,cex=1)

quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempfd$fd[17],xlab='day',ylab='temperature',cex.lab=2.5,cex.axis=2.5,col="black",lwd=4,ylim=c(-40,20))
lines(tempfd$fd[24],xlab='day',ylab='temperature',cex.lab=2.5,cex.axis=2.5,col="red",lwd=4)
lines(tempfd$fd[35],xlab='day',ylab='temperature',cex.lab=2.5,cex.axis=2.5,col="blue",lwd=4)
legend(-1, 15, c("Winnipeg", "Calgary","Resolute"),col = c("black","red","blue"),lty=1,lwd=4)

daily$place[c(17,24,35)]
# Remove the mean function
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempfd$fd[17]-temppca$meanfd,xlab='day',ylab='temperature',cex.lab=2.5,cex.axis=2.5,col="black",lwd=4,ylim=c(-40,20))
lines(tempfd$fd[24]-temppca$meanfd,xlab='day',ylab='temperature',cex.lab=2.5,cex.axis=2.5,col="red",lwd=4)
lines(tempfd$fd[35]-temppca$meanfd,xlab='day',ylab='temperature',cex.lab=2.5,cex.axis=2.5,col="blue",lwd=4)
legend(-1, 15, c("Winnipeg", "Calgary","Resolute"),col = c("black","red","blue"),lty=1,lwd=4)


