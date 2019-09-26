install.packages("fda")
library(fda)

tobs = seq(0,1,0.01)
nobs = length(tobs)
knots    = c(seq(0,1,0.1));
nknots   = length(knots);
norder   = 4;
nbasis   = length(knots) + norder - 2;
basis = create.bspline.basis(c(min(tobs),max(tobs)),nbasis,norder,knots);


# basis values at samplincurv points

basismat   = eval.basis(tobs, basis);
dim(basismat)
quartz()
plot(tobs,basismat[,1],type = "l",col=1,lwd=3)
lines(tobs,basismat[,2],type = "l",col=2,lwd=3)

quartz()
matplot(tobs,basismat,type='l',lwd=2,lty=1, xlab='day',ylab='basis',cex.lab=1.5,cex.axis=1.5)
for (i in 1:nknots)
{
  abline(v=knots[i],type="l", lty=2, lwd=3)
}

# Comments: If x(0) = 1, set the coefficient to the first basis function = 1; 

# evaluate the first derivative of the basis functions
Dbasismat   = eval.basis(tobs, basis,1);
quartz()
matplot(tobs,Dbasismat,type='l',lwd=2,lty=1, xlab='day',ylab='basis',cex.lab=1.5,cex.axis=1.5)
for (i in 1:nknots)
{
  abline(v=knots[i],type="l", lty=2, lwd=3)
}

# true curve
ytru = (tobs-0.3)^2
plot(tobs,ytru,type = "l")

# put noise to the true curve and generate noisy data
nobs = length(tobs)
noise = 0.03*rnorm(nobs)
yobs = ytru + noise
points(tobs,yobs)

# estimate basis coefficient
Mmat = ginv(t(basismat)%*%basismat)%*%t(basismat)
chat = Mmat%*%yobs

# fitted curve
yhat = basismat%*%chat;
lines(tobs,yhat,type = "l",col="red")

# estimate the variance of noise
SSE = t(yhat-yobs)%*%(yhat-yobs)
sigma2 = SSE/(nobs-nbasis)
sigma2
sqrt(sigma2)
# estimate the variance of the fitted curve
Smat = basismat%*%Mmat
varYhat = diag(Smat%*%Smat*matrix(sigma2,nobs,nobs))

# 95% confidence interval

yhat025 = yhat-1.96*sqrt(varYhat)
yhat975 = yhat+1.96*sqrt(varYhat)

lines(tobs,yhat025,type="l", lty=2, lwd=3,col="blue")
lines(tobs,yhat975,type="l", lty=2, lwd=3,col="blue")

###########################
# Smoothing Splines
#########################


# Use quadrature to get integral - Composite Simpson's Rule

delta <- 0.02
quadpts <- seq(0,1,delta)
nquadpts <- length(quadpts)
quadwts <- as.vector(c(1,rep(c(4,2),(nquadpts-2)/2),4,1),mode="any")
quadwts <- c(1,rep(c(4,2),(nquadpts-1)/2))
quadwts[nquadpts] <- 1
quadwts <- quadwts*delta/3


# Second derivative of basis functions at quadrature points

Q2basismat   = eval.basis(quadpts, basis,2);

# estimates for basis coefficients
Rmat = t(Q2basismat)%*%(Q2basismat*(quadwts%*%t(rep(1,nbasis))))

dim(Rmat)
basismat2 = t(basismat)%*%basismat;
lambda = 0.05   # smoothing parameter
Bmat                      = basismat2 + lambda*Rmat;
chat = ginv(Bmat)%*%t(basismat)%*%yobs;

# fitted value
yhat = basismat%*%chat;
yhat2 = basismat%*%ginv(t(basismat)%*%basismat)%*%t(basismat)%*%yobs;

quartz()
plot(tobs,ytru,type = "l")
points(tobs,yobs)
lines(tobs,yhat,type = "l",col="red")
lines(tobs,yhat2,type = "l",col="blue")

# degrees of freedom
Mmat = ginv(Bmat)%*%t(basismat)
Smat = basismat%*%Mmat
df = sum(diag(Smat))
c(df,nbasis)
# estimate the variance of noise
SSE = t(yhat-yobs)%*%(yhat-yobs)
sigma2 = SSE/(nobs-df)
sigma2
sqrt(sigma2)
# estimate the variance of the fitted curve
varYhat = diag(Smat%*%Smat*matrix(sigma2,nobs,nobs))

# 95% confidence interval

yhat025 = yhat-1.96*sqrt(varYhat)
yhat975 = yhat+1.96*sqrt(varYhat)

lines(tobs,yhat025,type="l", lty=2, lwd=3,col="blue")
lines(tobs,yhat975,type="l", lty=2, lwd=3,col="blue")
















