
#  ----------------- Climate zone effects for temperature -----------------

#
#  Section 10.1.1 Climate Region Effects on Temperature
#

library(fda)
#  set up the data for the analysis
CanadianWeather$dailyAv[,CanadianWeather$region==regions.[1],1]
regions.         = unique(CanadianWeather$region)
regions.

# find how many cities in each region
Atlantic = CanadianWeather$dailyAv[,CanadianWeather$region==regions.[1],1]
dim(Atlantic)
Continental = CanadianWeather$dailyAv[,CanadianWeather$region==regions.[2],1]
Pacific = CanadianWeather$dailyAv[,CanadianWeather$region==regions.[3],1]
Arctic = CanadianWeather$dailyAv[,CanadianWeather$region==regions.[4],1]

quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(Atlantic,xlab="Day", ylab="Mean Temperature",
        col="red",type="l",ylim=c(-40,25),cex.lab=2,cex.axis=2)
for (j in 1:ncol(Continental))
{lines(Continental[,j],xlab="Day", ylab="Mean Temperature",col="blue",type="l")}
for (j in 1:ncol(Pacific))
{lines(Pacific[,j],xlab="Day", ylab="Mean Temperature",col="black",type="l")}
for (j in 1:ncol(Arctic))
{lines(Arctic[,j],xlab="Day", ylab="Mean Temperature",col="brown",type="l")}
legend('topright', regions., 
       lty=1, col=c('red', 'blue', 'black','brown'), bty='n', cex=.75)


p                = length(regions.) + 1
regionList       = vector("list", p)
names(regionList)= c('Canada', regions.)
regionList[[1]]  = c(rep(1,35),0)
for (j in 2:p) {
  xj             = (CanadianWeather$region == regions.[j-1])
  regionList[[j]]= c(xj,1)
}

# Smoothing on the temperature data

Lcoef        = c(0,(2*pi/365)^2,0)
harmaccelLfd = vec2Lfd(Lcoef, c(0,365))

tempbasis    = create.fourier.basis(c(0, 365), 65)

lambda       = 1e6
tempfdPar65  = fdPar(tempbasis, harmaccelLfd, lambda)

temp  = daily$tempav
day = 1:365
tempSmooth65 = smooth.basis(day, temp, tempfdPar65)
tempfd       = tempSmooth65$fd
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempfd,xlab="Day", ylab="Mean Temperature",lwd=2,cex.lab=2,cex.axis=2)

#  augment tempfd by adding a 36th observation with temp(t) = 0

coef    = tempfd$coef
coef36  = cbind(coef,matrix(0,65,1))
temp36fd= fd(coef36,tempbasis,tempfd$fdnames)

#  set up the regression coefficient list

betabasis      = create.fourier.basis(c(0, 365), 11)
betafdPar      = fdPar(betabasis)
betaList       = vector("list",p)
names(betaList)= regions.
for (j in 1:p) betaList[[j]] = betafdPar

#  carry out the functional analysis of variance

fRegressList= fRegress(temp36fd, regionList, betaList)

#  extract the estimated regression coefficients and y-values

betaestList = fRegressList$betaestlist
regionFit   = fRegressList$yhatfd
regions     = c("Canada", regions.)

# Figure 10.1

quartz()
par(mfrow=c(2,3))
plot(betaestList[[1]]$fd, lwd=2, xlab="Day", ylab="Mean Temperature", main=regions[1],cex.lab=2,cex.axis=2)
for (j in 2:p) plot(betaestList[[j]]$fd, lwd=2, xlab="Day", ylab="Difference", main=regions[j],cex.lab=2,cex.axis=2)

#  ---------------  Testing for no effect of climate zone  ----------------

# temp36fd, regionList, betaList from Section 10.1.1 above

F.res = Fperm.fd(temp36fd, regionList, betaList)
names(F.res)
# Figure 10.14
# plot in black and white
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))

with(F.res,{
  q = 0.95
  ylims = c(min(c(Fvals, qval, qvals.pts)), max(c(Fobs,
                                                  qval)))
  plot(argvals, Fvals, type = "l", ylim = ylims, col = 1,
       lwd = 2, xlab = "day", ylab = "F-statistic",
       cex.lab=1.5,cex.axis=1.5)
  lines(argvals, qvals.pts, lty = 3, col = 1, lwd = 2)
  abline(h = qval, lty = 2, col = 1, lwd = 2)
  legendstr = c("Observed Statistic", paste("pointwise",
                                            1 - q, "critical value"), paste("maximum", 1 -
                                                                              q, "critical value"))
  legend(argvals[1], 1.4, legend = legendstr,
         col = c(1, 1, 1), lty = c(1, 3, 2), lwd = c(2,
                                                     2, 2))
}
)
