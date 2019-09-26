## A Functional Linear Model for Swedish Mortality
##


#  SwedeMat:  a dataframe object with 81 rows and 144 columns
#             containing the log hazard values for ages 0 through 80
#             and years 1751 through 1884
# The hazard rate refers to the rate of death for an item of a given age



library(fda)
#provide directory of the RData
setwd("/Users/cao/Dropbox/Teaching/FDA/SummerCourse2018/R")
mat=load("Sweden.Rdata")
mat=as.matrix(SwedeMat)
quartz()
matplot(mat,type="l")
SwedeMat = mat
dim(SwedeMat)

SwedeLogHazard = as.matrix(SwedeMat)

  dimnames(SwedeLogHazard)[[2]] <- paste('b', 1751:1894, sep='')

  Fig10.10data = cbind(SwedeLogHazard[, c('b1751', 'b1810', 'b1860')],
                       Swede1920)

  SwedeTime = 0:80;
  SwedeRng = c(0,80);
  quartz()
  matplot(SwedeTime, Fig10.10data,
          type='l',lwd=2,xlab='age',ylab='log Hazard',col=1,
          cex.lab=1.5,cex.axis=1.5)

#  smooth the log hazard observations

  nbasis = 85
  norder = 6
  SwedeBasis = create.bspline.basis(SwedeRng, nbasis, norder)

  D2fdPar = fdPar(SwedeBasis, lambda=1e-7)

  SwedeLogHazfd = smooth.basis(SwedeTime, SwedeLogHazard, D2fdPar)$fd

# The following requires manually clicking on the plot
# for each of 144 birth year cohorts

  plotfit.fd(SwedeLogHazard,SwedeTime,SwedeLogHazfd)

# Set up for the list of regression coefficient fdPar objects

  nbasis     = 23
  SwedeRng   = c(0,80)
  SwedeBetaBasis = create.bspline.basis(SwedeRng,nbasis)

  SwedeBeta0Par = fdPar(SwedeBetaBasis, 2, 1e-5)

  SwedeBeta1fd  = bifd(matrix(0,23,23), SwedeBetaBasis, SwedeBetaBasis)

  SwedeBeta1Par = bifdPar(SwedeBeta1fd, 2, 2, 1e3, 1e3)

  SwedeBetaList = list(SwedeBeta0Par, SwedeBeta1Par)

#  Define the dependent and independent variable objects

  NextYear = SwedeLogHazfd[2:144]
  LastYear = SwedeLogHazfd[1:143]

#  Do the regression analysis

  Swede.linmod = linmod(NextYear, LastYear, SwedeBetaList)

  Swede.ages = seq(0, 80, 2)
  Swede.beta1mat = eval.bifd(Swede.ages, Swede.ages, Swede.linmod$beta1estbifd)

  quartz()
  persp(Swede.ages, Swede.ages, Swede.beta1mat,
        xlab="age", ylab="age",zlab="beta(s,t)",
        cex.lab=1.5,cex.axis=1.5)
