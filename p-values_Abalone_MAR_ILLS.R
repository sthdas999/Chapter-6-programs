
## R version 4.0.5 ##

## Required libraries ##
library(MASS)
library(kedd)
library(mvtnorm)
library(Matrix)
library(truncnorm)
library(npmlda)
library(stats)
library(EnvStats)
library(lattice)
library(pracma)
library(readr)  # for read_csv
library(knitr)  # for kable
library(bbemkr) ## NadarayaWatsonkernel() ##

## Read the Abalone data from GitHub ##

Abl <- "https://raw.githubusercontent.com/sthdas999/Airfoil-Abalone-datsets/main/abalone.data.csv"
Abalone <- read_csv(Abl)
D = dim(Abalone)[1]
Abalone1 = head(Abalone, 100)
n = dim(Abalone1)[1]
abalone1<- Abalone1[-c(1, 9)]

## Renaming the columns ##
Len<- abalone1$Length
Diam<- abalone1$Diameter
Ht<- abalone1$Height  ## parametric covariate Z ##
Wt<- abalone1$`Whole weight`
Shuc.wt<- abalone1$`Shucked weight`
Visc.wt<- abalone1$`Viscera weight`
Sh.wt<- abalone1$`Shell weight`   ## nonparametric covariate X ##
Age<- abalone1$Age   ## response ##

## Scatterplot of Age vs Height ##

plot(x = abalone1$Height,y = abalone1$Age,
     xlab = "Height",
     ylab = "Age",
     xlim = c(min(Ht),max(Ht)),
     ylim = c(min(Age),max(Age)),		 
     main = "Scatterplot of Age vs Height"
)

## Scatterplot of Age vs Shell weight ##

plot(x = abalone1$`Shell weight`,y = abalone1$Age,
     xlab = "Height",
     ylab = "Age",
     xlim = c(min(Ht),max(Ht)),
     ylim = c(min(Age),max(Age)),		 
     main = "Scatterplot of Age vs Shell weight"
)

x = Sh.wt 
z = Ht
y = Age
h = 0.9 * sd(x) * (n^(-1/5))  ## bandwidth for estimation of regression function ##
## estimation of beta in Y=Z*beta+m(X)+epsilon ##
k<- function(u)  ## definition of epanechnikov kernel ##
{
  return((3/4)*(1-u^2)*(abs(u)<=1))
}
gy.hat<- function(t)  ## estimation of gy.hat = regression function of Y on X ##
{
  u<- (x-t)/h
  m<- k(u)*y
  if(sum(k(u)!=0))
  {
    return((sum(m)/sum(k(u))))
  }
  else
  {
    return(0)
  }
}
gx.hat<- function(t)  ## estimation of gy.hat = regression function of Y on Z ##
{
  u<- (x-t)/h
  m<- k(u)*z
  if(sum(k(u)!=0))
  {
    return((sum(m)/sum(k(u))))
  }
  else
  {
    return(0)
  }
}
ex.hat<- c()  ## estimation of error in the regression function of Y on Z ##
ey.hat<- c()  ## estimation of error in the regression function of Y on X ##
for(i in 1:n)
{
  ex.hat[i] = z[i] - gx.hat(x[i])
  ey.hat[i] = y[i] - gy.hat(x[i])
}
beta.hat = sum(ex.hat*ey.hat) / sum(ex.hat^2)  ## estimation of beta
beta.hat
y = y-z*beta.hat ## transformed non-missing response ##
p.values.MAR.ILLS = function(p) ## p being the proportion of missingness of responses ##
{
  n.hat = floor(n*p) ## number of observations in Y to make missing ##
  ## Now, n.hat number of observations on Y are to be made missing through MAR technique ##
  
  m.hat <- function(x, y, gridpoint)  ## Definition of NW estimator based on Epanechnikov kernel ##
  {
    ker = function(u) 0.75*(1-u^2) ## kernel = Epanechnikov ##
    n = length(y)
    mh = vector(,length(gridpoint))
    for(j in 1:length(gridpoint))
    {
      suma = sumb = vector(,n)
      for(i in 1:n)
      {
        suma[i] = ker((gridpoint[j] - x[i])/h) * y[i]
        sumb[i] = ker((gridpoint[j] - x[i])/h)
      }
      mh[j] = sum(suma)/sum(sumb)
    }
    return(list(gridpoint = gridpoint, mh = mh))
  }
  
  m.hat.x = function(x) m.hat(x,y,x)$mh 
  
  p.hat = function(x) exp(m.hat.x(x))*(1+exp(m.hat.x(x)))^-1 ## probabilities of missingness of X values as logit function ##
  
  round(p.hat(x),3) -> phat
  
  prob = sample(phat,n.hat,replace = F)  ## generation of n.hat number of probabilities ##
  
  count.1 = c() ## missing positions corresponding to the generated probabilities ##
  for(i in 1:n.hat)
  {
    for(j in 1:n)
    {
      if(prob[i]==phat[j])
      {
        count.1[i] = j
      }
    }
  }
  ## count.1 values need to be distinct ##
  y.miss = c()  ## the Y-values at the count.1 digited places are defined as NA's ##
  for(i in 1:n)
  {
    if(i %in% count.1)
    {
      y.miss[i] = NA
    }
    else
    {
      y.miss[i] = y[i]
    }
  }
  
  y.dash = y.miss[-c(count.1)]  ## Y-observations after removal of NA values from y.incomplete ##
  
  x.dash = x[-c(count.1)]  ## X-observations corresponding to y.dash ##
  
  ## Now, to estimate the unknown regression function using ILLS method at x.dash in first step as follows. ##
  
  ep.kernel = function(u) 0.75*(1-u^2) ## Epanechnikov kernel ##
  
  M1<- function(u)  
  {
    vec1<- c()
    for(i in 1:length(x.dash))
    {
      vec1[i]<- (u-x.dash[i])*ep.kernel((u-x.dash[i])/h)
    }
    return(sum(vec1))
  }
  
  M2<- function(v)  
  {
    vec2<- c()
    for(i in 1:length(x.dash))
    {
      vec2[i]<- (v-x.dash[i])^2*ep.kernel((v-x.dash[i])/h)
    }
    return(sum(vec2))
  }
  
  lls1<- function(c)
  {
    g1 = c()
    for(i in 1:length(x.dash))  { g1[i]<- (M2(x.dash[i])-(x.dash[i]-c)*M1(x.dash[i]))*ep.kernel((x.dash[i]-c)/h)*y.dash[i]  }
    g2 = c()
    for(i in 1:length(x.dash))  { g2[i]<- (M2(x.dash[i])-(x.dash[i]-c)*M1(x.dash[i]))*ep.kernel((x.dash[i]-c)/h)  }
    return(sum(g1)/sum(g2))
  }
  
  ## least square estimation of intercept based on non-missing (X,Y)-observations only ##
  
  t0.hat<- c()  ## t0.hat is also the first step estimator of the regression function based on non-missing observations of (X,Y) ##
  for(i in 1:length(x.dash))
  {
    t0.hat[i]<- lls1(x.dash[i])
  }
  x.miss = x[count.1]  ## X-observations corresponding to missing Y-values ##
  m.hat.miss = c()  ## estimation of regression function at the missing observations of Y ##
  for(i in 1:length(x.miss))
  {
    m.hat.miss[i] = lls1(x.miss[i])
  }
  m.hat.miss.arranged<- c()  
  for(i in 1:length(count.1))
  {
    m.hat.miss.arranged[i]<- m.hat.miss[rank(count.1)==i]
  }
  y.complete<- replace(y.miss,which(is.na(y.miss)==T),m.hat.miss.arranged)   ## complete data on Y after imputation ##
  ## Now, based on n paired observations (X, Y.complete), the regression function of the complete model needs to be estimated. ##
  M1.1<- function(u)  
  {
    v1<- c()
    for(i in 1:n)
    {
      v1[i]<- (u-x[i])*ep.kernel((u-x[i])/h)
    }
    return(sum(v1))
  }
  M2.1<- function(v)  
  {
    v2<- c()
    for(i in 1:n)
    {
      v2[i]<- (v-x[i])*ep.kernel((v-x[i])/h)
    }
    return(sum(v2))
  }
  lls2<- function(c)
  {
    r1 = c()
    for(i in 1:n)  { r1[i]<- (M2.1(x[i])-(x[i]-c)*M1.1(x[i]))*ep.kernel((x[i]-c)/h)*y.complete[i]  }
    r2 = c()
    for(i in 1:n)  { r2[i]<- (M2.1(x[i])-(x[i]-c)*M1.1(x[i]))*ep.kernel((x[i]-c)/h)  }
    return(sum(r1)/sum(r2))
  }
  ## least square estimation of intercept based on full set of (X,Y)-observations ##
  s0.hat<- c()  ## t0.hat is also the first step estimator of the regression function based on non-missing observations of (X,Y) ##
  for(i in 1:n)
  {
    s0.hat[i]<- lls2(x[i])
  }
  ## s0.hat is the secind step cum final estimator of the regression function of the model ##
  e.hat = y.complete - s0.hat ## estimation of errors ##
  e.cen = e.hat - mean(e.hat) ## estimation of centered errors ##
  e.cen.boot<- sample(e.cen, n, replace = T) ## resamples of centered errors from the empirical distribution function of centered error ##
  y.boot<- s0.hat+e.cen.boot ## resampled responses ##
  y3.boot<- c()  ## third difference of induced resampled responses ##
  for(i in 1:n)
  {
    if(i==1)
    {
      y3.boot[i] = -y.boot[2]
    }
    else if(i==2)
    {
      y3.boot[i] = -2*y.boot[1]+3*y.boot[2]-y.boot[3]
    }
    else if(i==n)
    {
      y3.boot[i] = y.boot[n-2]-3*y.boot[n-1]+2*y.boot[n]
    }
    else
    {
      y3.boot[i] = y.boot[i-2]-3*y.boot[i-1]+3*y.boot[i]-y.boot[i+1]
    }
  }
  ## Generation of datasets on X and 2nd, 3rd order differences for Y through resampling ##
  B = 200 ## number of resamples ##
  x.dt = vector("list", B)
  y3.dt = vector("list", B)
  for(j in 1:B)
  {
    x.dt[[j]] = x
    y3.dt[[j]] = sample(y3.boot, replace = T)
  }
  ################################ Test Statistics #####################################
  ####### T1 #########
  T1 = function(u, v)
  {
    temp = 0
    n<- length(u)
    for(i in 1:(n-1))
    {
      for(j in (i+1):n)
      {
        temp = temp + sign((u[i]-u[j])*(v[i]-v[j]))
      }
    }
    return(temp/choose(n, 2))
  }
  ####### T2 #########
  T2<-function(u,v)
  {
    x<-0
    for(i in 1:(n-3))
    {
      for(j in (i+1):(n-2))
      {
        for(k in (j+1):(n-1))
        {
          for(l in (k+1):n)
          {
            x<-x+(sign(abs(u[i]-u[j])+abs(u[k]-u[l])-abs(u[i]-u[k])-abs(u[j]-u[l]))*sign(abs(v[i]-v[j])+abs(v[k]-v[l])-abs(v[i]-v[k])-abs(v[j]-v[l])))
          }
        }
      }
    }
    return(x/choose(n,4))
  }
  ####### T3 #########
  T3<-function(u,v)
  {
    x<-0
    for(i in 1:(n-3))
    {
      for(j in (i+1):(n-2))
      {
        for(k in (j+1):(n-1))
        {
          for(l in (k+1):n)
          {
            x<-x+((abs(u[i]-u[j])+abs(u[k]-u[l])-abs(u[i]-u[k])-abs(u[j]-u[l]))*(abs(v[i]-v[j])+abs(v[k]-v[l])-abs(v[i]-v[k])-abs(v[j]-v[l])))
          }
        }
      }
    }
    return(x/choose(n,4))
  }  
  T1.boots.diff3<- c()
  for(j in 1:B)
  {
    T1.boots.diff3[j]<- T1(x.dt[[j]], y3.dt[[j]])
  }
  T1.crit.diff3<- T1(x,y3.boot)
  p.value.T1.diff3<- mean(T1.boots.diff3>T1.crit.diff3)
  T2.boots.diff3<- c()
  for(j in 1:B)
  {
    T2.boots.diff3[j]<- T2(x.dt[[j]], y3.dt[[j]])
  }
  T2.crit.diff3<- T2(x,y3.boot)
  p.value.T2.diff3<- mean(T2.boots.diff3>T2.crit.diff3)
  T3.boots.diff3<- c()
  for(j in 1:B)
  {
    T3.boots.diff3[j]<- T3(x.dt[[j]], y3.dt[[j]])
  }
  T3.crit.diff3<- T3(x.dt[[j]], y3.dt[[j]])
  p.value.T3.diff3<- mean(T3.boots.diff3>T3.crit.diff3)
  p.values.diff3 = cbind(p.value.T1.diff3,p.value.T2.diff3,p.value.T3.diff3)
  return(p.values.diff3)
}

