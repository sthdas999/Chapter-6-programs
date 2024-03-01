
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

## Read the Airfoil Self Noise data from GitHub ##

Arf <- "https://raw.githubusercontent.com/sthdas999/Airfoil-Abalone-datsets/main/Airfoil_Self_Noise_Data.csv"
Airfoil <- read_csv(Arf)
Airfoil
Airfoil1 = head(Airfoil, 100)
n = dim(Airfoil1)[1]

## Renaming the columns ##
freq<- Airfoil1$frequency ## to be taken as X ##
aoa<- Airfoil1$`angle of attack`
cl<- Airfoil1$`chord lengths`
fsv<- Airfoil1$`free-stream velocity`
dt<- Airfoil1$`displacement thickness`
ssp<- Airfoil1$`scaled sound pressure` ## to be taken as Y ##

## frequency vs scaled sound pressure scatterplot ##

plot(x = Airfoil$frequency,y = Airfoil$`scaled sound pressure`,
     xlab = "Frequency",
     ylab = "Scaled sound pressure",
     xlim = c(min(Airfoil$frequency),max(Airfoil$frequency)),
     ylim = c(min(Airfoil$`scaled sound pressure`),max(Airfoil$`scaled sound pressure`)),		 
     main = "Scatterplot of frequency vs scaled sound pressure"
)

x = sort(freq)  ## covariate ##
mx = 0.5*x^2-x^3  ## regression function ##
y = ssp[order(freq)]  ## responses ##
e = y-mx  ## errors ##
h = 0.9 * sd(x) * (n^(-1/5))  ## bandwidth for estimation of regression function ##
m.hat = NadarayaWatsonkernel(x, y, h = 0.1035, gridpoint = x)$mh ## to estimate the unknown regression function using NW method at x.dash in first step ##
m.hat
e.hat = y - m.hat ## estimation of errors ##
e.cen = e.hat - mean(e.hat) ## estimation of centered errors ##
e.cen.boot<- sample(e.cen, n, replace = T) ## resamples of centered errors from the empirical distribution function of centered error ##
y.boot<- m.hat+e.cen.boot ## resampled responses ##
y2.boot<- c()  ## second difference of induced resampled responses ##
for(i in 1:n)
{
  if(i==1)
  {
    y2.boot[i]<- y.boot[2]-y.boot[1]
  }
  else if (i==n)
  {
    y2.boot[i]<- y.boot[n-1]-y.boot[n]
  }
  else
  {
    y2.boot[i]<- y.boot[i+1]-2*y.boot[i]+y.boot[i-1]
  }
}
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
data.xy<- cbind(x,y2.boot,y3.boot) ## combined dataset ##
data.xy
## Generation of datasets on X and 2nd, 3rd order differences for Y through resampling ##
B = 200 ## number of resamples ##
x.dt = vector("list", B)
y2.dt = vector("list", B)
y3.dt = vector("list", B)
for(j in 1:B)
{
  x.dt[[j]] = x
  y2.dt[[j]] = sample(y2.boot, replace = T)
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
          x<-x+(1/4)*((abs(u[i]-u[j])+abs(u[k]-u[l])-abs(u[i]-u[k])-abs(u[j]-u[l]))*(abs(v[i]-v[j])+abs(v[k]-v[l])-abs(v[i]-v[k])-abs(v[j]-v[l])))
        }
      }
    }
  }
  return(x/choose(n,4))
}

#### p-values of the test statistics based on second order differences of Y ####

T1.boots.diff2<- c()
for(j in 1:B)
{
  T1.boots.diff2[j]<- T1(x.dt[[j]], y2.dt[[j]])
}
T1.crit.diff2<- T1(x,y2.boot)
p.value.T1.diff2<- mean(T1.boots.diff2>T1.crit.diff2)
p.value.T1.diff2

T2.boots.diff2<- c()
for(j in 1:B)
{
  T2.boots.diff2[j]<- T2(x.dt[[j]], y2.dt[[j]])
}
T2.crit.diff2<- T2(x,y2.boot)
p.value.T2.diff2<- mean(T2.boots.diff2>T2.crit.diff2)
p.value.T2.diff2

T3.boots.diff2<- c()
for(j in 1:B)
{
  T3.boots.diff2[j]<- T3(x.dt[[j]], y2.dt[[j]])
}
T3.crit.diff2<- T3(x.dt[[j]], y2.dt[[j]])
p.value.T3.diff2<- mean(T3.boots.diff2>T3.crit.diff2)
p.value.T3.diff2

#### p-values of the test statistics based on third order differences of Y ####

T1.boots.diff3<- c()
for(j in 1:B)
{
  T1.boots.diff3[j]<- T1(x.dt[[j]], y3.dt[[j]])
}
T1.crit.diff3<- T1(x,y3.boot)
p.value.T1.diff3<- mean(T1.boots.diff3>T1.crit.diff3)
p.value.T1.diff3

T2.boots.diff3<- c()
for(j in 1:B)
{
  T2.boots.diff3[j]<- T2(x.dt[[j]], y3.dt[[j]])
}
T2.crit.diff3<- T2(x,y3.boot)
p.value.T2.diff3<- mean(T2.boots.diff3>T2.crit.diff3)
p.value.T2.diff3

T3.boots.diff3<- c()
for(j in 1:B)
{
  T3.boots.diff3[j]<- T3(x.dt[[j]], y3.dt[[j]])
}
T3.crit.diff3<- T3(x.dt[[j]], y3.dt[[j]])
p.value.T3.diff3<- mean(T3.boots.diff3>T3.crit.diff3)
p.value.T3.diff3

p.values.diff2 = cbind(p.value.T1.diff2,p.value.T2.diff2,p.value.T3.diff2)
p.values.diff2
p.values.diff3 = cbind(p.value.T1.diff3,p.value.T2.diff3,p.value.T3.diff3)
p.values.diff3









