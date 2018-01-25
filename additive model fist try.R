###12/19/2017 ####
##first try to simulate the additive form model###
#list.of.packages <- c("rgl","ggplot2","knitr","rglwidget")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
##
#if(length(new.packages)) install.packages(new.packages)
##
## By now we have installed the requisite packages. Time to load them .
##

#lapply(list.of.packages,function(x){library(x,character.only=TRUE)})
#library(gam)
#library(mgcv)
#require(akima)
#library(visreg)
rm(list=ls())

library(truncnorm)
library("truncnorm")
library(ggplot2)
library("np")
require(grid)

set.seed(123)
### non-linear function form
### number of obervations####
n=10
beta <- 1
###sigma for the truncnorm distribution###
sigma.v <- 0.5
#### q here is the demension ####
q<-1
c0<-0


###simulation times ###
sims <-1

### choose the functions for use #######
dgp<-1

#### the true value of  m(x)############
re.fun<-function(x,dgp)
{
  if(dgp==1)
  {
    return(x*x)
  }
  else if(dgp==2)
  {
    return(exp(x))
  }

}
###inefficience function form###
scale.fun <- function(dgp,z){

  sf <- if(dgp==1){return(exp(z))}
  sf <- if(dgp==2){return(z^3-0.075*z^2)}

}

###kernel density function ####
ker<-function(x,data,h)
{
  k<-1/(sqrt(2*pi))*exp(-0.5*((x-data)/h)^2)
  #transfer<-(x-data)/h
  #k<-dnorm(transfer,0,1)
  return(k)
}

re<-function(x.it,z.it,y.it)
{
  ##### Bandwidth for the data #######
  h.x<-0.5*sd(x.it)*n^(-1/(4+q))

  h.z<-0.5*sd(z.it)*n^(-1/(4+q))
  ##### matrix store the number of the each observation ####
  k.x<-as.numeric()
  k.z<-as.numeric()
  m.x<-as.numeric()
  m.z<-as.numeric()
  n<-length(x.it)

  for( i in 1:n)
  {
    #### denstiy function for each of the observation in z.it ############
    f.z<-mean(ker(z.it[i],z.it,h.z)/h.z)

    #### density function for each of the observation in x.it ############
    f.x<-mean(ker(x.it[i],x.it,h.x)/h.x)

    #### density function for joint function the specific observation ########
    f.xz<-mean(ker(z.it[i],z.it,h.z)*ker(x.it[i],x.it,h.x))/(h.x*h.z)

    #### the estimation of the density for the y.it based on the x.it
    k.x[i]<-mean(ker(x.it[i],x.it,h.x)*(f.z/f.xz)*(y.it))/h.x

    #### the estimation of the density for the y.it based on the x.it
    k.z[i]<-mean(ker(z.it[i],z.it,h.z)*(f.x/f.xz)*(y.it))/h.z

    #### back-fit by put-back only step #################
    yx.it<-y.it-k.z[i]
    yz.it<-y.it-k.x[i]

    #### update the density function to obtain efficeint estimators
    m.x[i]<-mean(ker(x.it[i],x.it,h.x)*(yx.it))/mean(ker(x.it[i],x.it,h.x))
    m.z[i]<-mean(ker(z.it[i],z.it,h.z)*(yz.it))/mean(ker(z.it[i],z.it,h.z))
  }
  result<-list(m.x,m.z)
  return(result)

}




### the data set information##################################



mmx<-matrix(0,sims,n)
mmz<-matrix(0,sims,n)

MXSE<-matrix(0,sims,1)
MZSE<-matrix(0,sims,1)

### simulation
for( i in 1:sims)
{
  cat(paste(i," of ",sims,sep=""),"\r")
  z.it <- runif(n,-1,1)
  ## Now generate error terms
  v.it <- rnorm(n,sd=sigma.v)
  u <-rtruncnorm(n, a=0, b=Inf, mean = 0, sd = 1)

  ## Scaling function
  u.it <- scale.fun(dgp=dgp,z=z.it)*u
  x.it<-runif(n,-1,1)

  y.it <-c0+re.fun(x=x.it,dgp=dgp)+v.it-u.it

  result<-re(x.it,z.it,y.it)
  mmx[i,]<-result[[1]]
  mmz[i,]<-result[[2]]
  xx<-re.fun(x=x.it,dgp=dgp)
  true.m <- -sqrt(2/pi)*scale.fun(dgp=dgp,z=z.it)
  MXSE[i] <- mean((mmx[i,]-xx)^2)

  MZSE[i] <- mean((mmz[i,]-true.m)^2)





  flush.console()
}
mean(MXSE)
mean(MZSE)
plot(z.it,mmz[i,])
plot(x.it,mmx[i,])

plot3d(y.it, rx[i,], rz[i,], type="s", size=1, lit=TRUE, main = "Additive model form",sub="3-D Plot")

