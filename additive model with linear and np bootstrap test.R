####### 01/20/2018########### second try
rm(list=ls())


library(truncnorm)
library(np)

############ function to define the x############



###### function about the additive, the outcome is the tow additive functions
add.m<-function(y,x,z,h.x,h.z){

  Ybar <- mean(y)

  ## Construct S1 and S2 matrices as in Section 5
  ## Note to self no (1/(nh)) is appearing here
  S1 <- npksum(txdat=x,exdat=x,bws=h.x,return.kernel.weights=TRUE,bandwidth.divide=TRUE)$kw
  S2 <- npksum(txdat=z,exdat=z,bws=h.z,return.kernel.weights=TRUE,bandwidth.divide=TRUE)$kw

  ii<-as.matrix(rep(1,n),n,1)
  y<-as.matrix(y,n,1)

  #### the first stage estimator #########
  ## This follows from page 288 in KLH.
  gamma1.hat<-S1%*%(y*(S2%*%ii)/(n*(S1*S2)%*%ii))
  gamma2.hat<-S2%*%(y*(S1%*%ii)/(n*(S2*S1)%*%ii))

  ## Now calculate oracle/efficient estimator

  ## First neet to construct y^e -- see bottom of page 288 in KLH
  y1.e<-y-gamma2.hat-Ybar
  y2.e<-y-gamma1.hat-Ybar

  ### the efficient estimator ########
  gamma1.hate<-S1%*%y1.e/(S1%*%ii)
  gamma2.hate<-S2%*%y2.e/(S2%*%ii)

  ## Now construct regression function, see top of page 289 in KLH
  mtilde <- gamma1.hat+gamma2.hat-Ybar
  mtilde.e <- gamma1.hate+gamma2.hate-Ybar

  err   <- y-gamma1.hat-gamma2.hat-Ybar
  err.e <- y-gamma1.hate-gamma2.hate-Ybar
  err.c <- err-mean(err)

  ## Return out of function
  residuals<-err.e
  return(list(gamm1.ohat=gamma1.hate,gamm2.ohat=gamma2.hate,residuals=residuals))
}

m.fun <- function(x){return(x)}
g.fun <- function(z,a){return(z+a*exp(z))}

#### scale function that would help check the linear  ###
bootsl<-function(y,x,z,h.x,h.z,boot=100){

  ## Calculate Original Test Statistics
  par.model <-lm(y~x+z)

  ## the residual from the linear model: the null method
  uhat.p<-residuals(par.model)

  ### the new respond variable with error
  yhat <- fitted(par.model)

  sum.sqr.p<-sum(uhat.p^2)

  ### the additive model result, the alternative method

  gamm<-add.m(y,x,z,h.x,h.z)

  uhat.a <-residuals(gamm)
  sum.sqr.a<-sum(uhat.a^2)

  ## J value, the differen between the null and alternative model

  J<-(sum.sqr.p-sum.sqr.a)/sum.sqr.a



  ############# struct the new data sample set###########
  mean.uhat.p<-mean(uhat.p)
  uhat.p.star<-uhat.p-mean.uhat.p

  ## Wild resampling from parametric errors.
  boot.wild <- function(model.resid) {

    p1<-rbinom(n, 1, prob=(1+sqrt(5))/(2*sqrt(5)))
    p2<-1-p1

    ##### the new setting error ###############
    new.re<-(1-sqrt(5))/2*(model.resid)*p1+(1+sqrt(5))/2*(model.resid)*p2

    y.star <- yhat + new.re

    return(y.star)

  }

  ## the matrix pvalue the all the J value for the boostrap sample
  J.boot <- as.numeric()

  for (b in 1:boot){

    ### first got the new y value
    y.star <- boot.wild(uhat.p.star)

    resid.p <- residuals(lm(y.star~x+z))
    sum.sqr.p<-sum(resid.p^2)

    ### the additive regression result
    gamm<-add.m(y.star,x,z,h.x,h.z)

    resid.a <-residuals(gamm)
    sum.sqr.a<-sum(resid.a^2)

    J.boot[b] <- (sum.sqr.p-sum.sqr.a)/sum.sqr.a

  }

  ### calculate how many element exceed the J value from the original sample
  p.value<-length(subset(J.boot,J.boot>J))/boot
  return(list(J=J,J.boot=J.boot,p.value=p.value,boot.number=boot))
}


#### data part #################################
## alpha is the variable control the linear conditioin, the higher, the less linear
alpha<-c(0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9)

## do the bootstrap process trial times,
trial <- 15
pvalue<-matrix(0,trial,length(alpha))


c <- 0.1
g.x<-c*sd(x)*n^(-1/9)
g.z<-c*sd(z)*n^(-1/9)
## Loop over number of times to conduct test
for (j in 1:trial){
  ## Loop over alpha values (size vs. power)

  for(i in 1:length(alpha)){

    n=100
    x <- rnorm(n)
    z <- rnorm(n)
    e <- rnorm(n,sd=0.3)

    m.x <- m.fun(x)
    m.z <- g.fun(z,alpha[i])
    m.x<-m.x-mean(m.x)
    m.z<-m.z-mean(m.z)
    h.x<-0.5*sd(x)*n^(-1/5)
    h.z<-0.5*sd(z)*n^(-1/5)
    y <- m.x+m.z+e
    linear.test<-bootsl(y,x,z,h.x,h.z,boot=100)
    pvalue[j,i]<-linear.test$p.value

  }

}

### calculate the mean value for different alpha
colmean.p=colMeans(pvalue)

##########################################################
### dgp here is to choose different efficient function ###


m.x <- re.fun(x,dgp)
m.z <- scale.fun(dgp,z)
m.x<-m.x-mean(m.x)
m.z<-m.z-mean(m.z)

y <- m.x+m.z+e

gamm<-add.m(y,x,z,h.x,h.z)
gamm1.ohat<-gamm$gamm1.ohat
gamm2.ohat<-gamm$gamm2.ohat

################ plot show the difference between the estimation and the true function
id.x <- order(x)
id.z <- order(z)

par(mfrow=c(1,2))

plot(x[id.x],m.x[id.x],type="l",col="red",lty=1,lwd=2)
lines(x[id.x],gamm1.ohat[id.x],col="blue",lty=2,lwd=2)

plot(z[id.z],m.z[id.z],type="l",col="red",lty=1,lwd=2)
lines(z[id.z],sqrt(2/pi)*gamm2.ohat[id.z],col="blue",lty=2,lwd=2)


######################### the testing part ####################



############# struct the new data sample set#############
############# the bootstrap test to show the nonparametric and parametric difference #########


##########################################################



#### second test for the nonparametric method ##############
boolnp<-function(y,x,z,h.x,h.z,boot=100)
{

  ### null hypothesis: the additive model is right ####
  gamm<-add.m(y,x,z,h.x,h.z)

  uhat.a <-residuals(gamm)
  sum.sqr.a<-sum(uhat.a^2)
  yhat<-gamm$gamm1.ohat+gamm$gamm2.ohat
  mean.uhat.a<-mean(uhat.a)
  uhat.a.star<-uhat.a-mean.uhat.a
  #### alternative hypothesis: the normal nonparametric model ####
  npm<-npreg(y~x+z,bws=c(h.x,h.z),regtype="ll",bandwidth.compute
            =FALSE)

  resid.np<-resid(npm)
  sum.sqr.np<-sum(resid.np^2)

  J<-(sum.sqr.a-sum.sqr.np)/sum.sqr.np
  ##the efficient result of the regression ###


  a <- -0.6180339887499
  P.a <- 0.72360679774998
  b <- 1.6180339887499

  ############# struct the new data sample set###########


  ## Wild resampling from parametric errors.
  boot.wild <- function(model.resid) {
    p1<-rbinom(n, 1, prob=(1+sqrt(5))/(2*sqrt(5)))
    p2<-1-p1


    ##### the new setting error ###############
    new.re<-(1-sqrt(5))/2*(model.resid)*p1+(1+sqrt(5))/2*(model.resid)*p2

    y.star <- yhat + new.re

    return(y.star)

  }


  J.boot <- as.numeric()

  for (b in 1:boot){
    cat(paste(b," of ",boot,sep=""),"\r")
    y.star <- boot.wild(uhat.a.star)


    resid.a <- residuals(add.m(y.star,x,z,h.x,h.z))
    sum.sqr.a<-sum(resid.a^2)
    npm<-npreg(y.star~x+z,bws=c(h.x,h.z),regtype="ll",bandwidth.compute
               =FALSE)
    resid.np<-residuals(npm)
    sum.sqr.np<-sum(resid.np^2)

    J.boot[b] <- (sum.sqr.a-sum.sqr.np)/sum.sqr.np

    flush.console()

  }

  ###### the definition of the p.value here ###########

  p.value<-length(subset(J.boot,J.boot>J))/boot
  return(list(J=J,J.boot=J.boot,p.value=p.value,boot.number=boot))

}
#### the npparametric model test ###########




#### data part #################################
#### the alpha in the function is to decide the separation:
#### the higher alpha, the lesser separation condition

np.fun <- function(x,z,alpha){return(x^2+alpha*x*z+exp(z))}



alpha<-seq(0,5,by=0.4)
trial <- 30
pvalue<-matrix(0,trial,length(alpha))

## Loop over number of times to conduct test
for (j in 1:trial){
  ## Loop over alpha values (size vs. power)
  for(i in 1:length(alpha)){

    n<-200
    x <- rnorm(n)
    z <- rnorm(n)
    e <- rnorm(n,sd=0.3)

    y <- np.fun(x,z,alpha[i])+e
    h.x<-0.5*sd(x)*n^(-1/5)
    h.z<-0.5*sd(z)*n^(-1/5)
    boot<-200
    np.test<-boolnp(y,x,z,h.x,h.z,boot=boot)
    pvalue[j,i]<-np.test$p.value

  }
}
write.table(pvalue,file="pvalue.200.30")
nunp<-apply(pvalue,2,function(pvalue)length(which(pvalue<0.05))/trial)
nunp
## Size Plot w/ Confidence Bounds
a.seq <- seq(0,1,by=0.01)

nominal1 <- as.numeric()
se.nominal1 <- as.numeric()

nominal2 <- as.numeric()
se.nominal2 <- as.numeric()
num.sims=20

for (j in 1:length(a.seq)){

  nominal1[j] <- mean(ifelse(pvalue[,1]<=a.seq[j],1,0))
  se.nominal1[j] <- sqrt(a.seq[j]*(1-a.seq[j])/trial)

  nominal2[j] <- mean(ifelse(pvalue[,1]<=a.seq[j],1,0))
  se.nominal2[j] <- sqrt(a.seq[j]*(1-a.seq[j])/trial)
}

pdf(file="C:/Users/JZhou2/Desktop/modal regression method/additive form model/level_accuracy13.pdf")
plot(a.seq,a.seq,type="l",lwd=2,col="red",lty=1,
     xlab="Nominal Level",ylab="Actual Level",
     main="Level Accuracy")
lines(a.seq,nominal1,lwd=2,col="blue",lty=2)
lines(a.seq,nominal1+2*se.nominal1,lwd=2,col="blue",lty=3)
lines(a.seq,nominal1-2*se.nominal1,lwd=2,col="blue",lty=3)
abline(v=c(0.05,0.25,0.5,0.75,0.95),lwd=1,lty=4,col="black")

lines(a.seq,nominal2,lwd=2,col="black",lty=5)
lines(a.seq,nominal2+2*se.nominal2,lwd=2,col="black",lty=6)
lines(a.seq,nominal2-2*se.nominal2,lwd=2,col="black",lty=6)
dev.off()



colmean.p=colMeans(pvalue)


###### means, for the case that alpha=0, the nomial level, and the accuray level should be equal
### so pvalue[,1]<alpha[j],should be in the sequence follow the alpha

##### the size power ratio means, the speed that approximate to 1 would increase by the increasing of the
### size of the n=sample size.


