############1/15/2018 ##############
### the code with boostrap test 95% confidence ############
rm(list=ls())
n=100
set.seed(123)
x <- runif(n,-1.5,1.5)
z <- runif(n,-1.5,1.5)
e <- rnorm(n,sd=0.1)
y <- sin(x)+sin(z)+e
#### dgp here is to choose different efficient function ###
#dgp=1
#
#m.x <- re.fun(x,dgp)
#m.z <- scale.fun(dgp,z)
#m.x<-m.x-mean(m.x)
#m.z<-m.z-mean(m.z)
a <- -0.6180339887499
P.a <- 0.72360679774998
b <- 1.6180339887499

add.m<-function(y,x,z){


h.x<-0.5**sd(x)*n^(-1/5)
h.z<-0.5**sd(z)*n^(-1/5)
x.k <- npksum(txdat=x,exdat=x,bws=h.x,return.kernel.weights=TRUE,bandwidth.divide=TRUE)$kw
z.k <- npksum(txdat=z,exdat=z,bws=h.z,return.kernel.weights=TRUE,bandwidth.divide=TRUE)$kw
ii<-as.matrix(rep(1,n),n,1)
y<-as.matrix(y,n,1)

#### the normal estimator #########
rr1<-x.k%*%(y*(z.k%*%ii)/(n*(x.k*z.k)%*%ii))
rr2<-z.k%*%(y*(x.k%*%ii)/(n*(z.k*x.k)%*%ii))
y1.e<-y-mean(y)-rr2
y2.e<-y-mean(y)-rr1


### the effiiceint estimator ########
rr1e<-x.k%*%y1.e/(x.k%*%ii)
rr2e<-z.k%*%y2.e/(z.k%*%ii)
mm.add<-rr1e+rr2e+mean(y)
err<-y-mm.add
err.c<-err-mean(err)
return(list(r1=rr1,r2=rr2,r1e=rr1e,r2e=rr2e,err=err,err.e=err.c))
}


boot.wild <- function(model.resid,y) {
  p1<-rbinom(n, 1, prob=(1+sqrt(5))/(2*sqrt(5)))
  p2<-1-p1

  ##### the new setting error ###############
  new.re<-(1-sqrt(5))/2*(model.resid)*p1+(1+sqrt(5))/2*(model.resid)*p2

  y.star <- y + new.re

  return(y.star)

}



addb.m<-function(y,x,z,c=1){
  g.x<-c*sd(x)*n^(-1/9)
  g.z<-c*sd(z)*n^(-1/9)
  x.k <- npksum(txdat=x,exdat=x,bws=g.x,return.kernel.weights=TRUE,bandwidth.divide=TRUE)$kw
  z.k <- npksum(txdat=z,exdat=z,bws=g.z,return.kernel.weights=TRUE,bandwidth.divide=TRUE)$kw
  ii<-as.matrix(rep(1,n),n,1)
  y<-as.matrix(y,n,1)

  #### the normal estimator #########
  rr1<-x.k%*%(y*(z.k%*%ii)/(n*(x.k*z.k)%*%ii))
  rr2<-z.k%*%(y*(x.k%*%ii)/(n*(z.k*x.k)%*%ii))
  y1.e<-y-mean(y)-rr2
  y2.e<-y-mean(y)-rr1


  ### the effiiceint estimator ########
  rr1e<-x.k%*%y1.e/(x.k%*%ii)
  rr2e<-z.k%*%y2.e/(z.k%*%ii)
  mm.add<-rr1e+rr2e+mean(y)
  err<-y-mm.add
  err.c<-err-mean(err)
  return(list(r1=rr1,r2=rr2,r1e=rr1e,r2e=rr2e,err=err,err.e=err.c))
}

###### the result for additive smoother ############
result.o<-add.m(y,x,z)

###### err is the demean error for the whole additive model #######
err<-result.o$err.e

###### the new data for boostrap #########
y.star<-boot.wild(err,y)

### the boostrap number ######
b=200

### the four matrixs to store the boostrap result ###
rr1<-matrix(0,n,b)
rr2<-matrix(0,n,b)
rr1e<-matrix(0,n,b)
rr2e<-matrix(0,n,b)

for(i in 1:b)
{
  y.star<-boot.wild(err,y)
  result.b<-addb.m(y.star,x,z,c=0.1)
  rr1[,i]<-as.matrix(result.b$r1,n,1)
  rr2[,i]<-as.matrix(result.b$r2,n,1)
  rr1e[,i]<-as.matrix(result.b$r1e,n,1)
  rr2e[,i]<-as.matrix(result.b$r2e,n,1)

}









id.x <- order(x)
id.z <- order(z)

yy<-sin(x)
plot( x[id.x], yy[id.x],type="l",col="black" )

lines( x[id.x], result.o$r1e[id.x],type="l",col="blue" )
sse<-result.o$err
sde<-sd(sse)/(sqrt(100))
y1<-result.o$r1e-1.96*sde
y2<-result.o$r1e+1.96*sde
lines( x[id.x], y1[id.x],type="b",col="red" )
lines( x[id.x], y2[id.x],type="b",col="red" )



#################################
rr1b<-matrix(0,n,b)
rr1sd<-matrix(0,n,1)
for(i in 1:length(result.o$r1))
{
  rr1b[i,]<-(rr1e[i,]-mean(rr1e[i,]))/sd(rr1e[i,])
  rr1sd[i]<-sd(rr1e[i,])
}

rr1bl<-apply(rr1b,1,function(x) quantile(x,0.025))
rr1bh<-apply(rr1b,1,function(x) quantile(x,0.975))
rr1l<-result.o$r1e-rr1bl*rr1sd
rr1h<-result.o$r1e-rr1bh*rr1sd
yy<-sin(x)
plot( x[id.x], yy[id.x],type="l",col="black" )

lines( x[id.x], result.o$r1e[id.x],type="l",col="blue" )

lines( x[id.x], rr1l[id.x],type="b",col="grey" )
lines( x[id.x], rr1h[id.x],type="b",col="grey" )

