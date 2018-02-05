#### additive with restricted and unrestricted conditional mean###

rm(list=ls())
install.packages('stats')
library(stats)
library(lpSolve)
library(lpSolveAPI)
library(quadprog)
library(Deriv)
library(np)
n=100
x<-runif(n,-2,2)
u<-rnorm(n,0,1)
y<-x^3+u
h.x<-0.5*sd(x)*n^{-1/3}

### This is the function to find the kernel value for each element=
mx<-npksum(txdat=x,exdat=x,bws=h.x,return.kernel.weights=TRUE,bandwidth.divide=TRUE)$kw

### ii is to calculate the sum value of the function
ii<-as.matrix(rep(1,n),n,1)
### this is the function to calculate the first order derivatives
DD <- function(expr, name, order = 1) {
  if(order < 1) stop("'order' must be >= 1")
  if(order == 1) D(expr, name)
  else DD(D(expr, name), name, order - 1)
}

id.x<-order(x)
### the calculation with the assumpiton: weigths are uniformly distributed
mm<-mx%*%(y)/(mx%*%ii)

### this is the function to calculate the first order of the Gaussian
ddnorm <- function(x) eval(DD(expression(dnorm(x)), "x", order = 1))

### this is the function to build the matrixs for the optimization
yy<-t(matrix((rep(y,n)),n,n))

### the function coefficient matrixs
Amat<-(yy*n*as.vector(mx%*%ii)*(ddnorm(mx))-mx*yy*n*as.vector(ddnorm(mx)%*%ii))/(as.vector((mx%*%ii)^2))
### the corresponding matrix to calculate the optimization value
bvec<-rep(0,n)
ll<-rep(1,n)
### the first row is to make sure the sum of the probability equal to one
### the second rep vector is to make the inequation greater than zero
### the third rep vector is to make the element of the probabilities greater than zero
bvec<-c(1,rep(-0.0005,n),rep(-0.0005,n))

#### the basic struture of the varibales matrixs
#### which covers the inequation and the basic property of the probability greater than zero and sum equals to one
Amat<-matrix(rbind(ll,Amat),n+1,n)
ddm<-diag(1,n,n)
AAmat<-rbind(Amat,ddm)

### the D matrix to produce the objective function
Dmat<-diag(1,n,n)

### the function to cover the variable coefficients
dvec<-as.vector(rep(-2/n,n))

### to solve the objective function
sol <- solve.QP(Dmat,dvec,t(AAmat),bvec=bvec,meq=1)

### to solve the probability
p<-as.matrix(sol$solution)

### recalculate the estimation function plusing the probability : n*p
mmp<-(mx%*%(y*n*p))/((mx%*%ii))
plot(x[id.x],mm[id.x],col="blue",lty=2,lwd=2)
lines(x[id.x],mmp[id.x],type="l",col="red",lty=1,lwd=2)







