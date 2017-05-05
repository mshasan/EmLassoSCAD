## Glasso

library(glasso)
library(MASS)

mu<-rep(0,5)
sig<-diag(5)
data.sa<-mvrnorm(3,mu,sig)
cov.sa<-cov(data.sa)					# Used Covariance matrix
g<-glasso(cov.sa, rho=.01)				
gl<-glasso(cov.sa,rho=.02,w.init=a$w,wi.init=a$wi)	# Used value of first glasso


set.seed(100)
x<-matrix(rnorm(5*2),ncol=2)
s<-var(x)
g<-glasso(s,rho=.01)					# Used Covariance matrix
gl<-glasso(s,rho=.02,w.init=a$w,wi.init=a$wi)	# Used value of first glasso



s<-c(10,1,5,4,10,2,6,10,3,10)
S<-matrix(0,nrow=4,ncol=4)
S[row(S)>=col(S)]<-s
diag(S)<-10
zero<-matrix(c(1,3,2,4),ncol=2,byrow=TRUE)
g<-glasso(S,0,zero=zero)



## Adaptive Lasso

library(glasso)
library(MASS)

al1<-function(data.sa,lambda,gamma)
	{
	esti.L1<-NULL
	esti.aL1<-NULL

	n1<-dim(data.sa)[1]
	p1<-dim(data.sa)[2]

	cov.sa<-cov(data.sa)
	esti.L1.out<-glasso(cov.sa,lambda)
	rhomat<-lambda*matrix(1,p1,p1)/(pmax(abs(esti.L1.out$wi)^gamma,1e-5))
	esti.aL1<-glasso(cov.sa,rhomat)$wi
	print(esti.aL1)
	}

mu<-rep(0,5)
sig<-diag(5)
data.sa<-mvrnorm(3,mu,sig)
al1(data.sa,.01,0.5)



## SCAD


scadrightderv <- function(lamhat, a, lam)
{
 pmax(lam*((lamhat<=lam)+pmax(a*lam-lamhat,0)*(lamhat>lam)/(a-1)/lam),1e-4)
}


##  SCAD (one step)

library(glasso)
library(MASS)

scad<-function(data.sa, lam, a)
	{
	est.L1=NULL
	est.scad=NULL
	cov.sa=cov(data.sa)
	n=dim(data.sa)[1]
	p=dim(data.sa)[2]

		
		est.L1.out=glasso(cov.sa,lam)
		rhomat<-pmax(lam*((abs(est.L1.out$wi)<=lam)
				+pmax(a*lam-abs(est.L1.out$wi),0)
				*abs((est.L1.out$wi)>lam)/(a-1)/lam),1e-4)
		
		print(rhomat)
		est.scad=glasso(cov.sa,rhomat)$wi;
		print(est.scad)
		
	}

mu<-rep(0,5)
sig<-diag(5)
data.sa<-mvrnorm(3,mu,sig)
scad(data.sa,.1,3.7)




##  SCAD (multiple steps)

library(glasso)
library(MASS)

scad<-function(data.sa, lam, a, m)
	{
	est.L1=NULL
	est.scad=NULL
	cov.sa=cov(data.sa)
	n=dim(data.sa)[1]
	p=dim(data.sa)[2]

	est.scad<-c()
  	for (i in 1:m)
  		{
		est.L1.out=glasso(cov.sa,lam)
		rhomat<-pmax(lam*((abs(est.L1.out$wi)<=lam)
				+pmax(a*lam-abs(est.L1.out$wi),0)
				*abs((est.L1.out$wi)>lam)/(a-1)/lam),1e-4)
		
		print(rhomat)
		est.scad=glasso(cov.sa,rhomat)$wi;
		print(est.scad)
		}
 	est.scad

	}


mu<-rep(0,5)
sig<-diag(5)
data.sa<-mvrnorm(3,mu,sig)
scad(data.sa,.8,3.7,5)




##  SCAD (multiple steps with while statement and convergence)

library(glasso)
library(MASS)

scad<-function(data.sa, lam, a)
	{
	est.L1<-NULL
	est.scad.old<-NULL
	cov.sa<-cov(data.sa)
	n<-dim(data.sa)[1]
	p<-dim(data.sa)[2]

	est.L1.out<-glasso(cov.sa,lam)
	est.scad.new<-est.L1.out$wi
	eps<-1
	count<-1
	while(eps >1e-4)
		{
		est.scad.old<-est.scad.new
		rhomat<-pmax(lam*((abs(est.L1.out$wi)<=lam)
				+pmax(a*lam-abs(est.L1.out$wi),0)
				*abs((est.L1.out$wi)>lam)/(a-1)/lam),1e-4)
		
		est.scad.new<-glasso(cov.sa,rhomat)$wi
    		count<-count+1
  		eps<-max(abs(est.scad.old-est.scad.new))
		}
	cat("number of iteration is")
 	list(iteration=count,scad= est.scad.old)
	}
mu<-rep(0,5)
sig<-diag(5)
data.sa<-mvrnorm(3,mu,sig)
scad(data.sa,.8,3.7)


##  SCAD (multiple steps with while statement and convergence
	and without convergence of diagonal elements)

library(glasso)
library(MASS)

scad<-function(data.sa, lam, a)
	{
	est.L1<-NULL
	est.scad.old<-NULL
	cov.sa<-cov(data.sa)
	n<-dim(data.sa)[1]
	p<-dim(data.sa)[2]

	est.L1.out<-glasso(cov.sa,lam,penalize.diagonal=FALSE)
	est.scad.new<-est.L1.out$wi
	eps<-1
	count<-1
	while(eps >1e-4)
		{
		est.scad.old<-est.scad.new
		rhomat<-pmax(lam*((abs(est.L1.out$wi)<=lam)
				+pmax(a*lam-abs(est.L1.out$wi),0)
				*abs((est.L1.out$wi)>lam)/(a-1)/lam),1e-4)
		
		est.scad.new<-glasso(cov.sa,rhomat,penalize.diagonal=FALSE)$wi
    		count<-count+1
  		eps<-max(abs(est.scad.old-est.scad.new))
		}
	cat("number of iteration is")
 	list(iteration=count,scad= est.scad.old)
	}
mu<-rep(0,5)
sig<-diag(5)
data.sa<-mvrnorm(3,mu,sig)
scad(data.sa,.8,3.7)






