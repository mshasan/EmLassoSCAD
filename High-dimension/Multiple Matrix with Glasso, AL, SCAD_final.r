##############################################################################
##############################################################################
# 				Multiple matrix with Graphical lasso (glasso)
# Input	
# library(MASS)	: Basic R library
# library(glasso)	: Graphical lasso library
# data.sa		: Data matrix
# p & q		: Columns of sliced matrices
# n			: Number of rows
# it 			: Number of Iterations
# Output
# shi & omega0	: Sliced matrices			 
# kr			:Kronecker product
##############################################################################

library(MASS)
library(glasso)

gl<-function(data.sa,p,q,n,it)
	{  
	y<-matrix(0,n,p*q)
   	m<-list(n)
  
   	x0<-numeric(p*q)
   	for(i in 1:n)
   		{  
		y[i,]<-data.sa[i,]   		
      	x0<-x0+y[i,]
      	m[[i]]<-matrix(y[i,],p,q)
  	 	}

   	xbar<-x0/n
   	xic<-list(n)
   	for(i in 1:n)
   		{  
		xic[[i]]<-m[[i]]-xbar
   		}

   	omega0<-diag(p)
	a<-glasso(omega0,rho=.01)
   	
   	for(i in 1:it)
   		{  
		shi<-matrix(0,q,q)
      	for(j in 1:n)
       		{ 
			shi<-shi+(1/(n*p))*(t(xic[[j]])%*%solve(omega0)%*%xic[[j]])
       		}
      	shi<-(1/shi[q,q])*shi
      	aa<-glasso(shi,rho=.01)    	

      	omega<-matrix(0,p,p)
      	for(k in 1:n)
       		{ 
			omega<-omega+((1/(n*q))*(xic[[k]]%*%solve(shi)%*%t(xic[[k]])))
       		}
      	omega0<-omega
		aaa<-glasso(omega0,rho=.01)

      	}    
   	kr<-kronecker(shi,omega0)
   	it<-i
   	list(xbar=xbar,psi=shi,omega=omega0,kron=kr,Iteration=it)
	} 
##############################################################################
#				Example of glasso
p<-2;q<-3
n<-4;it<-25
mu<-c(rep(0,p*q))
sigma<-diag(p*q)
data.sa<-mvrnorm(n,mu,sigma)
gl(data.sa,2,3,4,10)

##############################################################################
##############################################################################
# 			Multiple matrix with Addative glasso
# Input	
# library(MASS)	: Basic R library
# library(glasso)	: Graphical lasso library
# data.sa		: Data matrix
# p & q		: Columns of sliced matrices
# n			: Number of rows
# it 			: Number of Iterations
# lamda & gamma	: Two parameters
# Output
# shi & omega0	: Sliced matrices			 
# kr			: Kronecker product
##############################################################################

library(MASS)
library(glasso)

agl<-function(data.sa,p,q,n,it,lamda,gamma)
	{  
	y<-matrix(0,n,p*q)
   	m<-list(n)
  
   	x0<-numeric(p*q)
   	for(i in 1:n)
   		{  
		y[i,]<-data.sa[i,]   		
      	x0<-x0+y[i,];x0
      	m[[i]]<-matrix(y[i,],p,q)
  	 	}

   	xbar<-x0/n
   	xic<-list(n)
   	for(i in 1:n)
   		{  
		xic[[i]]<-m[[i]]-xbar
   		}

	omega0<-diag(p)
	est.om.L1.out<-glasso(omega0,lamda)
	rhomat.om<-lamda*matrix(1,p,p)/(pmax(abs(est.om.L1.out$wi)^.5,1e-5))
	omega0<-glasso(omega0,rhomat.om)$w
   	
	for(i in 1:it)
   		{  
		shi<-matrix(0,q,q)
      	for(j in 1:n)
       		{ 
			shi<-shi+(1/(n*p))*(t(xic[[j]])%*%solve(omega0)%*%xic[[j]])
       		}
		rhomat.shi<-lamda*matrix(1,q,q)/(pmax(abs(solve(shi))^gamma,1e-5))
		shi<-glasso(shi,rhomat.shi)$w
		shi<-(1/shi[q,q])*shi
      	
      	omega<-matrix(0,p,p)
      	for(k in 1:n)
       		{ 
			omega<-omega+((1/(n*q))*(xic[[k]]%*%solve(shi)%*%t(xic[[k]])))
       		}
		rhomat.om1<-lamda*matrix(1,p,p)/(pmax(abs(solve(omega))^gamma,1e-5))
      	omega0<-glasso(omega,rhomat.om1)$w
      	}    
   	kr<-kronecker(shi,omega0)
   	it<-i
   	list(xbar=xbar,psi=shi,omega=omega0,kron=kr,Iteration=it)
	} 
#################################################################################
#				Example of Addaptive glasso
p<-2;q<-3
n<-4;it<-25
mu<-c(rep(0,p*q))
sigma<-diag(p*q)
data.sa<-mvrnorm(n,mu,sigma)
agl(data.sa,2,3,4,10,12,.5)

#################################################################################
#################################################################################
# 		Multiple matrix with Smoothly Clipped Absolute Deviation (SCAD)
# Input	
# library(MASS)	: Basic R library
# library(glasso)	: Graphical lasso library
# data.sa		: Data matrix
# p & q		: Columns of sliced matrices
# n			: Number of rows
# it 			: Number of Iterations
# lam	& a		: Two parameters
# Output
# shi & omega0	: Sliced matrices			 
# kr			: Kronecker product
#################################################################################

library(MASS)
library(glasso)

scad<-function(data.sa,p,q,n,it,lam,a)
	{  
	y<-matrix(0,n,p*q)
   	m<-list(n)
  
   	x0<-numeric(p*q)
   	for(i in 1:n)
   		{  
		y[i,]<-data.sa[i,]   		
      	x0<-x0+y[i,]
      	m[[i]]<-matrix(y[i,],p,q)
  	 	}

   	xbar<-x0/n
   	xic<-list(n)
   	for(i in 1:n)
   		{  
		xic[[i]]<-m[[i]]-xbar
   		}

	omega0<-diag(p)
	est.om.L1.out<-glasso(omega0,lam)
	rhomat.om<-pmax(lam*((abs(est.om.L1.out$wi)<=lam)
				+pmax(a*lam-abs(est.om.L1.out$wi),0)
				*abs((est.om.L1.out$wi)>lam)/(a-1)/lam),1e-4)
	omega0<-glasso(omega0,rhomat.om)$w
   	
	for(i in 1:it)
   		{  
		shi<-matrix(0,q,q)
      	for(j in 1:n)
       		{ 
			shi<-shi+(1/(n*p))*(t(xic[[j]])%*%solve(omega0)%*%xic[[j]])
       		}
		rhomat.shi<-pmax(lam*((abs(solve(shi))<=lam)
				+pmax(a*lam-abs(solve(shi)),0)
				*abs((solve(shi))>lam)/(a-1)/lam),1e-4)

		shi<-glasso(shi,rhomat.shi)$w
		shi<-(1/shi[q,q])*shi
      	
      	omega<-matrix(0,p,p)
      	for(k in 1:n)
       		{ 
			omega<-omega+((1/(n*q))*(xic[[k]]%*%solve(shi)%*%t(xic[[k]])))
       		}
		rhomat.om1<-pmax(lam*((abs(solve(omega))<=lam)
				+pmax(a*lam-abs(solve(omega)),0)
				*abs((solve(omega))>lam)/(a-1)/lam),1e-4)

      	omega0<-glasso(omega,rhomat.om1)$w
      	}    
   	kr<-kronecker(shi,omega0)
   	it<-i
   	list(xbar=xbar,psi=shi,omega=omega0,kron=kr,Iteration=it)
	} 
#################################################################################
# 					Exapmle of SCAD
p<-2;q<-3
n<-4;it<-25
mu<-c(rep(0,p*q))
sigma<-diag(p*q)
data.sa<-mvrnorm(n,mu,sigma)
scad(data.sa,2,3,4,10,12,3.7)

#################################################################################
#################################################################################









