## Multiple matrix with glasso

library(MASS)
library(glasso)
p<-2;q<-3
n<-4;it<-25

mu<-c(rep(0,p*q));mu
sigma<-diag(p*q);sigma
data.sa<-mvrnorm(n,mu,sigma);data.sa
rho<-glasso(omega0,rho=.01)$loglik

gl<-function(data.sa,p,q,n,it,rho)
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
	omega<-matrix(0,p,p)
	omega.old<-matrix(0,p,p)
	shi<-matrix(0,q,q) 	
	shi.old<-matrix(0,q,q)
	   	
   	for(i in 1:it)
   		{  
		for(j in 1:n)
       		{ 
			shi1<-shi+(t(xic[[j]])%*%solve(omega0)%*%xic[[j]])
       		}
		shi<-(1/(n*p))*shi1
		d1<-(shi-shi.old)
         	shi<-glasso(shi,rho=rho)$w
		shi.old<-shi   	

      	for(k in 1:n)
       		{ 
			omega<-omega+(xic[[k]]%*%solve(shi)%*%t(xic[[k]]))
       		}
		omega<-(1/(n*q))*omega
		d2<-(omega-omega.old)
		omega.old<-omega
      	omega0<-omega
		omega<-glasso(omega0,rho)$w
		rho<-glasso(omega0,rho=rho)$loglik

	if ((max(abs(d1))<0.001)&&(max(abs(d2))<0.001))
        break 
      else 
        next

      	}    
   	kr<-kronecker(shi,omega0)
   	it<-i
   	list(y=y,m=m,xbar=xbar,xic=xic,shi=shi,omega=omega0,rho=rho,kron=kr,Iteration=it)
	} 
 
gl(data.sa,2,3,4,10,rho)




## Multiple matrix with Addative glasso

library(MASS)
library(glasso)
p<-2;q<-3
n<-4;it<-25

mu<-c(rep(0,p*q))
sigma<-diag(p*q)
data.sa<-mvrnorm(n,mu,sigma)

agl<-function(data.sa,p,q,n,it,lamda,gamma)
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
	omega<-matrix(0,p,p)
	shi<-matrix(0,q,q)
	
	for(i in 1:it)
   		{  
		for(j in 1:n)
       		{ 
			shi<-shi+(t(xic[[j]])%*%solve(omega0)%*%xic[[j]])
       		}
		lamda<-nk*log(abs(omega0/)) - (1/shi[q,q])*shi
		shi<-(1/(n*p))*shi
		rhomat.shi<-lamda*matrix(1,q,q)/(pmax(abs(solve(shi))^gamma,1e-5))
		shi<-glasso(shi,rhomat.shi)$w
		
      	for(k in 1:n)
       		{ 
			omega<-omega+((xic[[k]]%*%solve(shi)%*%t(xic[[k]])))
       		}
		omega<-(1/(n*q))*omega
      	rhomat.om1<-lamda*matrix(1,p,p)/(pmax(abs(solve(omega))^gamma,1e-5))
      	omega0<-glasso(omega,rhomat.om1)$w
      	} 
   
   	kr<-kronecker(shi,omega0)
   	it<-i
   	list(y=y,m=m,xbar=xbar,xic=xic, shi=shi,omega=omega0,kron=kr,Iteration=it)
	} 
 
agl(data.sa,2,3,4,10,12,.5)




## Multiple matrix with SCAD


library(MASS)
library(glasso)
p<-2;q<-3
n<-4;it<-25

mu<-c(rep(0,p*q))
sigma<-diag(p*q)
data.sa<-mvrnorm(n,mu,sigma)

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
	omega<-matrix(0,p,p)
	shi<-matrix(0,q,q)
   	
	for(i in 1:it)
   		{  
		for(j in 1:n)
       		{ 
			shi<-shi+(t(xic[[j]])%*%solve(omega0)%*%xic[[j]])
       		}
		shi<-(1/(n*p))*shi
	#	shi<-(1/shi[q,q])*shi
      	rhomat.shi<-pmax(lam*((abs(solve(shi))<=lam)
				+pmax(a*lam-abs(solve(shi)),0)
				*abs((solve(shi))>lam)/(a-1)/lam),1e-4)
		shi<-glasso(shi,rhomat.shi)$w

		for(k in 1:n)
       		{ 
			omega<-(xic[[k]]%*%solve(shi)%*%t(xic[[k]]))
       		}
		omega<-(1/(n*q))*omega
      	rhomat.om1<-pmax(lam*((abs(solve(omega))<=lam)
				+pmax(a*lam-abs(solve(omega)),0)
				*abs((solve(omega))>lam)/(a-1)/lam),1e-4)

      	omega0<-glasso(omega,rhomat.om1)$w

		} 
   
   	kr<-kronecker(shi,omega0)
   	it<-i
   	list(y=y,m=m,xbar=xbar,xic=xic,shi=shi,omega=omega0,kron=kr,Iteration=it)
	} 
 
scad(data.sa,2,3,4,10,12,3.7)









