
				## Glasso Final Cross Validation

library(glasso)
library(MASS)

set.seed(100)
cv.gl<-function(n,p,k,m,r,d) 		# n rows,p columns,k fold,m rows in each fold
	{					#  r rho(s) and d diff. between two rho(s)
	rho<-0
	cv.tot<-numeric(r)
	
	for (i in 1:r)
  		{
		rho<-rho + d
		y<-matrix(rnorm(n*p),ncol=p)
		
		cv.s<-numeric(k)
		for(j in 1:k)
			{								# Elimination of kth-fold
			if (j == 1) {x<-y[(m*j+1):n, ]} else
			{if (j > 1 && j < k) {x<-rbind(y[1:(m*(j-1)), ] , y[(m*j+1):20, ])}} 
			{if (j == k) {x<-y[1:(m*(j-1)), ]}}
			s<-var(x)
			cv.s[j]<-glasso(s,rho)$loglik				# Single Cross Validation
			}
		sum(cv.s)					# Sum of all corss validation for one rho
		cv.tot[i]<-sum(cv.s)		# Crosso Validation for different rho
		}
	list(data=y,fold=x,var=s,ind.cv=cv.s,ind.sum=sum(cv.s),CV=cv.tot)
	}

cv.gl(20,3,5,4,10,.1)



##rough work


rho<-0
cv.tot<-numeric(10)
for (i in 1:10)
  
	{
	y<-matrix(rnorm(20*2),ncol=2)
	rho<-rho + .1	
	cv.s<-numeric(5)
	for(k in 1:5)
		{
		n<-5
		k<-4
		if (k == 1) {x<-y[(n*k+1):20, ]} else
		{if (k > 1 && k < 4) {x<-rbind(y[1:(n*(k-1)), ] , y[(n*k+1):20, ])}} 
		{if (k == 4) {x<-y[1:(n*(k-1)), ]}}
		s<-var(x)
		cv.s[k]<-glasso(s,rho)$loglik		
		}
	sum(cv.s)
	cv.tot[i]<-sum(cv.s)
	
	}
sum(cv.tot)
list(data=y,fold=x,var=s,ind.cv=cv.s,ind.sum=sum(cv.s), CV=sum(cv.tot))





			## Adaptive Glasso Final Cross Validation



library(glasso)
library(MASS)

cv.al1<-function(n,p,k,m,r,d,gamma)
	{
	y<-matrix(rnorm(n*p),ncol=p)

	esti.L1.out<-NULL
	
	n<-dim(y)[1]
	p<-dim(y)[2]

	lambda<-0
	cv.tot<-numeric(r)

	for(i in 1:r)
		{
		lambda<-lambda+d
		cov.sa<-cov(y)
		esti.L1.out<-glasso(cov.sa,lambda)
		rhomat<-lambda*matrix(1,p,p)/(pmax(abs(esti.L1.out$wi)^gamma,1e-5))
		
		cv.s.aL1<-numeric(k)
		for(j in 1:k)
			{								# Elimination of kth-fold
			if (j == 1) {x<-y[(m*j+1):n, ]} else
			{if (j > 1 && j < k) {x<-rbind(y[1:(m*(j-1)), ] , y[(m*j+1):20, ])}} 
			{if (j == k) {x<-y[1:(m*(j-1)), ]}}
			s<-cov(x)
			cv.s.aL1[j]<-glasso(s,rhomat)$loglik		# Single Cross Validation
			}
		cv.tot[i]<-sum(cv.s.aL1)	# Crosso Validation for different rho

		}
	list(data=y,ind.cv=cv.s.aL1,CV=cv.tot)
	}

cv.al1(20,3,5,4,10,.01,0.5)



						## SCAD Final Cross Validation




library(glasso)
library(MASS)

cv.scad<-function(n,p,k,m,r,d,a)
	{
	y<-matrix(rnorm(n*p),ncol=p)

	est.L1.out<-NULL
		
	n<-dim(y)[1]
	p<-dim(y)[2]

	lam<-0
	cv.tot<-numeric(r)

	for(i in 1:r)
		{
		lam<-lam+d
		cov.sa<-cov(y)

		est.L1.out<-glasso(cov.sa,lam)
		rhomat<-pmax(lam*((abs(est.L1.out$wi)<=lam)
				+pmax(a*lam-abs(est.L1.out$wi),0)
				*abs((est.L1.out$wi)>lam)/(a-1)/lam),1e-4)

		cv.s.scad<-numeric(k)
		for(j in 1:k)
			{								# Elimination of kth-fold
			if (j == 1) {x<-y[(m*j+1):n, ]} else
			{if (j > 1 && j < k) {x<-rbind(y[1:(m*(j-1)), ] , y[(m*j+1):20, ])}} 
			{if (j == k) {x<-y[1:(m*(j-1)), ]}}
			s<-cov(x)
			cv.s.scad[j]<-glasso(s,rhomat)$loglik		# Single Cross Validation
			}
		cv.tot[i]<-sum(cv.s.scad)			# Crosso Validation for different rho
		}
	list(data=y,ind.cv=cv.s.scad,CV=cv.tot)
	}

cv.scad(20,3,5,4,10,.01,3.7)




				## SCAD Final Cross Validation(multiple steps)


library(glasso)
library(MASS)

cv.scad<-function(n,p,k,m,r,d,a,g)
	{
	y<-matrix(rnorm(n*p),ncol=p)

	est.L1.out<-NULL
	cv.scad<-NULL
	
	n<-dim(y)[1]
	p<-dim(y)[2]

	lam<-0
	CV<-numeric(r)

	for(i in 1:r)
		{
		lam<-lam+d
		cov.sa<-cov(y)

		cv.scad<-numeric(g)
  		for (l in 1:g)
  			{
			est.L1.out<-glasso(cov.sa,lam)
			rhomat<-pmax(lam*((abs(est.L1.out$wi)<=lam)
					+pmax(a*lam-abs(est.L1.out$wi),0)
					*abs((est.L1.out$wi)>lam)/(a-1)/lam),1e-4)

			cv.s.scad<-numeric(k)
			for(j in 1:k)
				{								# Elimination of kth-fold
				if (j == 1) {x<-y[(m*j+1):n, ]} else
				{if (j > 1 && j < k) {x<-rbind(y[1:(m*(j-1)), ] , y[(m*j+1):20, ])}} 
				{if (j == k) {x<-y[1:(m*(j-1)), ]}}
				s<-cov(x)
				cv.s.scad[j]<-glasso(s,rhomat)$loglik		# Single Cross Validation
				}
			cv.scad[l]<-cv.s.scad[j]
			}
		CV[i]<-sum(cv.scad)			# Cross Validation for different rho
		}
	list(data=y,cv.s.scad=cv.s.scad,cv.scad=cv.scad,CV=CV)
	}

cv.scad(20,3,5,4,10,.01,3.7,5)




##  SCAD Final Cross Validation (multiple steps with while statement and convergence)


library(glasso)
library(MASS)

cv.scad<-function(n,p,k,m,r,d,a)
	{
	y<-matrix(rnorm(n*p),ncol=p)

	est.L1.out<-NULL
		
	n<-dim(y)[1]
	p<-dim(y)[2]

	lam<-0
	cv.tot<-numeric(r)

	for(i in 1:r)
		{
		lam<-lam+d
		cov.sa<-cov(y)

		est.L1.out<-glasso(cov.sa,lam)
		rhomat<-pmax(lam*((abs(est.L1.out$wi)<=lam)
				+pmax(a*lam-abs(est.L1.out$wi),0)
				*abs((est.L1.out$wi)>lam)/(a-1)/lam),1e-4)

		cv.s.scad<-numeric(k)
		for(j in 1:k)
			{								# Elimination of kth-fold
			if (j == 1) {x<-y[(m*j+1):n, ]} else
			{if (j > 1 && j < k) {x<-rbind(y[1:(m*(j-1)), ] , y[(m*j+1):20, ])}} 
			{if (j == k) {x<-y[1:(m*(j-1)), ]}}
			s<-cov(x)
			
			est.scad.old<-NULL
			est.scad.new<-est.L1.out$wi
			eps<-1
			count<-1
			while(eps >1e-4)
				{
				est.scad.old<-est.scad.new
				est.scad.new<-glasso(s,rhomat)$wi
    				count<-count+1
  				eps<-max(abs(est.scad.old-est.scad.new))
				}
			cat(" Iterration number ")
		
			cv.s.scad[j]<-glasso(est.scad.old, rhomat)$loglik	# Single Cross Validation
			}
		cv.tot[i]<-sum(cv.s.scad)			# Crosso Validation for different rho
		}
	list(data=y,ind.cv=cv.s.scad,CV=cv.tot,iteration=count)
	}

cv.scad(20,3,5,4,10,.01,3.7)







