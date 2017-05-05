				## Cross Validation for Glasso
## Row Elimination (every row considered as a fold) for one rho
		

set.seed(100)
cv.gl<-function(p,q)					# p rows and q culumns
	{
	rho<-.01
	y<-matrix(rnorm(p*q),ncol=q)
		
	cv.s<-numeric(p)
	for(k in 1:p)
		{
		if (k >= p) {x<-y[0:(k-1), ]} else		# Elimination of kth-fold
		{if (k < p) {x<-rbind(y[0:(k-1), ] , y[(k+1):p, ])}} 
		s<-var(x)
		cv.s[k]<-glasso(s,rho)$loglik			# Single Cross Validation

		}
	sum(cv.s)				# Sum of all corss validation for one rho

	list(data=y,fold=x,var=s,ind.cv=cv.s,sum=sum(cv.s))
	}
cv.gl(5,2)





				## Fold elimination for one rho


set.seed(100)
cv.gl<-function(p,q,m,n) 	# p rows,q columns,m fold and n rows in each fold
	{
	rho<-.01
	y<-matrix(rnorm(p*q),ncol=q)
		
	cv.s<-numeric(m)
	for(k in 1:m)
		{							# Elimination of kth-fold
		if (k == 1) {x<-y[(n*k+1):20, ]} else
		{if (k > 1 && k < m) {x<-rbind(y[1:(n*(k-1)), ] , y[(n*k+1):20, ])}} 
		{if (k == m) {x<-y[1:(n*(k-1)), ]}}
		s<-var(x)
		cv.s[k]<-glasso(s,rho)$loglik			# Single Cross Validation
		}
	sum(cv.s)					# Sum of all corss validation for one rho

	list(data=y,fold=x,var=s,ind.cv=cv.s,sum=sum(cv.s))
	}
cv.gl(20,2,4,5)





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
		
		cv.tot[i]<-sum(cv.s)		# Crosso Validation for different rho
		}
	list(data=y,fold=x,var=s,ind.cv=cv.s,CV=cv.tot)
	}

cv.gl(20,3,5,4,10,.1)





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



mu0<-c(3, 1, 4);mu0
sig0<-matrix(c(6,1,-2,1,13,4,-2,4,4), nrow=3, byrow=T);sig0
Y<-mvrnorm(20,mu0,sig0);Y

##K-fold cross-validation
library(DAAG)
cv.lm(df=mydata, fit, m=4) # 3 fold cross-validation





