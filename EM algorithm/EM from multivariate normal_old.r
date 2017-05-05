library(Matrix)
library(MASS)
library(ggm)

## Reverse sweep function

RSW<-function(g, k, p) 	# kk sweep and p variable
	{
	h<-matrix(0,p,p)
	for(j in 1:p)
		{
		for(l in 1:p)
			{
			if (j == k && l == k) {h[j,l] = -1/g[k,k]} else
			{if (j != k && l == k) {h[j,l] = -g[j,k]/g[k,k]}} 
			{if (j == k && l != k) {h[j,l] = -g[k,l]/g[k,k]}}
			{if (j != k && l != k ) {h[j,l] = (g[j,k]*g[k,l])/g[k,k]}}
			}
		}
	h
	}



em.norm<-function(Y)
	{
	Yobs<-Y[apply(Y, 1, function(x)!any(is.na(x))), , drop=F];Yobs	# Obs withour missing values
	Ymis<-Y[apply(Y, 1, function(x)any(is.na(x))), , drop=F];Ymis	# With missing values
	Y<-rbind(Yobs[ ,1:3], Ymis[ ,1:3]);Y

	n<-dim(Y)[1];n
	r<-dim(Yobs)[1];r
	pi<-22/7

	mut<-apply(Yobs, 2, mean);mut						# initial values
	sigmat<-cov(Yobs)*(r-1)/r;sigmat
	
	ln.L<-function(y, mu, sigma2, n)					# log-likelihood function
		{
		dbar<-apply(Yobs,1,'-',mut);dbar
		-.5*sum(log(det(sigma2)))-.5*(t(dbar)%*%solve(sigmat)%*%dbar)
		}

	lltm1<-ln.L(Yobs, mut, sigmat, n);lltm1				# Initial log-lik 
	repeat{					
		EY<-apply(Yobs, 2, sum) + (n-r)*mut ;EY			# E-step
 		mut1<-EY/n;mut1	
											# M-step
		dbar1<-apply(Yobs,1,'-',mut1);dbar1
		sigmat1<-(1/n)*(dbar1%*%t(dbar1) + RSW.12[2:4,2:4]);sigmat1

		mut<-mut1;mut							# Update parameters
		sigmat<-sigmat1;sigmat	
		
		llt<-ln.L(Yobs, mut, sigmat, n);llt				# current log-like
		cat(mut, sigmat, llt, "\n")					# Print all values
		
		if(abs(lltm1-llt)< 0.00001) break				# Stop converged
		lltm1<-llt
		}
	list(iteration=count,mu=mut, sigma=sigmat)			# Fill missing values with new
	}

	





mu0<-c(3, 1, 4);mu0
sig0<-matrix(c(6,1,-2,1,13,4,-2,4,4), nrow=3, byrow=T);sig0
Y<-mvrnorm(5,mu0,sig0);Y
Y[4, ]<-NA;Y
em.norm(Y)



