#################################################################################
#################################################################################
#				Reverse Sweep Function
# Input	 
# RSW		: Name of the Reverse sweep function
# G		: Matrix that has to be reverse swept
# k		: Reverse swept on row k and column k
# q		: Total number of columns in G
# Output	
# H		: Reverse swept matrix

#################################################################################

RSW<-function(G, k, q) 					
	{
	H<-matrix(0,q,q)
	for(j in 1:q)
		{
		for(l in 1:q)
			{
			if (j == k && l == k) {H[j,l] = -1/G[k,k]} else
			{if (j != k && l == k) {H[j,l] = -G[j,k]/G[k,k]}} 
			{if (j == k && l != k) {H[j,l] = -G[k,l]/G[k,k]}}
			{if (j != k && l != k ) {H[j,l] = (G[j,k]*G[k,l])/G[k,k]}}
			}
		}
	H
	}
#################################################################################
#################################################################################
#			Function of Conditional Covariance
# Input
# library(ggm)	: Library needed for Sweep(swp) function
# C			: Name of the function
# X			: Matrix without missing values
# p			: Total number of columns in the data matrix
# Output
# SWP			: Swept matrix
# Rev.swp		: Reverse swept
# Rev.swp[]		: Conditional Covariance matrix

#################################################################################

library(ggm)
	
library(ggm)
C<-function(X)
	{
	p<-ncol(X)
	mu1<-mean(X[,1])
	sig1<-var(X[,1])
	sw1<-cbind(rbind(-1,mu1), c(t(mu1),sig1))
	SWP.1<-swp(sw1, 2)
	SWP<-SWP.1

	for(i in 2:p)
		{
		fit<-lm(X[,i]~X[,1:(i-1)])
		beta<-coef(fit)
		res<-residuals(fit)
		sigma<-var(res)
		sw<-cbind(rbind(SWP, beta), c(t(beta), sigma))
		SWP<-swp(sw, (i+1))
		}
	sw
		
	for(j in 2:p)
		{
		Rev.swp<-RSW(sw, j,(p+1))
		sw<-Rev.swp
		}
	Rev.swp
	Rev.swp[2:(p+1), 2:(p+1)]
	}

#################################################################################
#################################################################################
#					 Log likelihood function
# Input
# ln.L	: Function's name	
# y		: Data matrix
# mu		: Multivariate mean of the data
# sigma	: Multivariate covariance matrix
# n		: Number of observations or rows
# log(-1)	: 1.36

#################################################################################
#ln.L<-function(y, mu, sigma, n)				
#	{
#	dbar<-apply(y,1,'-',mu)
#	if (det(sigma) < 0) {-1.36 -.5*n*(log(abs(det(sigma))))-.5*(t(dbar)%*%solve(sigma)%*%dbar)} else
#	{if(det(sigma) >= 0){-.5*n*(log(abs(det(sigma))))-.5*(t(dbar)%*%solve(sigma)%*%dbar)}} 
#	}

ln.L<-function(y, mu, sigma, n)				
		{
		dbar<-apply(y,1,'-',mu)
		-.5*n*(log(abs(det(sigma))))-.5*(t(dbar)%*%solve(sigma)%*%dbar)
		}
#################################################################################
#################################################################################
# 					EM algorithm function
# em.norm	: Function's name
# Y		: Given dataset or matrix with missing values
# mut		: Initial multivariate mean
# sigmat	: Initial covariane matrix
# p		: Number of columns
# Yobs	: Data without missing values
# Ymis	: Data with only missing values
# EY		: E-step
# Output
# Iteration	: Number of iteration needed for convergence
# mu		: Estimated mean vector
# sigma	: Estimated covariance matrix

#################################################################################
	
em.norm<-function(Y)
	{
	p<-dim(Y)[2]
	n<-dim(Y)[1]
	pi<-22/7

	Yobs<-Y[apply(Y, 1, function(x)!any(is.na(x))), , drop=F]	
	Ymis<-Y[apply(Y, 1, function(x)any(is.na(x))), , drop=F]	
	Y<-rbind(Yobs[ ,1:p], Ymis[ ,1:p])

	r<-dim(Yobs)[1]
	mut<-apply(Yobs, 2, mean)					
	sigmat<-cov(Yobs)*(r-1)/r
#	mut<-c(sample(1:p))
#	sigmat<-matrix(c(sample(1:(p*p))), p, p)
	
	lltm1<-ln.L(Yobs, mut, sigmat, n)
	lltm1.new<-lltm1
	eps<-1
	count<-1
	while(eps >1e-4)
		{
		lltm1.old<-lltm1.new
		EY<-apply(Yobs, 2, sum) + (n-r)*mut 		
 		mut1<-EY/n

		for(j in 1:p)
			{
			Y[,j][is.na(Y[,j])]<-mean(Y[,j], na.rm=T)
			}
		Y
					
		dbar1<-apply(Yobs,1,'-',mut1)
		sigmat1<-(1/n)*(dbar1%*%t(dbar1) + C(Yobs))

		mut<-mut1							
		sigmat<-sigmat1	
		
		lltm1.new<-ln.L(Yobs, mut, sigmat, n)
		count<-count+1
		eps<-max(abs(lltm1.old-lltm1.new))
		}
	list(iteration=count,mu=mut, sigma=sigmat)		
	}
#################################################################################
#################################################################################
# 					Example-1
# Input		
# library(MASS)	: Basic R library
# library(Matrix)	: Library needed for multivariate data
# library(mvnmle)	: Library need to find out MLE of Multivriate data
# mu0			: Multivariate mean 
# sigma0		: Covariance matrix
# Y			: Multivariate data before missing values
# Output		
# Y[4:7, ]<-NA	: Previous multivariate data with some missing values
# em.norm()		: Created function
# mlest()		: R built-in function

#################################################################################

library(MASS)
library(Matrix)
library(mvnmle)

mu0<-c(3, 1, 4)
sigma0<-matrix(c(6,1,-2,1,13,4,-2,4,4), nrow=3, byrow=T)
Y<-mvrnorm(10,mu0,sigma0)
Y<-Y
Y[4:7, ]<-NA;Y
em.norm(Ym)

mlest(Ym)		# R built-in function

#################################################################################
#			Estimate by using original data
colMeans(Y)
cov(Y)

################################################################################
################################################################################
#					Example-2
x1<-c(7,1,11,11,7,11,3,1,2,NA,NA,NA,NA)
x2<-c(26,29,56,31,52,55,71,31,54,NA,NA,NA,NA)
x3<-c(6,15,8,8,6,9,17,22,18,4,23,9,8)
x4<-c(60,52,20,47,33,22,NA,NA,NA,NA,NA,NA,NA)
x5<-c(78.5,74.3,104.3,87.6,95.9,109.2,102.7,72.5,93.1,115.9,83.8,113.3,109.4)
Wmis<-as.matrix(cbind(x1,x2,x3,x4,x5))

em.norm(Wmis)

mlest(Wmis)		# R built-in function

#################################################################################
#			Estimate by using original data

x1<-c(7,1,11,8,7,11,3,1,2,21,1,11,10)
x2<-c(26,29,30,31,52,55,71,31,54,47,40,66,68)
x3<-c(6,15,8,8,7,9,17,16,18,4,23,9,8)
x4<-c(60,52,20,25,33,22,6,44,22,26,34,12,12)
x5<-c(78.5,74.3,104.3,87.6,95.9,109.2,102.7,72.5,93.1,115.9,83.8,113.3,109.4)
Wobs<-as.matrix(cbind(x1,x2,x3,x4,x5))

colMeans(Wobs)
cov(Wobs)

################################################################################
################################################################################
#					Example-3

Ymis<-as.matrix(read.table("EMmiss.txt", header=F))
em.norm(Ymis)

mlest(Ymis)		# R built-in function

###############################################################################
#			Estimate by using original data

Yobs<-as.matrix(read.table("EMdata.txt", header=F))
colMeans(Yobs)
cov(Yobs)

################################################################################
################################################################################
#				Hotelling T est
# Input
# Hypothesis	: Two means are equal


alpha<-.05
p<-6
n1<-nrow(Ymis)
n2<-nrow(Yobs)
X1.bar<-em.norm(Ymis)$mu
X2.bar<-colMeans(Yobs)
D.bar<-X1.bar-X2.bar
S1<-em.norm(Ymis)$sigma
S2<-cov(Yobs)
Spl<-((n1-1)/(n1+n2-2))*S1 + ((n2-1)/(n1+n2-2))*S2

H.T<-((n1*n2)/(n1+n2))*t(D.bar)%*%solve(Spl)%*%D.bar;H.T
F.val<-(((n1+n2-2)*p)/(n1+n2-p-1))*qf(1-alpha, p, n1+n2-p-1);F.val
reject<-H.T > F.val; reject

##############################################################################
#				Box's M test
g<-2
M<-(n1+n2-2)*log(det(Spl))-(n1-1)*log(det(S1))-(n2-1)*log(det(S2));M
U<-(1/(n1-1)+1/(n2-1)-1/(n1+n2-2))*((2*p^2+3*p-1)/(6*(p+1)*(g-1)));U
C<-(1-U)*M;C
Chi.sq<-qchisq(1-alpha, p*(p+1)*(g-1)/2);Chi.sq




