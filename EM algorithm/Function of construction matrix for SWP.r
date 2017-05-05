#################################################################################
#################################################################################
#				Reverse Sweep Function
# Input	 
# RSW		: Name of the Reverse sweep function
# G		: Matrix that has to be reverse swept
# k		: Reverse swept on row k and column k
# p		: Total number of columns in G
# Output	
# H		: Reverse swept matrix

#################################################################################

RSW<-function(G, k, p) 					
	{
	H<-matrix(0,p,p)
	for(j in 1:p)
		{
		for(l in 1:p)
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
#				Construction of Conditional Covariance
# Input
# library(ggm)	: Library needed for Sweep(swp) function
# C.sub		: Conditional matrix
# p			: Total number of columns in G or H
# Yobs		: Y matrix without missing values
# sw1, sw2 & sw3	: 2x2, 3x3 & 4x4 matrices have to be swept respectively
# swp			: Sweep(swp) function
# Output
# RS			: Reverse sweep and sweep matrix
# C			: Conditional Covariance matrix

#################################################################################

library(ggm)
	
C.sub<-function(p)						
	{
	i<-p-3
	Yobs<-Y[apply(Y, 1, function(x)!any(is.na(x))), , drop=F]
		mu<-mean(Yobs[ ,1])
		sig1<-var(Yobs[ ,1])
		sw1<-cbind(rbind(-1,mu), c(t(mu),sig1))
		SWP.1<-swp(sw1, (i+1))
	i<-p-2
		fit<-lm(Yobs[ ,i] ~ Yobs[ ,(i-1)])
			beta1<-coef(fit)
			R1<-residuals(fit)
			sig22.1<-var(R1)
			sw2<-cbind(rbind(SWP.1, beta1), c(t(beta1), sig22.1))
			SWP.2<-swp(sw2, (i+1))
	i<-p-1		
		fit<-lm(Yobs[ ,i] ~ Yobs[ ,(i-1)] + Yobs[ ,(i-2)])
			beta2<-coef(fit)
			R2<-residuals(fit)
			sig33.12<-var(R2)
			sw3<-cbind(rbind(SWP.2, beta2), c(t(beta2), sig33.12))
			SWP.1<-swp(sw3, (i+1))
			RSW.1<-RSW(sw3, i, p)
	    	      RS<-RSW(RSW.1, (i+1), p)
	RS
	}
p<-4
C<-C.sub(p)[2:p,2:p];C

#################################################################################
#################################################################################






