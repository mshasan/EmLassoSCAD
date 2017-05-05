library(ggm)

## Reverse sweep function

RSW<-function(g, k, q) 	# kk sweep and q variable
	{
	h<-matrix(0,q,q)
	for(j in 1:q)
		{
		for(l in 1:q)
			{
			if (j == k && l == k) {h[j,l] = -1/g[k,k]} else
			{if (j != k && l == k) {h[j,l] = -g[j,k]/g[k,k]}} 
			{if (j == k && l != k) {h[j,l] = -g[k,l]/g[k,k]}}
			{if (j != k && l != k ) {h[j,l] = (g[j,k]*g[k,l])/g[k,k]}}
			}
		}
	h
	}

## Sweep 1

y1<-c(35,35,40,6,20,35,35);y1
mu1<-mean(y1);mu1
sig1<-var(y1);sig1
sw1<-cbind(rbind(-1,mu1), c(t(mu1),sig1));sw1
SWP.1<-swp(sw1, 2);SWP.1


##Sweep 2

y2<-c(3.5,4.9,30,2.7,2.8,10.9,8.0);y2
fit2 <- lm(y2 ~ y1);fit2
beta2<-coef(fit2);beta2
R2<-residuals(fit2);R2 
sig22.1<-var(R2);sig22.1
sw2<-cbind(rbind(SWP.1, beta2), c(t(beta2), sig22.1));sw2
SWP.2<-swp(sw2, 3);SWP.2
#RSW.1<-RSW(sw2, 2, 3);RSW.1

##Sweep 3

y3<-c(2.8,2.7,3.21,2.81,2.88,2.90,3.28);y3
fit3 <- lm(y3 ~ y1+y2);fit3
beta3<-coef(fit3);beta3
R3<-residuals(fit3);R3 
sig33.12<-var(R3);sig33.12
sw3<-cbind(rbind(SWP.2, beta3), c(t(beta3), sig33.12));sw3
#SWP.3<-swp(sw3, 4);SWP.3
RSW.1<-RSW(sw3, 2, 4);RSW.1
RSW.12<-RSW(RSW.1, 3, 4);RSW.12



y1<-c(35,35,40,6,20,35,35)
y2<-c(3.5,4.9,30,2.7,2.8,10.9,8.0)
y3<-c(2.8,2.7,3.21,2.81,2.88,2.90,3.28)

Y<-as.data.frame(cbind(y1, y2, y3))
X<-as.matrix(Y);X

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
C(X)
	

