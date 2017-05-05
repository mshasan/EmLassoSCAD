em.norm<-function(Y)
	{
	Yobs<-Y[!is.na(Y)];Yobs					# Obs withour missing values
	Ymis<-Y[is.na(Y)];Ymis					# With missing values
	n<-length(c(Yobs, Ymis));n
	r<-length(Yobs);r
	pi<-22/7

	mut<-mean(Yobs);mut					# initial values
	sigmat<-var(Yobs)*(r-1)/r;sigmat
	
	ln.L<-function(y, mu, sigma2, n)			# log-likelihood function
		{
		-.5*n*log(2*pi*sigma2)-.5*sum((y-mu)^2)/sigma2
		}

	lltm1<-ln.L(Yobs, mut, sigmat, n);lltm1		# Initial log-lik 
	repeat{					
		EY<-sum(Yobs)+(n-r)*mut;EY			# E-step
		EY2<-sum(Yobs^2)+(n-r)*(mut^2+sigmat);EY2
		
		mut1<-EY/n;mut1					# M-step
		sigmat1<-EY2/n - mut1^2;sigmat1

		mut<-mut1;mut					# Update parameters
		sigmat<-sigmat1;sigmat	
		
		llt<-ln.L(Yobs, mut, sigmat, n);llt		# current log-like
		cat(mut, sigmat, llt, "\n")			# Print all values
		
		if(abs(lltm1-llt)< 0.00001) break		# Stop converged
		lltm1<-llt
		}
	list(mu=mut, sigma=sigmat)				# Fill missing values with new
	}

x<-rnorm(20,5)
x[16:20]<-NA;x
em.norm(x)

mean(x, na.rm=T)		
mean(x^2, na.rm=T)-mean(x, na.rm=T)^2			# variance
var(x, na.rm=T)*(n-1)/n						# variance another way
crossprod(x[!is.na(x)]-mean(x, na.rm=T))[1,1]/15	# yet another way