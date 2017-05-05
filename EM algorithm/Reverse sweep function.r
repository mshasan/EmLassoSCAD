library(ggm)
g<-matrix(c(2,3,5,3,4,7,5,7,8), nrow=3);g

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


RSW(g, 3, 3)

a<-matrix(c(2,3,4,5,2,6,7,8,9,2,3,5,3,4,7,5),4,byrow=T);a
b<-matrix(c(5,3,7,2),2,byrow=T);b
	
RSW(a, 3, 4)



		