Mat <- matrix(sample(0:5, 20, replace = TRUE), ncol = 4);Mat
is.na(Mat) <- Mat == 0; Mat
#Mat[is.na(Mat)] <- apply(Mat,2,mean, na.rm=T);Mat
mean<-apply(Mat,2,mean, na.rm=T);mean
Mat[is.na(Mat[,2])]<-mean[2];Mat

p<-4
n<-5
Mat <- matrix(sample(0:5, 20, replace = TRUE), ncol = 4)
is.na(Mat) <- Mat == 0; Mat
mean<-apply(Mat,2,mean, na.rm=T);mean

for(j in 1:p)
	{
	Mat[,j][is.na(Mat[,j])]<-mean(Mat[,j], na.rm=T)
	}
Mat



