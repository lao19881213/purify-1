library(lattice)
library(gridExtra)
path<-"/Users/lwolz/Documents/purify/data/test/chirp5pc/"

fileext<-c("Gmatconv_test_noshift.txt", "Gmatconv_noshift.txt")
skiprows<-0
D<-c()
for(i in 1:length(fileext)){
	file<-paste0(path, fileext[i])
	data<-read.table(file, nrow=1, skip=skiprows)

	mat<-matrix((abs(unlist(data))), 128, 128)
	#mat[which(mat<(-40))]=0
	p<-levelplot(mat)
	if(i==1){ mat1<-mat
		plot1<-p}
	else{ mat2<-mat
		plot2<-p}
}

plot3<-levelplot(mat1-mat2)
grid.arrange(plot1, plot2, plot3)