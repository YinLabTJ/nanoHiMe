data <- read.table(file = "coefficients.input", sep = "\t", header = F)
x1<-data[,5]
x2<-data[,6]
y1<-data[,2]
y2<-data[,3]
lm1.model<-lm(y1~x1)
lm2.model<-lm(y2~x2-1)
coefficients(lm1.model)
coefficients(lm2.model)
