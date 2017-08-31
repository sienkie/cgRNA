setwd("~/PycharmProjects/cgRNA")
library(ggplot2)

file <- read.csv(file = "Kyu_dists.csv",header = TRUE, sep=",")
file2 <- read.csv(file = "Ding_dists.csv",header = TRUE, sep=",")

pp <- c(file$distPP,file2$distPP)
cp <- c(file$distCP,file2$distCP, file$distPC,file2$distPC)
cc <- c(file$distCC,file2$distCC)
ccpp <- c(file$distPP,file2$distPP,file$distCC,file2$distCC)

# PP

m<-mean(pp)
std<-sqrt(var(pp))
hist(pp, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", ylim=c(0, 0.7), 
     main="P-P distance in Ding&Kyu data set")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m
m+2*std
m-2*std

# CP

m2<-mean(cp)
std2<-sqrt(var(cp))
hist(cp, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", ylim=c(0, 2), 
     main="C-P distance in Ding&Kyu data set")
curve(dnorm(x, mean=m2, sd=std2), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m2, col = "red")
#abline(v = m2+2*std2, col = "red")
#abline(v = m2-2*std2, col = "red")
m2
m2+2*std2
m2-2*std2

# CC

m3<-mean(cc)
std3<-sqrt(var(cc))
hist(cc, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", ylim=c(0, 2), 
     main="C-C distance in Ding&Kyu data set")
curve(dnorm(x, mean=m3, sd=std3), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m2, col = "red")
#abline(v = m2+2*std2, col = "red")
#abline(v = m2-2*std2, col = "red")
m3
m3+2*std3
m3-2*std3

# CC&PP

m4<-mean(ccpp)
std4<-sqrt(var(ccpp))
hist(ccpp, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", ylim=c(0, 1), 
     main="C-C/P-P distance in Ding&Kyu data set")
curve(dnorm(x, mean=m4, sd=std4), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m2, col = "red")
#abline(v = m2+2*std2, col = "red")
#abline(v = m2-2*std2, col = "red")
m4
std4
m4+2*std4
m4-2*std4
