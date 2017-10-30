setwd("~/PycharmProjects/cgRNA")
library(ggplot2)

file <- read.csv(file = "cross_products_all.csv",header = TRUE, sep=",")
head(file)
cp <- c(file$cross_products)

m<-mean(cp, na.rm=TRUE)
std<-sqrt(var(cp, na.rm=TRUE))
hist(cp, density=20, breaks=20, prob=TRUE, xlab="cross products distance to [0,0,0]", ylab="Probability", freq=F, 
     main=" C4'-N1/N9 vectors and lattice vectors\n cross products 'parallelity'")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
m
m+2*std
m-2*std

##########
file <- read.csv(file = "RMSDs_all.csv",header = TRUE, sep=",")
head(file)
intra <- c(file$intra_structure_rmsd)
inter <- c(file$mean_inter_structure_rmsd)

m<-mean(intra, na.rm=TRUE)
std<-sqrt(var(intra, na.rm=TRUE))
hist(intra, density=20, breaks=20, prob=TRUE, xlab="RMSDs", ylab="Probability", freq=F, 
     main="RMSDs between random P-N1/N9-P triangles from structures\n (all from structure vs all)")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
m
m+2*std
m-2*std

m<-mean(inter, na.rm=TRUE)
std<-sqrt(var(inter, na.rm=TRUE))
hist(inter, density=20, breaks=20, prob=TRUE, xlab="RMSDs", ylab="Probability", freq=F, 
     main="Mean RMSDs between P-N1/N9-P triangles in every structure")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
m
m+2*std
m-2*std

######
file <- read.csv(file = "DingKyu_dist_two_ways.csv",header = TRUE, sep=",")
head(file)
pur <- c(file$p_mc_pur_AG)
pir <- c(file$p_mc_pir_CU)
all_pmc <- c(file$p_mc_pir_CU,file$p_mc_pur_AG)
pur_pn <- c(file$p_n_pur_AG)
pir_pn <- c(file$p_n_pir_CU)
all_pn <- c(file$p_n_pir_CU, file$p_n_pur_AG)

pur_c4 <- c(file$c4_mc_pur_AG)
pir_c4 <- c(file$c4_mc_pir_CU)
all_c4mc <- c(file$c4_mc_pir_CU,file$c4_mc_pur_AG)
pur_c4n <- c(file$c4_n_pur_AG)
pir_c4n <- c(file$c4_n_pir_CU)
all_c4n <- c(file$c4_n_pir_CU, file$c4_n_pur_AG)

pur_c4cg <- c(file$c4_cg_pur)
pir_c4cu <- c(file$c4_cu_pir)
all <- c(file$c4_cg_pur, file$c4_cu_pir)

pur_pcg <- c(file$p_cg_pur)
pir_pcu <- c(file$p_cu_pir)
all_pcgu <- c(file$p_cg_pur, file$p_cu_pir)

c_c <- c(file$c4_c4)

# PURYNY A/G

m<-mean(pur, na.rm=TRUE)
std<-sqrt(var(pur, na.rm=TRUE))
hist(pur, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="P-MassCenter distance in Purines")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m
m+2*std
m-2*std

#> m
#[1] 5.429752
#> m+2*std
#[1] 6.40083
#> m-2*std
#[1] 4.458673

# PIRYMIDYNY C/U

m2<-mean(pir, na.rm=TRUE)
std2<-sqrt(var(pir,na.rm=TRUE))
hist(pir, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="P-MassCenter distance in Pyrimidines")
curve(dnorm(x, mean=m2, sd=std2), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m2
m2+2*std2
m2-2*std2

#> m2
#[1] 6.209618
#> m2+2*std2
#[1] 7.337574
#> m2-2*std2
#[1] 5.081663

# ALL

m3<-mean(all_pmc, na.rm=TRUE)
std3<-sqrt(var(all_pmc,na.rm=TRUE))
hist(all_pmc, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="P-MassCenter distances")
curve(dnorm(x, mean=m3, sd=std3), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m3
m3+2*std3
m3-2*std3

#> m3
#[1] 5.802005
#> m3+2*std3
#[1] 7.108527
#> m3-2*std3
#[1] 4.495483


# PURYNY A/G P-N

m<-mean(pur_pn, na.rm=TRUE)
std<-sqrt(var(pur_pn, na.rm=TRUE))
hist(pur_pn, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="P-N9 distance in Purines")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m
m+2*std
m-2*std

#ONE WAY DISTS
#> m
#[1] 5.429752
#> m+2*std
#[1] 6.40083
#> m-2*std
#[1] 4.458673
#TWO WAY DISTS
#> m
#[1] 5.500791
#> m+2*std
#[1] 6.597009
#> m-2*std
#[1] 4.404574

# PIRYMIDYNY C/U

m2<-mean(pir_pn, na.rm=TRUE)
std2<-sqrt(var(pir_pn,na.rm=TRUE))
hist(pir_pn, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="P-N1 distance in Pyrimidines")
curve(dnorm(x, mean=m2, sd=std2), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m2
m2+2*std2
m2-2*std2

#ONE WAY DISTS
#> m2
#[1] 5.38635
#> m2+2*std2
#[1] 6.225597
#> m2-2*std2
#[1] 4.547104
#TWO WAY DISTS
#> m2
#[1] 5.475633
#> m2+2*std2
#[1] 6.260225
#> m2-2*std2
#[1] 4.691042

# ALL

m3<-mean(all_pn, na.rm=TRUE)
std3<-sqrt(var(all_pn,na.rm=TRUE))
hist(all_pn, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="P-N1/N9 distances")
curve(dnorm(x, mean=m3, sd=std3), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m3
m3+2*std3
m3-2*std3

#ONE WAY DISTS
#> m3
#[1] 5.409035
#> m3+2*std3
#[1] 6.3205
#> m3-2*std3
#[1] 4.49757
#TWO WAY DISTS
#> m3
#[1] 5.488952
#> m3+2*std3
#[1] 6.451449
#> m3-2*std3
#[1] 4.526455

####################################################3

# PURYNY A/G

m<-mean(pur_c4, na.rm=TRUE)
std<-sqrt(var(pur_c4, na.rm=TRUE))
hist(pur_c4, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="C4-MassCenter distance in Purines")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m
m+2*std
m-2*std

#> m
#[1] 3.391569
#> m+2*std
#[1] 3.5908
#> m-2*std
#[1] 3.192339

# PIRYMIDYNY C/U

m2<-mean(pir_c4, na.rm=TRUE)
std2<-sqrt(var(pir_c4,na.rm=TRUE))
hist(pir_c4, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="C4-MassCenter distance in Pyrimidines")
curve(dnorm(x, mean=m2, sd=std2), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m2
m2+2*std2
m2-2*std2

#> m2
#[1] 4.918714
#> m2+2*std2
#[1] 5.152942
#> m2-2*std2
#[1] 4.684487

# ALL

m3<-mean(all_c4mc, na.rm=TRUE)
std3<-sqrt(var(all_c4mc,na.rm=TRUE))
hist(all_c4mc, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="C4-MassCenter distances")
curve(dnorm(x, mean=m3, sd=std3), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m3
m3+2*std3
m3-2*std3

#> m3
#[1] 4.120521
#> m3+2*std3
#[1] 5.661569
#> m3-2*std3
#[1] 2.579472


# PURYNY A/G C4-N

m<-mean(pur_c4n, na.rm=TRUE)
std<-sqrt(var(pur_c4n, na.rm=TRUE))
hist(pur_c4n, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="C4-N9 distance in Purines")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m
m+2*std
m-2*std

#> m
#[1] 3.391569
#> m+2*std
#[1] 3.5908
#> m-2*std
#[1] 3.192339

# PIRYMIDYNY C/U

m2<-mean(pir_c4n, na.rm=TRUE)
std2<-sqrt(var(pir_c4n,na.rm=TRUE))
hist(pir_c4n, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="C4-N1 distance in Pyrimidines")
curve(dnorm(x, mean=m2, sd=std2), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m2
m2+2*std2
m2-2*std2

#> m2
#[1] 3.409576
#> m2+2*std2
#[1] 3.580421
#> m2-2*std2
#[1] 3.238731

# ALL

m3<-mean(all_c4n, na.rm=TRUE)
std3<-sqrt(var(all_c4n,na.rm=TRUE))
hist(all_c4n, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="C4-N1/N9 distances")
curve(dnorm(x, mean=m3, sd=std3), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m3
m3+2*std3
m3-2*std3

#> m3
#[1] 3.400164
#> m3+2*std3
#[1] 3.587232
#> m3-2*std3
#[1] 3.213096

####################################################################################
# PURYNY

m<-mean(pur_c4cg, na.rm=TRUE)
std<-sqrt(var(pur_c4cg,na.rm=TRUE))
hist(pur_c4cg, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="Purines: C4-CG distances")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m
m+2*std
m-2*std

#> m
#[1] 3.81845
#> m+2*std
#[1] 4.322397
#> m-2*std
#[1] 3.314503

# PIRYMIDYNY

m<-mean(pir_c4cu, na.rm=TRUE)
std<-sqrt(var(pir_c4cu,na.rm=TRUE))
hist(pir_c4cu, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="Pyrimidines: C4-CU distances")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m
m+2*std
m-2*std

#> m
#[1] 3.646458
#> m+2*std
#[1] 4.004011
#> m-2*std
#[1] 3.288905

# ALL

m<-mean(all, na.rm=TRUE)
std<-sqrt(var(all,na.rm=TRUE))
hist(all, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="All C4-CG and C4-CU distances")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m
m+2*std
m-2*std

#> m
#[1] 3.736353
#> m+2*std
#[1] 4.208843
#> m-2*std
#[1] 3.263863

##################################### P - CG/CU ###############################################
# PURYNY

m<-mean(pur_pcg, na.rm=TRUE)
std<-sqrt(var(pur_pcg,na.rm=TRUE))
hist(pur_pcg, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="Purines: P-CG distances")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m
m+2*std
m-2*std

#> m
#[1] 6.001569
#> m+2*std
#[1] 10.38241
#> m-2*std
#[1] 1.620732

# PIRYMIDYNY

m<-mean(pir_pcu, na.rm=TRUE)
std<-sqrt(var(pir_pcu,na.rm=TRUE))
hist(pir_pcu, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="Pyrimidines: P-CU distances")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m
m+2*std
m-2*std

#> m
#[1] 6.010625
#> m+2*std
#[1] 11.43423
#> m-2*std
#[1] 0.5870153

# ALL

m<-mean(all_pcgu, na.rm=TRUE)
std<-sqrt(var(all_pcgu,na.rm=TRUE))
hist(all_pcgu, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="All P-CG and P-CU distances")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m
m+2*std
m-2*std

#> m
#[1] 6.004013
#> m+2*std
#[1] 10.68865
#> m-2*std
#[1] 1.319373

# INTER DISTS
#> m
#[1] 4.750205
#> m+2*std
#[1] 6.273313
#> m-2*std
#[1] 3.227097

################## c4-c4 ####################

m<-mean(c_c, na.rm=TRUE)
std<-sqrt(var(c_c,na.rm=TRUE))
hist(c_c, density=20, breaks=20, prob=TRUE, xlab="distance", ylab="Probability", freq=F, 
     main="All c4-c4 distances")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = m, col = "red")
#abline(v = m+2*std, col = "red")
#abline(v = m-2*std, col = "red")
m
m+2*std
m-2*std

#> m
#[1] 7.733181
#> m+2*std
#[1] 14.29276
#> m-2*std
#[1] 1.173598