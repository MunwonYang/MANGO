rm(list = ls())
library(rlist)
library(foreach)
library(doParallel)
main.path <- "~/Real data analysis"
setwd(main.path)
load("EEG_dataset.RData")

x <- aperm(eegdata$X.mat,c(2,1,3))
y <- eegdata$y.vec
dimen <- dim(x)[1:2] # dimension of x
n <- dim(x)[3] # sample size
K <- 2 # order of tensor
x <- x[,31:256,]
x1 <- x[,,which(y==0)]
x2 <- x[,,which(y==1)]
x1 <- x1 - replicate(dim(x1)[K+1],apply(x1,c(1,2),mean))
x2 <- x2 - replicate(dim(x2)[K+1],apply(x2,c(1,2),mean))

cat("Dimension of Non alcoholic dataset:",dim(x1),'\n')
cat("Dimension of Alcoholic dataset:",dim(x2),'\n')
cat('\n')

source("~/Separate.R")
source("~/Basic Setting.R")

con1 <- Matrix.Surrogate(x1,type = "Constant");con2 <- Matrix.Surrogate(x2, type = "Constant")
w1 <- sweep(x1,3,con1,"/");w2 <- sweep(x2,3,con2,"/")

mat1 <- Matrix.Surrogate(x1,type = "Matrix");mat2 <- Matrix.Surrogate(x2,type = "Matrix")
w.mat1 <- array(0,dim(x1));w.mat2 <- array(0,dim(x2))
for(i in 1:dim(x1)[K+1]){
  eval(parse(text = paste0("w.mat1[",paste0(rep(",",K),collapse = ""),"i] = solve(diag(mat1[,i])) %*% x1[",paste0(rep(",",K),collapse = ""),'i]')))
}
for(i in 1:dim(x2)[K+1]){
  eval(parse(text = paste0("w.mat2[",paste0(rep(",",K),collapse = ""),"i] = solve(diag(mat2[,i])) %*% x2[",paste0(rep(",",K),collapse = ""),'i]')))
}

# Precision matrix estimation under MGGM
fit1 <- cv.Separate(x = w.mat1)
fit2 <- cv.Separate(x = w.mat2)
fit.adapt1.matrix <- cv.Separate.penalty(x = w.mat1, Omega = fit1$Omegahat, penalty = "adapt")
fit.adapt2.matrix <- cv.Separate.penalty(x = w.mat2, Omega = fit2$Omegahat, penalty = "adapt")

# Precision matrix estimation under U-MNAGO
fit1 <- cv.Separate(x = w1)
fit2 <- cv.Separate(x = w2)
fit.adapt1.constant <- cv.Separate.penalty(x = w1, Omega = fit1$Omegahat, penalty = "adapt")
fit.adapt2.constant <- cv.Separate.penalty(x = w2, Omega = fit2$Omegahat, penalty = "adapt")

# Precision matrix estimation under H-MANGO
fit1 <- cv.Separate(x = x1)
fit2 <- cv.Separate(x = x2)
fit.adapt1.original <- cv.Separate.penalty(x = x1, Omega = fit1$Omegahat, penalty = "adapt")
fit.adapt2.original <- cv.Separate.penalty(x = x2, Omega = fit2$Omegahat, penalty = "adapt")

# Hypothesis test(MGGM vs MANGO model)
cat("Hypothesis test(MGGM vs MANGO model)",'\n')
x1.stand <- atrans(x1, lapply(fit.adapt1.original$Omegahat, function(x) sqrtm(x)))
x2.stand <- atrans(x2, lapply(fit.adapt2.original$Omegahat, function(x) sqrtm(x)))

Type1.Non <- Matrix.tail.test(x1.stand, type = "type1", alternative = "two.sided")
Type1.Alc <- Matrix.tail.test(x2.stand, type = "type1", alternative = "two.sided")

cat("Test Statistic for Non Alcoholic in Global test:",Type1.Non$statistic,'\n')
cat("Test Statistic for Alcoholic in Global test:",Type1.Alc$statistic,'\n')
cat('\n')

# Hypothesis test(U-MANGO vs H-MANGO)
cat("Hypothesis test(U-MANGO vs H-MANGO)",'\n')
con1 <- Matrix.Surrogate(x1, type = "Constant")
con2 <- Matrix.Surrogate(x2, type = "Constant")
w1 <- sweep(x1,3,con1,"/");w2 <- sweep(x2,3,con2,"/")
w.stand1 <- atrans(w1,lapply(fit.adapt1.constant$Omegahat, function(x) sqrtm(x)))
w.stand2 <- atrans(w2,lapply(fit.adapt2.constant$Omegahat, function(x) sqrtm(x)))
w.mat1 <- Matrix.Surrogate(w.stand1,type = "Matrix");w.mat2 <- Matrix.Surrogate(w.stand2, type=  "Matrix")
w.con1 <- Matrix.Surrogate(w.stand1, type = "Constant");w.con2 <- Matrix.Surrogate(w.stand2, type = "Constant")

ncores=strtoi(Sys.getenv("SLURM_NTASKS"))
cl = makeCluster(ncores, type = "SOCK")
registerDoParallel(cl)

# Parameter generation for the Local test.
source("~/Basic Setting.R")
para1 <- Hypothesis.test.parameter(n = dim(x1)[K+1],dimen = dimen)
para2 <- Hypothesis.test.parameter(n = dim(x2)[K+1],dimen = dimen)
stopCluster(cl)

est1 <- mean((sweep(w.mat1^2,2,w.con1^2,"/")-1)^2);test.stat1 = (est1 - para1$mu) / para1$sigma
est2 <- mean((sweep(w.mat2^2,2,w.con2^2,"/")-1)^2);test.stat2 = (est2 - para2$mu) / para2$sigma
cat("Test statistic for Non Alcoholic in Local test:",test.stat1,'\n')
cat("Test statistic for Alcoholic in Local test:",test.stat2,'\n')
cat('\n')

if(test.stat1 < 0){
  pvalue1 <- pnorm(test.stat1) * 2
}else{
  pvalue1 <- pnorm(test.stat1,lower.tail = F) * 2
}
if(test.stat2 < 0){
  pvalue2 <- pnorm(test.stat2) * 2
}else{
  pvalue2 <- pnorm(test.stat2,lower.tail = F) * 2
}
cat("p-value for Non Alcoholic in Local test:",pvalue1,'\n')
cat("p-value for Alcoholic in Local test:",pvalue2,'\n')
cat('\n')
