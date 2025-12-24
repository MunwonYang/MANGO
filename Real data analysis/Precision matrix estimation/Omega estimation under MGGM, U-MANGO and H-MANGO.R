rm(list = ls())
main.path <- "~/Real data analysis"
setwd(main.path)
library(igraph)
load("EEG_dataset.RData")

x <- aperm(eegdata$X.mat,c(2,1,3))
x <- x[,31:256,]
y <- eegdata$y.vec
dimen <- dim(x)[1:2] # dimension of x
n <- dim(x)[3] # sample size
K <- 2 # order of tensor

x1 <- x[,,which(y==0)]
x2 <- x[,,which(y==1)]
x1 <- x1 - replicate(dim(x1)[K+1],apply(x1,c(1,2),mean))
x2 <- x2 - replicate(dim(x2)[K+1],apply(x2,c(1,2),mean))

source("~/Separate.R")
source("~/Basic Setting.R")

con1 <- Matrix.Surrogate(x1,type = "Constant");con2 <- Matrix.Surrogate(x2, type = "Constant")
w1 <- sweep(x1,3,con1,"/");w2 <- sweep(x2,3,con2,"/")

# Eventually, I will conduct the precision matrix under H-MANGO.
mat1 <- Matrix.Surrogate(x1,type = "Matrix");mat2 <- Matrix.Surrogate(x2,type = "Matrix")
w.mat1 <- array(0,dim(x1));w.mat2 <- array(0,dim(x2))
for(i in 1:dim(x1)[K+1]){
  eval(parse(text = paste0("w.mat1[",paste0(rep(",",K),collapse = ""),"i] = solve(diag(mat1[,i])) %*% x1[",paste0(rep(",",K),collapse = ""),'i]')))
}
for(i in 1:dim(x2)[K+1]){
  eval(parse(text = paste0("w.mat2[",paste0(rep(",",K),collapse = ""),"i] = solve(diag(mat2[,i])) %*% x2[",paste0(rep(",",K),collapse = ""),'i]')))
}

fit1 <- cv.Separate(x = w.mat1)
fit2 <- cv.Separate(x = w.mat2)
fit.adapt1.matrix <- cv.Separate.penalty(x = w.mat1, Omega = fit1$Omegahat, penalty = "adapt")
fit.adapt2.matrix <- cv.Separate.penalty(x = w.mat2, Omega = fit2$Omegahat, penalty = "adapt")


fit1 <- cv.Separate(x = w1)
fit2 <- cv.Separate(x = w2)
fit.adapt1.constant <- cv.Separate.penalty(x = w1, Omega = fit1$Omegahat, penalty = "adapt")
fit.adapt2.constant <- cv.Separate.penalty(x = w2, Omega = fit2$Omegahat, penalty = "adapt")

fit1 <- cv.Separate(x = x1)
fit2 <- cv.Separate(x = x2)
fit.adapt1.original <- cv.Separate.penalty(x = x1, Omega = fit1$Omegahat, penalty = "adapt")
fit.adapt2.original <- cv.Separate.penalty(x = x2, Omega = fit2$Omegahat, penalty = "adapt")

separate.first <- array(list(),c(3,2))
Assumption <- c("Normal","Constant","Matrix")
type <- c("Non.Alcoholic","Alcoholic")

dimnames(separate.first) <- list(
  Assumption, type
)
separate.first[1,1][[1]] <- fit.adapt1.original$Omegahat[[1]];separate.first[1,2][[1]] <- fit.adapt2.original$Omegahat[[1]]
separate.first[2,1][[1]] <- fit.adapt1.constant$Omegahat[[1]];separate.first[2,2][[1]] <- fit.adapt2.constant$Omegahat[[1]]
separate.first[3,1][[1]] <- fit.adapt1.matrix$Omegahat[[1]];separate.first[3,2][[1]] <- fit.adapt2.matrix$Omegahat[[1]]
list.save(separate.first,"Separate.rdata")