rm(list = ls())
library(rlist)
library(foreach)
library(doParallel)
setwd("~/Real data analysis/")
load("EEG_dataset.RData")
source("~/Basic_Setting.R")

x <- aperm(eegdata$X.mat,c(2,1,3))
x <- x[,31:256,]
y <- eegdata$y.vec
dimen <- dim(x)[1:2] # dimension of x
n <- dim(x)[3] # sample size
K <- 2 # order of tensor
x <- x - replicate(dim(x)[K+1],apply(x,c(1,2),mean))

source("~/Separate.R")
source("~/Basic Setting.R")

mat1 <- Matrix.Surrogate(x,mode = 1,type = "Matrix")
mat2 <- Matrix.Surrogate(x,mode = 2,type = "Matrix")

library(latex2exp)
lwd = 1;cex.main = 3; cex.lab =  2.5;cex.axis = 2
par(mfrow = c(1, 2), mar = c(5.5, 6, 4, 2), oma = c(0, 0, 0, 0),mgp   = c(4, 1, 0))  

# First plot
plot(density(mat1[1,]), main = TeX("$X = \\diag(\\Gamma)W$, with $\\Gamma \\in R^{64}$"),
     xlim = c(0, 20), ylab = "Density", xlab = TeX("$\\widehat{\\delta}_j$ for j = 1,...,64"),
     ylim = c(0, max(apply(mat1, 1, function(x) max(density(x)$y))) + 0.05),
     lwd = lwd, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis)
for (i in 2:dimen[1]) lines(density(mat1[i,]), lwd = lwd)

# Second plot
plot(density(mat2[1,]), main = TeX("$X = W\\diag(\\Gamma)$ with $\\Gamma \\in R^{226}$"),
     xlim = c(0, 12), ylab = "", xlab = TeX("$\\widehat{\\delta}_j$ for j = 31,...,256"),
     ylim = c(0, max(apply(mat2, 1, function(x) max(density(x)$y))) + 0.05),
     lwd = lwd, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis)
for (i in 2:dimen[2]) lines(density(mat2[i,]), lwd = lwd)



