rm(list = ls())
library(MASS)
library(expm)
library(rTensor)
library(tensr)
library(ramify)
library(stats)
library(doParallel)
library(foreach)

ncores=strtoi(Sys.getenv("SLURM_NTASKS"))
cl = makeCluster(ncores, type = "SOCK")
registerDoParallel(cl)

source("~/Basic Setting.R")
source("~/Simulation/Hypothesis test/Test on Estimation.R")

dist.type <- "Normal";Test.type <- "Global";Simulation.type <- "Size"

dimen <- array(0,c(4,2))
dimen[1:2,1] = dimen[3:4,2] <- 80
dimen[1:2,2] = dimen[3:4,1] <- c(32,256)
sam.size <- c(10,25,50,100,200)

Omega.list = Omega.Tilde.list <- array(list(),nrow(dimen))
for(i in 1:nrow(dimen)){
  Omega.list[[i]] <- array(list(),2)
  Omega.list[[i]][[1]] <- NeighborOmega(dimen[i,1]);Omega.list[[i]][[2]] <- TR(dimen[i,2])
  Omega.Tilde.list[[i]] <- Matrix.Identifiability(Omega.list[[i]])
}

n1 <- Test.Estimation.simul(n = sam.size[1],dimen = dimen, ncores = ncores,
                            dist.type = dist.type, Test.type = Test.type, Simulation.type = Simulation.type,
                            Omega.list = Omega.Tilde.list)
n2 <- Test.Estimation.simul(n = sam.size[2],dimen = dimen, ncores = ncores,
                            dist.type = dist.type, Test.type = Test.type, Simulation.type = Simulation.type,
                            Omega.list = Omega.Tilde.list)
n3 <- Test.Estimation.simul(n = sam.size[3],dimen = dimen, ncores = ncores,
                            dist.type = dist.type, Test.type = Test.type, Simulation.type = Simulation.type,
                            Omega.list = Omega.Tilde.list)
n4 <- Test.Estimation.simul(n = sam.size[4],dimen = dimen, ncores = ncores,
                            dist.type = dist.type, Test.type = Test.type, Simulation.type = Simulation.type,
                            Omega.list = Omega.Tilde.list)
n5 <- Test.Estimation.simul(n = sam.size[5],dimen = dimen, ncores = ncores,
                            dist.type = dist.type, Test.type = Test.type, Simulation.type = Simulation.type,
                            Omega.list = Omega.Tilde.list)
stopCluster(cl)