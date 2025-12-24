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

n <- 200
dist.type <- c(rep("Chi",7),"Contaminated")
nu <- c(4,5,7,10,20,50,100)
value <- c(1,5);prob <- c(0.9,0.1);different = variation = i.i.d= T
Test.type <- "Local";Simulation.type <- "Power"
dimen <- array(0,c(4,2))
dimen[1:2,1] = dimen[3:4,2] <- 80
dimen[1:2,2] = dimen[3:4,1] <- c(32,256)
K <- ncol(dimen)

Omega.list = Omega.Tilde.list <- array(list(),nrow(dimen))
for(i in 1:nrow(dimen)){
  Omega.list[[i]] <- array(list(),2)
  Omega.list[[i]][[1]] <- NeighborOmega(dimen[i,1]);Omega.list[[i]][[2]] <- TR(dimen[i,2])
  Omega.Tilde.list[[i]] <- Matrix.Identifiability(Omega.list[[i]])
}

Row1 <- Test.Estimation.simul(n = n,dimen = dimen, Omega.list = Omega.Tilde.list, dist.type = dist.type[1],
                              Simulation.type = Simulation.type, Test.type = Test.type,
                              nu = nu[1],variation = variation,i.i.d = i.i.d)
Row2 <- Test.Estimation.simul(n = n,dimen = dimen, Omega.list = Omega.Tilde.list, dist.type = dist.type[2],
                              Simulation.type = Simulation.type, Test.type = Test.type,
                              nu = nu[2],variation = variation,i.i.d = i.i.d)
Row3 <- Test.Estimation.simul(n = n,dimen = dimen, Omega.list = Omega.Tilde.list, dist.type = dist.type[3],
                              Simulation.type = Simulation.type, Test.type = Test.type,
                              nu = nu[3],variation = variation,i.i.d = i.i.d)
Row4 <- Test.Estimation.simul(n = n,dimen = dimen, Omega.list = Omega.Tilde.list, dist.type = dist.type[4],
                              Simulation.type = Simulation.type, Test.type = Test.type,
                              nu = nu[4],variation = variation,i.i.d = i.i.d)
Row5 <- Test.Estimation.simul(n = n,dimen = dimen, Omega.list = Omega.Tilde.list, dist.type = dist.type[5],
                              Simulation.type = Simulation.type, Test.type = Test.type,
                              nu = nu[5],variation = variation,i.i.d = i.i.d)
Row6 <- Test.Estimation.simul(n = n,dimen = dimen, Omega.list = Omega.Tilde.list, dist.type = dist.type[6],
                              Simulation.type = Simulation.type, Test.type = Test.type,
                              nu = nu[6],variation = variation,i.i.d = i.i.d)
Row7 <- Test.Estimation.simul(n = n,dimen = dimen, Omega.list = Omega.Tilde.list, dist.type = dist.type[7],
                              Simulation.type = Simulation.type, Test.type = Test.type,
                              nu = nu[7],variation = variation,i.i.d = i.i.d)
Row8 <- Test.Estimation.simul(n = n,dimen = dimen, Omega.list = Omega.Tilde.list, dist.type = dist.type[8],
                              Simulation.type = Simulation.type, Test.type = Test.type,
                              value = value, prob = prob, different = different)
end <- Sys.time();end - start

stopCluster(cl)