rm(list = ls())
library(MASS)
library(expm)
library(rTensor)
library(tensr)
library(ramify)
library(rlist)
library(stats)
library(glasso)
library(Tlasso)
library(doParallel)
library(foreach)

Matrix.Normal = FALSE;Chi.Square = T;Contaminated <- FALSE;Laplace <- FALSE
nlambda <- 50;norm.off.diag <- F;ratio <- 0.8
gam <- 1;penalty <- 'adapt'
Separate.rate <- c(1e-5,1e-6)
Tlasso.rate <- c(1e-5,1e-6)
Gemini.rate <- c(1e-5,1e-6)
thr = 1e-4;maxit <- 1e3;thres <- 1e-5
nu = 4;i.i.d = TRUE;variation=TRUE;nu.min = NULL;nu.max <- NULL
source("~/Basic Setting.R")
source("~/Simulation Evaluation.R")
source("~/Simulation/Data Generation.R")
source("~/Simulation/Evaluation/Simulation Code/Methods/Separate.R")
source("~/Simulation/Evaluation/Simulation Code/Methods/Gemini.R")
source("~/Simulation/Evaluation/Simulation Code/Methods/Tlasso.R")
source("~/Simulation/Evaluation/Simulation Code/Chi Simulation.R")

ncores=strtoi(Sys.getenv("SLURM_NTASKS"))
cl = makeCluster(ncores, type = "SOCK")
registerDoParallel(cl)
##################################################################################################################################
n <- 200 # sample size
dimen <- c(80,32) # dimension of tensor
nvars <- prod(dimen) # Number of variables
K <- length(dimen) # Order of tensor

Omega <- array(list(),K)
Omega[[1]] <- NeighborOmega(dimen[1])
Omega[[2]] <- TR(dimen[2])
sparse = c(1,2)
Omega_Tilde <- Matrix.Identifiability(Omega)

# Separate
Heavy <- FALSE
start <- Sys.time()
C1.Separate <- Separate.simul(n = n, dimen = dimen,
                              Omega = Omega_Tilde,norm.off.diag = norm.off.diag,
                              gam = gam, penalty = penalty,
                              Heavy = Heavy, nu = nu,
                              rate = Separate.rate,
                              i.i.d = i.i.d, variation = variation,
                              nu.min = nu.min, nu.max = nu.max,
                              sparse = sparse)
end <- Sys.time()
end - start

# Cyclic
Heavy <- FALSE
start <- Sys.time()
C1.Cyclic <- Cyclic.simul(n = n, dimen = dimen, t = 5,
                          Omega = Omega_Tilde,norm.off.diag = norm.off.diag,
                          gam = gam, penalty = penalty,
                          Heavy = Heavy, nu = nu,
                          rate = Tlasso.rate, thres = thres,
                          i.i.d = i.i.d, variation = variation,
                          nu.min = nu.min, nu.max = nu.max,
                          sparse = sparse)
end <- Sys.time()
end - start

# Gemini
Heavy <- FALSE
start <- Sys.time()
C1.Gemini <- Gemini.simul(n = n, dimen = dimen,gam = gam,
                          Omega = Omega_Tilde,
                          norm.off.diag = norm.off.diag,
                          penalty = penalty, 
                          rate = Gemini.rate,
                          Heavy = Heavy, nu = nu, i.i.d = i.i.d,
                          nu.min = nu.min,nu.max = nu.max, variation = variation,
                          sparse = sparse)
end <- Sys.time()
end - start

# Heavy Separate(Each Row)
Heavy = "Matrix"
start <- Sys.time()
C1.Heavy.Matrix.Separate <- Separate.simul(n = n, dimen = dimen,
                                           Omega = Omega_Tilde,norm.off.diag = norm.off.diag,
                                           gam = gam, penalty = penalty,
                                           Heavy = Heavy, nu = nu,
                                           rate = Separate.rate,
                                           i.i.d = i.i.d, variation = variation,
                                           nu.min = nu.min, nu.max = nu.max,
                                           sparse = sparse)
end <- Sys.time()
end - start


# Heavy Separate(Each Matrix)
Heavy = "Constant"
start <- Sys.time()
C1.Heavy.Constant.Separate <-  Separate.simul(n = n, dimen = dimen,
                                              Omega = Omega_Tilde,norm.off.diag = norm.off.diag,
                                              gam = gam, penalty = penalty,
                                              Heavy = Heavy, nu = nu,
                                              rate = Separate.rate,
                                              i.i.d = i.i.d, variation = variation,
                                              nu.min = nu.min, nu.max = nu.max,
                                              sparse = sparse)
end <- Sys.time()
end - start
###################################################################################################################################
dimen <- c(80,256) # dimension of tensor
nvars <- prod(dimen) # Number of variables
K <- length(dimen) # Order of tensor

Omega <- array(list(),K)
Omega[[1]] <- NeighborOmega(dimen[1])
Omega[[2]] <- TR(dimen[2])
Omega_Tilde <- Matrix.Identifiability(Omega)

# Separate
Heavy <- FALSE
start <- Sys.time()
C2.Separate <- Separate.simul(n = n, dimen = dimen,
                              Omega = Omega_Tilde,norm.off.diag = norm.off.diag,
                              gam = gam, penalty = penalty,
                              Heavy = Heavy, nu = nu,
                              rate = Separate.rate,
                              i.i.d = i.i.d, variation = variation,
                              nu.min = nu.min, nu.max = nu.max,
                              sparse = sparse)
end <- Sys.time()
end - start

# Cyclic
Heavy <- FALSE
start <- Sys.time()
C2.Cyclic <- Cyclic.simul(n = n, dimen = dimen, t = 5,
                          Omega = Omega_Tilde,norm.off.diag = norm.off.diag,
                          gam = gam, penalty = penalty,
                          Heavy = Heavy, nu = nu,
                          rate = Tlasso.rate, thres = thres,
                          i.i.d = i.i.d, variation = variation,
                          nu.min = nu.min, nu.max = nu.max,
                          sparse = sparse)
end <- Sys.time()
end - start

# Gemini
Heavy <- FALSE
start <- Sys.time()
C2.Gemini <- Gemini.simul(n = n, dimen = dimen,gam = gam,
                          Omega = Omega_Tilde,
                          norm.off.diag = norm.off.diag,
                          penalty = penalty, 
                          rate = Gemini.rate,
                          Heavy = Heavy, nu = nu, i.i.d = i.i.d,
                          nu.min = nu.min,nu.max = nu.max, variation = variation,
                          sparse = sparse)
end <- Sys.time()
end - start

# Heavy Separate(Each Row)
Heavy = "Matrix"
start <- Sys.time()
C2.Heavy.Matrix.Separate <- Separate.simul(n = n, dimen = dimen,
                                           Omega = Omega_Tilde,norm.off.diag = norm.off.diag,
                                           gam = gam, penalty = penalty,
                                           Heavy = Heavy, nu = nu,
                                           rate = Separate.rate,
                                           i.i.d = i.i.d, variation = variation,
                                           nu.min = nu.min, nu.max = nu.max,
                                           sparse = sparse)
end <- Sys.time()
end - start

# Heavy Separate(Each Matrix)
Heavy = "Constant"
start <- Sys.time()
C2.Heavy.Constant.Separate <-  Separate.simul(n = n, dimen = dimen,
                                              Omega = Omega_Tilde,norm.off.diag = norm.off.diag,
                                              gam = gam, penalty = penalty,
                                              Heavy = Heavy, nu = nu,
                                              rate = Separate.rate,
                                              i.i.d = i.i.d, variation = variation,
                                              nu.min = nu.min, nu.max = nu.max,
                                              sparse = sparse)
end <- Sys.time()
end - start
##################################################################################################################################
dimen <- c(32,80) # dimension of tensor
nvars <- prod(dimen) # Number of variables
K <- length(dimen) # Order of tensor

Omega <- array(list(),K)
Omega[[1]] <- NeighborOmega(dimen[1])
Omega[[2]] <- TR(dimen[2])
Omega_Tilde <- Matrix.Identifiability(Omega)

# Separate
Heavy <- FALSE
start <- Sys.time()
R1.Separate <- Separate.simul(n = n, dimen = dimen,
                              Omega = Omega_Tilde,norm.off.diag = norm.off.diag,
                              gam = gam, penalty = penalty,
                              Heavy = Heavy, nu = nu,
                              rate = Separate.rate,
                              i.i.d = i.i.d, variation = variation,
                              nu.min = nu.min, nu.max = nu.max,
                              sparse = sparse)
end <- Sys.time()
end - start

# Cyclic
Heavy <- FALSE
start <- Sys.time()
R1.Cyclic <- Cyclic.simul(n = n, dimen = dimen, t = 5,
                          Omega = Omega_Tilde,norm.off.diag = norm.off.diag,
                          gam = gam, penalty = penalty,
                          Heavy = Heavy, nu = nu,
                          rate = Tlasso.rate, thres = thres,
                          i.i.d = i.i.d, variation = variation,
                          nu.min = nu.min, nu.max = nu.max,
                          sparse = sparse)
end <- Sys.time()
end - start

# Gemini
Heavy <- FALSE
start <- Sys.time()
R1.Gemini <- Gemini.simul(n = n, dimen = dimen,gam = gam,
                          Omega = Omega_Tilde,
                          norm.off.diag = norm.off.diag,
                          penalty = penalty, 
                          rate = Gemini.rate,
                          Heavy = Heavy, nu = nu, i.i.d = i.i.d,
                          nu.min = nu.min,nu.max = nu.max, variation = variation,
                          sparse = sparse)
end <- Sys.time()
end - start

# Heavy Separate(Each Row)
Heavy = "Matrix"
start <- Sys.time()
R1.Heavy.Matrix.Separate <- Separate.simul(n = n, dimen = dimen,
                                           Omega = Omega_Tilde,norm.off.diag = norm.off.diag,
                                           gam = gam, penalty = penalty,
                                           Heavy = Heavy, nu = nu,
                                           rate = Separate.rate,
                                           i.i.d = i.i.d, variation = variation,
                                           nu.min = nu.min, nu.max = nu.max,
                                           sparse = sparse)
end <- Sys.time()
end - start

# Heavy Separate(Each Matrix)
Heavy = "Constant"
start <- Sys.time()
R1.Heavy.Constant.Separate <-  Separate.simul(n = n, dimen = dimen,
                                              Omega = Omega_Tilde,norm.off.diag = norm.off.diag,
                                              gam = gam, penalty = penalty,
                                              Heavy = Heavy, nu = nu,
                                              rate = Separate.rate,
                                              i.i.d = i.i.d, variation = variation,
                                              nu.min = nu.min, nu.max = nu.max,
                                              sparse = sparse)
end <- Sys.time()
end - start
###################################################################################################################################
dimen <- c(256,80) # dimension of tensor
nvars <- prod(dimen) # Number of variables
K <- length(dimen) # Order of tensor

Omega <- array(list(),K)
Omega[[1]] <- NeighborOmega(dimen[1])
Omega[[2]] <- TR(dimen[2])
Omega_Tilde <- Matrix.Identifiability(Omega)

# Separate
Heavy <- FALSE
start <- Sys.time()
R2.Separate <- Separate.simul(n = n, dimen = dimen,
                              Omega = Omega_Tilde,norm.off.diag = norm.off.diag,
                              gam = gam, penalty = penalty,
                              Heavy = Heavy, nu = nu,
                              rate = Separate.rate,
                              i.i.d = i.i.d, variation = variation,
                              nu.min = nu.min, nu.max = nu.max,
                              sparse = sparse)
end <- Sys.time()
end - start

# Cyclic
Heavy <- FALSE
start <- Sys.time()
R2.Cyclic <- Cyclic.simul(n = n, dimen = dimen, t = 5,
                          Omega = Omega_Tilde,norm.off.diag = norm.off.diag,
                          gam = gam, penalty = penalty,
                          Heavy = Heavy, nu = nu,
                          rate = Tlasso.rate, thres = thres,
                          i.i.d = i.i.d, variation = variation,
                          nu.min = nu.min, nu.max = nu.max,
                          sparse = sparse)
end <- Sys.time()
end - start

# Gemini
Heavy <- FALSE
start <- Sys.time()
R2.Gemini <- Gemini.simul(n = n, dimen = dimen,gam = gam,
                          Omega = Omega_Tilde,
                          norm.off.diag = norm.off.diag,
                          penalty = penalty, 
                          rate = Gemini.rate,
                          Heavy = Heavy, nu = nu, i.i.d = i.i.d,
                          nu.min = nu.min,nu.max = nu.max, variation = variation,
                          sparse = sparse)
end <- Sys.time()
end - start

# Heavy Separate(Each Row)
Heavy = "Matrix"
start <- Sys.time()
R2.Heavy.Matrix.Separate <- Separate.simul(n = n, dimen = dimen,
                                           Omega = Omega_Tilde,norm.off.diag = norm.off.diag,
                                           gam = gam, penalty = penalty,
                                           Heavy = Heavy, nu = nu,
                                           rate = Separate.rate,
                                           i.i.d = i.i.d, variation = variation,
                                           nu.min = nu.min, nu.max = nu.max,
                                           sparse = sparse)
end <- Sys.time()
end - start

# Heavy Separate(Each Matrix)
Heavy = "Constant"
start <- Sys.time()
R2.Heavy.Constant.Separate <-  Separate.simul(n = n, dimen = dimen,
                                              Omega = Omega_Tilde,norm.off.diag = norm.off.diag,
                                              gam = gam, penalty = penalty,
                                              Heavy = Heavy, nu = nu,
                                              rate = Separate.rate,
                                              i.i.d = i.i.d, variation = variation,
                                              nu.min = nu.min, nu.max = nu.max,
                                              sparse = sparse)
end <- Sys.time()
end - start
stopCluster(cl)