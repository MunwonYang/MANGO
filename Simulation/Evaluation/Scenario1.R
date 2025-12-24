rm(list = ls())
library(MASS)
library(expm)
library(rTensor)
library(rlist)
library(tensr)
library(ramify)
library(stats)
library(glasso)
library(Tlasso)
library(doParallel)
library(foreach)

### Matrix Normal Case
Matrix.Normal = TRUE;Chi.Square = FALSE;Contaminated <- FALSE;Laplace <- FALSE
ratio <- 0.8;nlambda <- 50;norm.off.diag <- F
gam <- 1;penalty <- 'adapt'
Separate.rate <- c(1e-5,1e-6)
Tlasso.rate <- c(1e-5,1e-6)
Gemini.rate <- c(1e-5,1e-6)
thr <- 1e-4;maxit <- 1e3

source("~/Basic Setting.R")
source("~/Simulation Evaluation.R")
source("~/Simulation/Data Generation.R")
source("~/Simulation/Evaluation/Simulation Code/Methods/Separate.R")
source("~/Simulation/Evaluation/Simulation Code/Methods/Gemini.R")
source("~/Simulation/Evaluation/Simulation Code/Methods/Tlasso.R")
source("~/Simulation/Evaluation/Simulation Code/Normal Simulation.R")


ncores=strtoi(Sys.getenv("SLURM_NTASKS"))
cl = makeCluster(ncores, type = "SOCK")
registerDoParallel(cl)
####################################################################################################################################
n <- 200 # sample size
dimen <- c(80,32) # dimension of tensor
nvars <- prod(dimen) # Number of variables
K <- length(dimen) # Order of tensor

Omega <- array(list(),K)
Omega[[1]] <- TR(dimen[1])
Omega[[2]] <- AR(ratio,dimen[2])
Omega_Tilde <- Matrix.Identifiability(Omega)

sparse <- 1

# Separate
Heavy <- FALSE
start <- Sys.time()
C1.Separate <- Separate.simul(n = n, dimen = dimen,Omega = Omega_Tilde,
                              gam = gam, penalty = penalty,
                              norm.off.diag = norm.off.diag,Heavy = Heavy,
                              rate = Separate.rate,sparse = sparse)
end <- Sys.time()
end - start

# Cyclic
Heavy <- FALSE
start <- Sys.time()
C1.Cyclic <- Cyclic.simul(n = n, dimen = dimen, Omega = Omega_Tilde,
                          gam = gam, penalty = penalty,
                          t = 5, rate = Tlasso.rate,
                          sparse = sparse, Heavy = Heavy)
end <- Sys.time()
end - start

# Gemini
Heavy <- FALSE
start <- Sys.time()
C1.Gemini <- Gemini.simul(n = n, dimen = dimen,gam = gam,
                          norm.off.diag = norm.off.diag,
                          penalty = penalty, rate = Gemini.rate, 
                          Omega = Omega_Tilde, sparse = sparse,
                          Heavy = Heavy)
end <- Sys.time()
end - start

# Heavy Separate(Each Row)
Heavy = "Matrix"
start <- Sys.time()
C1.Heavy.Matrix.Separate <- Separate.simul(n = n, dimen = dimen,Omega = Omega_Tilde,
                                           gam = gam, penalty = penalty,
                                           norm.off.diag = norm.off.diag,Heavy = Heavy,
                                           rate = Separate.rate,sparse = sparse)
end <- Sys.time()
end - start

# Heavy Separate(Each Matrix)
Heavy = "Constant"
start <- Sys.time()
C1.Heavy.Constant.Separate <- Separate.simul(n = n, dimen = dimen,Omega = Omega_Tilde,
                                             gam = gam, penalty = penalty,
                                             norm.off.diag = norm.off.diag,Heavy = Heavy,
                                             rate = Separate.rate,sparse = sparse)
end <- Sys.time()
end - start
###################################################################################################################################
dimen <- c(80,256) # dimension of tensor
nvars <- prod(dimen) # Number of variables
K <- length(dimen) # Order of tensor

Omega <- array(list(),K)
Omega[[1]] <- TR(dimen[1])
Omega[[2]] <- AR(ratio,dimen[2])
Omega_Tilde <- Matrix.Identifiability(Omega)

cat("First Precision matrix is TR",'\n')
cat("And Second Precision matrix is AR",'\n')

# Separate
Heavy <- FALSE
start <- Sys.time()
C2.Separate <- Separate.simul(n = n, dimen = dimen,Omega = Omega_Tilde,
                              gam = gam, penalty = penalty,
                              norm.off.diag = norm.off.diag,Heavy = Heavy,
                              rate = Separate.rate,sparse = sparse)
end <- Sys.time()
end - start

# Cyclic
Heavy <- FALSE
start <- Sys.time()
C2.Cyclic <- Cyclic.simul(n = n, dimen = dimen, Omega = Omega_Tilde,
                          gam = gam, penalty = penalty,
                          t = 5, rate = Tlasso.rate,
                          sparse = sparse, Heavy = Heavy)
end <- Sys.time()
end - start

# Gemini
Heavy <- FALSE
start <- Sys.time()
C2.Gemini <- Gemini.simul(n = n, dimen = dimen,gam = gam,
                          norm.off.diag = norm.off.diag,
                          penalty = penalty, rate = Gemini.rate,
                          Omega = Omega_Tilde, sparse = sparse,
                          Heavy = Heavy)
end <- Sys.time()
end - start

# Heavy Separate(Each Row)
Heavy = "Matrix"
start <- Sys.time()
C2.Heavy.Matrix.Separate <- Separate.simul(n = n, dimen = dimen,Omega = Omega_Tilde,
                                           gam = gam, penalty = penalty,
                                           norm.off.diag = norm.off.diag,Heavy = Heavy,
                                           rate = Separate.rate,sparse = sparse)
end <- Sys.time()
end - start

# Heavy Separate(Each Matrix)
Heavy = "Constant"
start <- Sys.time()
C2.Heavy.Constant.Separate <- Separate.simul(n = n, dimen = dimen,Omega = Omega_Tilde,
                                             gam = gam, penalty = penalty,
                                             norm.off.diag = norm.off.diag,Heavy = Heavy,
                                             rate = Separate.rate,sparse = sparse)
end <- Sys.time()
end - start
#######################################################################################################################################################################
dimen <- c(32,80) # dimension of tensor
nvars <- prod(dimen) # Number of variables
K <- length(dimen) # Order of tensor

Omega <- array(list(),K)
Omega[[1]] <- TR(dimen[1])
Omega[[2]] <- AR(ratio,dimen[2])
Omega_Tilde <- Matrix.Identifiability(Omega)

# Separate
Heavy <- FALSE
start <- Sys.time()
R1.Separate <- Separate.simul(n = n, dimen = dimen,Omega = Omega_Tilde,
                              gam = gam, penalty = penalty,
                              norm.off.diag = norm.off.diag,Heavy = Heavy,
                              rate = Separate.rate,sparse = sparse)
end <- Sys.time()
end - start

# Cyclic
Heavy <- FALSE
start <- Sys.time()
R1.Cyclic <- Cyclic.simul(n = n, dimen = dimen, Omega = Omega_Tilde,
                          gam = gam, penalty = penalty,
                          t = 5, rate = Tlasso.rate,
                          sparse = sparse, Heavy = Heavy)
end <- Sys.time()
end - start

# Gemini
Heavy <- FALSE
start <- Sys.time()
R1.Gemini <- Gemini.simul(n = n, dimen = dimen,gam = gam,
                          norm.off.diag = norm.off.diag,
                          penalty = penalty, rate = Gemini.rate,
                          Omega = Omega_Tilde, sparse = sparse,
                          Heavy = Heavy)
end <- Sys.time()
end - start

# Heavy Separate(Each Row)
Heavy = "Matrix"
start <- Sys.time()
R1.Heavy.Matrix.Separate <- Separate.simul(n = n, dimen = dimen,Omega = Omega_Tilde,
                                           gam = gam, penalty = penalty,
                                           norm.off.diag = norm.off.diag,Heavy = Heavy,
                                           rate = Separate.rate,sparse = sparse)
end <- Sys.time()
end - start

# Heavy Separate(Each Matrix)
Heavy = "Constant"
start <- Sys.time()
R1.Heavy.Constant.Separate <- Separate.simul(n = n, dimen = dimen,Omega = Omega_Tilde,
                                             gam = gam, penalty = penalty,
                                             norm.off.diag = norm.off.diag,Heavy = Heavy,
                                             rate = Separate.rate,sparse = sparse)
end <- Sys.time()
end - start
###################################################################################################################################
dimen <- c(256,80) # dimension of tensor
nvars <- prod(dimen) # Number of variables
K <- length(dimen) # Order of tensor

Omega <- array(list(),K)
Omega[[1]] <- TR(dimen[1])
Omega[[2]] <- AR(ratio,dimen[2])
Omega_Tilde <- Matrix.Identifiability(Omega)

# Separate
Heavy <- FALSE
start <- Sys.time()
R2.Separate <- Separate.simul(n = n, dimen = dimen,Omega = Omega_Tilde,
                              gam = gam, penalty = penalty,
                              norm.off.diag = norm.off.diag,Heavy = Heavy,
                              rate = Separate.rate,sparse = sparse)
end <- Sys.time()
end - start

# Cyclic
Heavy <- FALSE
start <- Sys.time()
R2.Cyclic <- Cyclic.simul(n = n, dimen = dimen, Omega = Omega_Tilde,
                          gam = gam, penalty = penalty,Heavy = Heavy,
                          t = 5, sparse = sparse, 
                          rate = Tlasso.rate)
end <- Sys.time()
end - start

# Gemini
Heavy <- FALSE
start <- Sys.time()
R2.Gemini <- Gemini.simul(n = n, dimen = dimen,sparse = sparse,
                          gam = gam,Omega = Omega_Tilde,
                          norm.off.diag = norm.off.diag,Heavy = Heavy,
                          penalty = penalty, rate = Gemini.rate)
end <- Sys.time()
end - start

# Heavy Separate(Each Row)
Heavy = "Matrix"
start <- Sys.time()
R2.Heavy.Matrix.Separate <- Separate.simul(n = n, dimen = dimen,Omega = Omega_Tilde,
                                           gam = gam, penalty = penalty,
                                           norm.off.diag = norm.off.diag,Heavy = Heavy,
                                           rate = Separate.rate,sparse = sparse)
end <- Sys.time()
end - start

# Heavy Separate(Each Matrix)
Heavy = "Constant"
start <- Sys.time()
R2.Heavy.Constant.Separate <-  Separate.simul(n = n, dimen = dimen,Omega = Omega_Tilde,
                                              gam = gam, penalty = penalty,
                                              norm.off.diag = norm.off.diag,Heavy = Heavy,
                                              rate = Separate.rate,sparse = sparse)
end <- Sys.time()
end - start

stopCluster(cl)