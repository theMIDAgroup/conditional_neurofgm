setwd("/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Graph_estimation")
rm(list=ls(all=TRUE))

library(glasso)
library(glassoFast)
library(sparseMatEst)
library(psych)
library("cglasso")
library("sparsegl")
library(parallel)
library(doParallel)
library(foreach)
library(stabs)

func.path <- "/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Previous_litt/Boxin_Zhao_FGM_Neighboorhood/Functions"
save.path <- "/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Previous_litt/Boxin_Zhao_FGM_Neighboorhood/Simulations_with_other_method_comparison/Model_B"
runtime.path <- "/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Previous_litt/Boxin_Zhao_FGM_Neighboorhood/Simulations_with_other_method_comparison/Model_B"

# Packages
library(fda)
library(matrixcalc)
library(MASS)
library(mvtnorm)

# Global Parameter Settings
p <- 10
mu <- 15 # number of basis used to generate data
M <- 5 # number of basis used for neighborhood selection
n <- 100
tau <- 100 # number of observations
thres.ctrl <- 0 # recognition threshold epsilon_n = thres.ctrl * lambda_n
tol.abs <- 1e-4  # Tolerance (absolute) in ADMM
tol.rel <- 1e-4  # Tolerance (relative) in ADMM
L <- 100 # number of lambdas

# Read R files
source(paste(func.path,"ADMM.new.R", sep="/"))      # ADMM Function
source(paste(func.path,"prec.rec.R", sep="/"))      # Precision & Recall
source(paste(func.path,"auc.R", sep="/"))           # Computing AUC from ROC
source(paste(func.path,"FPCA.score.R", sep="/"))    # Computing FPCA scores
source(paste(func.path,"bases.func.R", sep="/"))    # Basis functions for FPCA
source(paste(func.path,"B.prelim.func.R", sep="/")) # For Model B generation
source(paste(func.path,"ProxAlg_FGM.R", sep="/"))   # FGLasso

total.time.start <- proc.time()


####################################
#     PART 0.1: DATA GENERATION      #
####################################

#    Generate Random Functions     
#      and Observation h_ijk       

# h_ijk = Fourier basis func(t_k) %*% delta_ij + iid error

# 0. Generate precision matrix and real adjacency matrix
set.seed(123)
Theta <- cov.mat.model.B(p, mu) # p*mu by p*mu large square matrix
G.true <- matrix(0, p, p) # p by p adjacency matrix
for(i in 1:p){
  for(j in 1:p){
    if(sum(abs(Theta[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])) > 0)
      G.true[i,j] <- 1
  }
}

# 1. Generating delta
delta <- rmvnorm(n, sigma = solve(Theta))

# 2. Observation time
obs.time <- seq(1/tau, 1, 1/tau) # vector of observation time points of delta

# 3. Fourier basis function for data generation
b.mat.list <- list()
for(j in 1:p){
  b.mat.list[[j]] <- fda.fourier.mat(obs.time, mu)
}

# 4. Observations h_ijk
h <- array(0, c(n, p, tau))
for(i in 1:n){
  for(j in 1:p){
    h[i,j,] <- b.mat.list[[j]] %*% 
      matrix(delta[i, ((j-1)*mu+1) : (j*mu)], ncol=1) + rnorm(tau, 0, 0.5)
  }
}

# Reserved part for PSKL Method
y.list <- list()
for(j in 1:p){
  y.list[[j]] <- h[,j,]
}
####################################
#     PART 0.2: GAIN FPC SCORE       #
####################################

# For the use of gX group
time.start.fpc <- proc.time()
fpc.score <- numeric(0)
for(j in 1:p){
  obs.val.matrix <- matrix(0, nrow=tau, ncol=n)
  for (i in c(1:n)){
    obs.val.vec <- as.vector(h[i, j, ])
    obs.val.matrix[, i] <- obs.val.vec
  }
  bspline.basis <- create.bspline.basis(rangeval=c(0, 1), nbasis=M)
  # Construct a functional data object from the observation
  # bspline basis is the default setting of this function
  # It does not mean that the basis function is bspline!
  fd.object.array <- Data2fd(argvals=obs.time, y=obs.val.matrix, basisobj=bspline.basis)
  # FPCA process
  fpc.score <- cbind(fpc.score, pca.fd(fd.object.array, nharm=M)$scores)
}
time.end.fpc <- proc.time()
runtime.fpc <- (time.end.fpc - time.start.fpc)[3]

lambda.max.gX <- lambda.sup.global.gX(fpc.score, p)
lambdas.gX <- seq(lambda.max.gX, 0, length.out=L)

####################################
#     PART 1: DATA PREPROCESSING  #
####################################

scores <- as.data.frame(fpc.score)
covariates <- data.frame(group= c(rep(0,50),rep(1,50)))
covariates$group <- as.factor(covariates$group)
names <- rep(NA, 50)
for(l in 1:ncol(scores)){
  names[l] <- paste("f",ceiling(l/M) ,".",l%%M, sep ="")
}
colnames(scores) <- names
full_data <- cbind(scores,covariates)
save(full_data,lambdas.gX, G.true, M, p, n, tau, thres.ctrl, tol.abs, tol.rel, L, file = "cFGGM_code_debug_data.RData")


res <- FGGReg_diff_two_groups(scores, # functional score on a defined basis, nrow: subjects; ncol: functions*n_basis.
                                      n_basis = 5, #Number of bases considered 
                                      covariates  = covariates, #Additional covariates to regress on
                                      L = 100, # How many penalization term in the Lasso to try
                                      thres.ctrl = 0, # recognition threshold epsilon_n = thres.ctrl * lambda_n,
                                      verbose = TRUE,
                                      tol.abs =1e-4 ,
                                      tol.rel = 1e-4,
                                      eps = 1e-08)
