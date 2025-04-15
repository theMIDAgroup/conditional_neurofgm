rm(list=ls(all=TRUE))
setwd("/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1")
library(fda)
library(RcppCNPy)

n = 50
tau = 10000
T_end = 100 ##
time <- seq(from = 0, to = T_end, length.out = tau)
p  = 64 ##
M=5

h <- array(0, c(n, p, tau))
for(i in 1:n){
  print( paste('Processing patient ', i))
  signal <- npyLoad(file = paste('Data_fold/fd_data_', as.character(i-1), '.npy', sep=''))
  h[i, ,] <- signal
}


fpc.score <- numeric(0)
for(j in 1:p){
  obs.val.matrix <- matrix(0, nrow=tau, ncol=n)
  for (i in c(1:n)){
    obs.val.vec <- as.vector(h[i, j, ])
    obs.val.matrix[, i] <- obs.val.vec
  }
  bspline.basis <- create.bspline.basis(rangeval=c(0, T_end), nbasis=M)
  fd.object.array <- Data2fd(argvals=time, y=obs.val.matrix, basisobj=bspline.basis)
  # FPCA process
  fpc.score <- cbind(fpc.score, pca.fd(fd.object.array, nharm=M)$scores)
}

write.csv(fpc.score, "fpca_scores.csv")


