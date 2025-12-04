####################################
#  PART 0: PRELIMINARY FUNCTIONS   #
####################################

# Generating M*M Tridiagonal
tridiag <- function(M){
  result <- diag(M)
  for(i in 1:M){
    for(j in 1:M){
      if(abs(i-j)==1) result[i,j] <- 0.5
    }
  }
  return(result)
}

toeplitz <- function(M){
  result <- diag(M)
  for(i in 1:M){
    for(j in 1:M){
      if(i != j){
        result[i,j] <- 0.5 / abs(i-j)
      }
    }
  }
  return(result)
}

# 0.1. A function generating precision matrix of delta
cov.mat.model.A <- function(p,M){
  # Input:
  #   p, number of covariates
  #   M, number of basis functions
  # Output:
  #   the precision matrix that generates delta via MVN
  Theta <- matrix(nrow=p*M, ncol=p*M)
  for(i in 1:p){
    for(j in 1:p){
      if(i==j) Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- toeplitz(M)
      else if(abs(i-j)==1) Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- tridiag(M) * 0.4
      else if(abs(i-j)==2) Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- diag(M) * 0.2
      else Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- matrix(0, M, M)
    }
  }
  eig.val <- eigen(Theta)$values
  if(min(eig.val)<0.05){
    Theta = Theta + (abs(min(eig.val)) + 0.05) * diag(p*M)
  }
  
  return(Theta)
}

cov.mat.model.A.red <- function(p,M,diff_nodes){
  # Input:
  #   p, number of covariates
  #   M, number of basis functions
  # Output:
  #   the precision matrix that generates delta via MVN
  Theta.full <- cov.mat.model.A(p,M)
  Theta.red <- Theta.full
  for(i in 1:diff_nodes){
    for(j in c(i-2, i-1, i+1, i+2)){
      if (j >= 1 && j <= diff_nodes && j != i) {  # Ensure valid index and exclude self
        Theta.red[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)] <- 0
      }
    }
  }
  return(Theta.red)
}