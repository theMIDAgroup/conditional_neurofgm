rm(list=ls(all=TRUE))
time.start <- Sys.time()

#################################################
## 0. USER DEFINED PARAMETERS (MODIFY THIS PART)
#################################################
score_path = ""
output_path = ""
name_output = ""
func.path = ""
n_basis = 
L = 100
K = 5
thres.ctrl <- c(0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0)
tol.abs =1e-4
tol.rel = 1e-4
eps = 1e-08
verbose = FALSE
##############################################

#################################################
## Include parameter info in the logs
cat("Score source: ", score_path,"\n")
cat("Parameters: \n")
cat("n_basis: ", n_basis ,"\n")
cat("L: ", L ,"\n")
cat("K : ", K  ,"\n")
cat("thres.ctrl : ", thres.ctrl  ,"\n")
#################################################



#############################################
###### 1. Upload all the data needed
#############################################
args <- commandArgs(trailingOnly = TRUE)
j = as.numeric(args[[1]])
scores <- read.csv(score_path)[, -1]
n_nodes <- ncol(scores)/n_basis
n_samples <- nrow(scores)
names <- rep(NA, ncol(scores))
for(l in 1:ncol(scores)){
  names[l] <- paste("f",ceiling(l/n_basis) ,".",l%%n_basis, sep ="")
}
colnames(scores) <- names


#############################################
###### 2.Define functions needed for the computation
#############################################

# Read R files
source(paste(func.path,"ADMM.new.R", sep="/"))      # ADMM Function
source(paste(func.path,"prec.rec.R", sep="/"))      # Precision & Recall
source(paste(func.path,"auc.R", sep="/"))           # Computing AUC from ROC
source(paste(func.path,"FPCA.score.R", sep="/"))    # Computing FPCA scores
source(paste(func.path,"bases.func.R", sep="/"))    # Basis functions for FPCA



##############################################
##### 3.COMPUTATION
##########################################

len.t <- length(thres.ctrl)
n <- nrow(scores)
M <- n_basis
Mp <- ncol(scores)
p <- ceiling(ncol(scores)/M)

# Q: A cosa servono d.array e  d ? Quali sono le dimensioni corrette?
# Per il momento li ho sostuiti com dei vettori di uni della dimensione coretta (?)
d.array <- matrix(1, nrow=p, ncol=(p-1)*M)
d_out <- list()
norm.adj <- rep(NA,p)
for(k in 1:p){
  d_out[[k]] <- d.array[k,]
  norm.adj[k] <- norm(d_out[[k]],"2")
}

P.def <- matrix(0, (p-1)*M, M)
Q.def <- matrix(0.1, (p-1)*M, M)
U.def <- matrix(0.01, (p-1)*M, M)

cat(paste("Processing node ", j,"\n"))
jth.range.y <- (j-1)*M+(1:M)
A.Y <- as.matrix(scores[, jth.range.y])
A.X <- as.matrix(scores[, -jth.range.y])
  
P <- P.def; Q <- Q.def; U <- U.def
  
SCV.mat <- matrix(NA, L , len.t )
  
lambdas <- seq(lambda.sup(A.X, A.Y), 0, length.out=L)
for(l in 1:L){
    lambda <- lambdas[l]
    # Step 1. Use the full dataset to estimate B.hat(lambda)
    if(l%%10 == 0){
      cat(paste("Lambda ",l," of node ",j,"\n",sep=""))}
    
    grp.lasso.result <- tryCatch({
      ADMM.grplasso(A.X = A.X, A.Y = A.Y, d = d_out[[j]],
                               lambda=lambda, rho.init=1,
                               P.in=P, Q.in=Q, U.in=U,
                               tol.abs = tol.abs, tol.rel = tol.rel,
                               maxiter = 400)
    }, error = function(e) {
      message("Error in ADMM.grplasso.two.groups: ", e$message)
      return(NULL) })
    if (is.null(grp.lasso.result)) next
    P <- grp.lasso.result$P
    Q <- grp.lasso.result$Q
    U <- grp.lasso.result$U
    
    
    n_blocks <- nrow(P)/M
    # Process the estimated P into a neighborhood selection vector
    P.frob <- list()
    seq <- 1:p
    seq <- seq[-j]
    for(k in 1:n_blocks){
      key <- seq[k]
      P.frob[[key]] <- norm(P[(k-1)*M + (1:M), ], "F")
    }
    len.t.random <- length(thres.ctrl)
    ind.t <- 1
    for(ind.t in 1:len.t.random){
      threshold <- thres.ctrl[ind.t] * lambda
      # Step 2. finding N.hat.j according to each threhold defined by a combination of lambda  and thres.ctrl
      N.hat.jlt <- rep(FALSE, length(P.frob))
      for(n_block in 1:length(P.frob)){
        if(!(is.null(P.frob[[n_block]]))){
          if(P.frob[[n_block]] > threshold){
            N.hat.jlt[n_block] <- TRUE} 
        }
      }
      if(j == p){N.hat.jlt <- c(N.hat.jlt,FALSE)}
      card.N.hat.jlt <- sum(N.hat.jlt)
      slice.pM <- rep(N.hat.jlt, each=M)
      # Check that takes the right columns
      A.X.eff <- as.matrix(scores[, slice.pM]) # select only the columns of A.X with selected features
      
      SCV.single <- rep(NA, K)
      for(k in 1:K){
        if(card.N.hat.jlt > 0){  # Circumstance A. if |N.hat| >=1
          
          # Step A3: slicing out kth fold as CV set, the rest as training set
          fold.k.ind <- seq(floor((k-1)*n/K) + 1, floor(k*n/K)) # index set for the kth fold of SCV
          fold.k.size <- length(fold.k.ind) # cardinality of Ik under the paper's notation
          A.X.cv <- A.X.eff[fold.k.ind, ]
          A.X.train <- A.X.eff[-fold.k.ind, ]
          A.Y.cv <- A.Y[fold.k.ind, ]
          A.Y.train <- A.Y[-fold.k.ind, ]
          
          # Step A4: calculate B.tilde from a ridge regression
          B.tilde <- solve(t(A.X.train) %*% A.X.train + 0.1*diag(card.N.hat.jlt * M)) %*% t(A.X.train) %*% A.Y.train
          
          # Step A5: evaluate SCV for this single k
          residual <- A.Y.cv - A.X.cv %*% B.tilde # residual is n*M matrix
          emp.cov.residual <- t(residual) %*% residual / fold.k.size # empirical covariance of residual
          #SCV.single[k] <- norm(residual,"F")^2 + log(fold.k.size) * card.N.hat.jlt
          SCV.single[k] <- norm(residual,"F")^2
          
        }else{ # Circumstance B. if cardinality of N.hat is 0
          # Step B3: slicing out kth fold as CV set, the rest as training set
          fold.k.ind <- seq(floor((k-1)*n/K) + 1, floor(k*n/K)) # index set for the kth fold of SCV
          fold.k.size <- length(fold.k.ind) # cardinality of Il under the paper's notation
          A.Y.cv <- A.Y[fold.k.ind, ]
          A.Y.train <- A.Y[-fold.k.ind, ]
          
          # Step B5: evaluate SCV for this single k
          residual <- A.Y.cv # residual is n*M matrix
          emp.cov.residual <- t(residual) %*% residual / fold.k.size # empirical covariance of residual
          SCV.single[k] <- norm(residual,"F")^2
        }
      }# end of the k'th fold
      # Step 6: averaging SCV over K folds for this specific lambda and epsilon
      SCV.mat[l, ind.t] <- mean(SCV.single)
    } # end of all thresholds (epsilon) 1:len.t
  } # end of lambda 1:L
scv.min <- which(SCV.mat == min(SCV.mat), arr.ind=T)
index.optimal <- scv.min[dim(scv.min)[1], ] # there could be multiple. Take the last one
l.optimal <- index.optimal[1]
ind.t.optimal <- index.optimal[2]
lambda.optimal <- lambdas[l.optimal]
t.optimal <- thres.ctrl[ind.t.optimal]
cat("Oprimal lambda: ", lambda.optimal, "\n")
cat("Oprimal thershold: ", t.optimal, "\n")
  
cat("Computing optimal neighboors \n")
grp.lasso.result <- tryCatch({
    ADMM.grplasso(A.X = A.X, A.Y = A.Y, d = d_out[[j]],
                             lambda=lambda.optimal, rho.init=1,
                             P.in=P, Q.in=Q, U.in=U,
                             tol.abs = tol.abs, tol.rel = tol.rel,
                             maxiter = 400)
}, error = function(e) {
message("Error in ADMM.grplasso.two.groups: ", e$message)
return(NULL) })
P <- grp.lasso.result$P
Q <- grp.lasso.result$Q
U <- grp.lasso.result$U
  
  
n_blocks <- nrow(P)/M
# Process the estimated P into a neighborhood selection vector
P.frob <- list()
for(k in 1:n_blocks){
  key <- seq[k]
  P.frob[[key]] <- norm(P[(k-1)*M + (1:M), ], "F")
}
  
N.hat.optimal <- rep(FALSE, length(P.frob))
for(n_block in 1:length(P.frob)){
    if(!(is.null(P.frob[[n_block]]))){
      threshold <- t.optimal*lambda.optimal
      if(P.frob[[n_block]] > threshold){
        N.hat.optimal[n_block] <- TRUE} 
    }
  }
if(j == p){N.hat.optimal <- c(N.hat.optimal,FALSE)}
cat("Optimal neighboors of node ", j, ":\n")
cat(N.hat.optimal)
cat("\n Computational time of: ")
cat( Sys.time() - time.start )
full_result_path = paste(output_path, name_output, "_", j, ".rda", sep ="")
save(N.hat.optimal, file =full_result_path)
cat("Output saved to: ", full_result_path,"\n" )

