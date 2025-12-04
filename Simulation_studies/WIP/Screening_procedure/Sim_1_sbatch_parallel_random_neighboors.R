rm(list=ls(all=TRUE))
time.start <- Sys.time()

#################################################
## 0. USER DEFINED PARAMETERS (MODIFY THIS PART)
#################################################
score_path = "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Model A/fpc_score_25_balanced.csv"
grouping_path = "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Model A/true_grouping_reverse_25.csv"
#distance_path =""
output_path = "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Model A/sample_size_study/Screening_procedure/25_balanced_rand_k15/"
name_output = "rand_hyper_search"
n_basis = 5
L = 100
K = 5
thres.ctrl = c(0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0)
tol.abs =1e-4
tol.rel = 1e-4
eps = 1e-08
verbose = FALSE
p.rand.lam = 0.5
p.rand.thr = 1
pre_screen_corr = FALSE
pre_screen_corr_threshold = NULL
pre_screen_dist = FALSE
pre_screen_dist_neigh = 15
pre_screen_random = TRUE
pre_screen_random_neigh = 15
##############################################

#################################################
## Include parameter info in the logs
cat("Score source: ", score_path,"\n")
cat("Grouping source: ", grouping_path,"\n")
cat("Parameters: \n")
cat("n_basis: ", n_basis ,"\n")
cat("L: ", L ,"\n")
cat("K : ", K  ,"\n")
cat("thres.ctrl : ", thres.ctrl  ,"\n")
cat("p.rand.lam:", p.rand.lam, "\n")
cat("p.rand.thr:", p.rand.thr, "\n")
if (pre_screen_corr == TRUE) {cat("The pre-screening procedure based on correlation is active \n")}
if(pre_screen_dist == TRUE) {cat("The pre-screening procedure based on distance is active \n")}
if(pre_screen_random == TRUE) {cat("The pre-screening procedure based on random is active \n")}
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
covariates <- data.frame( group = as.factor(read.csv(grouping_path)[,-1]))
full_data <- cbind(covariates,scores)
if(pre_screen_dist == TRUE){
  #dist_matrix <- read.csv(dist_path)[, -1]
  dist_matrix <- outer(1:50, 1:50, FUN = function(x, y) abs(x - y))
}

#############################################
###### 2.Define functions needed for the computation
#############################################

###### SAME lambda.sup
lambda.sup <- function(A.X, A.Y){
  n <- nrow(A.X)
  M <- ncol(A.Y)
  p <- ncol(A.X)/M
  candidates <- rep(0,p)
  for(j in 1:p){
    A.X.j <- A.X[,(j-1)*M + (1:M)]
    candidates[j] <- norm(t(A.X.j) %*% A.Y,"F") / n
  }
  return(max(candidates))
}


###### UPDATE grp.soft.thres
grp.soft.thres.two.groups <- function(V.Grp, d, lambda, rho){
  # V.Grp: a group of p matrices (Vk), each dim M * M
  # d: vector, length pM.
  # lambda: penalty param of group lasso
  # rho: penalty param of augmented Lagrangian
  M <- nrow(V.Grp[[1]])
  n_blocks <- length(V.Grp)
  D.Grp <- list()
  P.Grp <- list() # output group of p matrices (Pk), each dim M * n
  
  for(k in 1:n_blocks){
    D.Grp[[k]] <- diag(d[(k-1)*M+(1:M)])
    left.plus <- matrix(0, M, M)
    norm.F <- norm(D.Grp[[k]] %*% V.Grp[[k]], "F")
    for(l in 1:M){
      left.plus[l,l] <- max(0, 1 - lambda / rho * D.Grp[[k]][l,l]^2 / norm.F )
    }
    
    P.Grp[[k]] <- left.plus %*% V.Grp[[k]]
  }
  
  return(P.Grp)
}

######### UPDATE objective.function
objective.function.two.groups <- function(A.X, A.Y, d, B, lambda){
  # A.X: matrix, n * M*(q+1)*num_nodes_included_j
  # A.Y: matrix, n * M
  # d: vector, length M*(q+1)*num_nodes_included_j
  # B: matrix, M*(q+1)*num_nodes_included_j * M
  # Q: matrix, M*(q+1)*num_nodes_included_j * M
  # lambda: penalty parameter
  
  n <- nrow(A.X)
  M <- ncol(A.Y)
  n_blocks <- ncol(A.X) / M
  
  D <- diag(d)
  term.1 <- norm(A.X %*% D %*% B - A.Y, "F")^2 / (2*n)
  term.2 <- 0
  for(k in 1:n_blocks){
    B.k <- B[(k-1)*M + (1:M), ]
    Dkk <- diag(d[(k-1)*M + (1:M)])
    term.2 <- term.2 + lambda * norm(Dkk %*% B.k ,"F")
  }
  obj <- term.1 + term.2
  return(obj)
}

########### UPDATE ADMM.grplasso -> miss a parantesis or one more

ADMM.grplasso.two.groups <- function(A.X, A.Y, d, lambda,
                                     rho.init=1, maxiter=2000, mu=10, tau=2,
                                     P.in, Q.in, U.in, tol.rel=1e-4, tol.abs=1e-4){
  
  # Input:
  #   A.X, n*2(p-1)M matrix.
  #   A.Y, n*M matrix
  #   d, a 2(p-1)M long vector
  # 1.1. Reading in the parameters
  n <- nrow(A.Y)
  M <- ncol(A.Y)
  n_blocks <- ncol(A.X) / M # corresponds to 2(p-1) in theory
  
  # 1.2. Recover matrices used in computation
  D <- diag(d)
  AX.D <- A.X %*% D
  DXXD.n <- (t(AX.D) %*% AX.D) / n
  DXY.n <- (t(AX.D) %*% A.Y)/n
  
  # 1.3. Setting initial values for P,Q,U
  P.old <- P.in
  Q.old <- Q.in
  U.old <- U.in
  
  # 1.4. Initialize output quantities
  prim.resid.vec <- numeric(0)
  dual.resid.vec <- numeric(0)
  rho.path <- rho.init
  tol.prim.path <- numeric(0)
  tol.dual.path <- numeric(0)
  obj.func.path <- numeric(0)
  
  # 1.5. Cached quantities -- initial
  # 1.5.1. Initial Vk's
  V <- Q.old - U.old # pM * M
  V.Grp <- list()
  for(k in 1:n_blocks){
    V.Grp[[k]] <- V[(k-1)*M + (1:M), ] # each Vk dim M * M
  }
  
  # 1.5.2. indicator of rho's change
  rho <- rho.init
  rho.changed <- TRUE
  ##############################
  ## 2. Main part: Iterations ##
  ##############################
  for(iter in 1:maxiter){
    
    #############################
    ## 7.1. Updating Variables ##
    #############################
    # 7.1.1. Update P with group soft thresholding
    P.new <- matrix(NA, nrow=n_blocks*M, ncol=M)
    P.Grp <- grp.soft.thres.two.groups(V.Grp, d, lambda, rho) # Group of p, each M * M
    for(k in 1:n_blocks){
      P.new[(k-1)*M + (1:M), ] <- P.Grp[[k]]
    }
    
    # 7.1.2. Update Q
    if(rho.changed){
      left.inv <- tryCatch(
        solve(DXXD.n + rho * diag(nrow(DXXD.n))),
        error = function(e) {
          stop("Matrix inversion error: Check dimensions of DXXD.n and rho term")
        }
      )
    }
    Q.new <- left.inv %*% (DXY.n + rho * P.new + rho * U.old)
    
    # 7.1.3. Update U
    U.new <- U.old + P.new - Q.new
    
    
    #####################################
    ## 7.2. Checking Stopping Criteria ##
    #####################################
    
    # 7.2.1. Calculate primal and dual residuals
    prim.resid <- norm( (P.new - Q.new) , "F")
    dual.resid <- norm( rho * (Q.new - Q.old) , "F")
    prim.resid.vec <- c(prim.resid.vec,prim.resid)
    dual.resid.vec <- c(dual.resid.vec,dual.resid)
    
    # 7.2.2. Calculate tolerance
    tol.prim <- tol.abs * sqrt(n_blocks*M*M) + tol.rel * max(norm(P.new, "F"), 
                                                             norm(Q.new, "F"))
    tol.dual <- tol.abs * sqrt(n_blocks*M*M) + tol.rel * norm(U.new, "F")
    # sqrt(n_blocksM^2) is because the Frob norms are in R^{n_blocksM * M}.
    
    # 7.2.3. Compare residuals with tolerance
    if(prim.resid < tol.prim && dual.resid < tol.dual){
      break
    }
    
    #####################################
    ## 7.3. Prepare for next iteration ##
    #####################################
    
    # 7.3.1. Update V
    V <- Q.new - U.new # pM * M
    V.Grp <- list()
    for(k in 1:n_blocks){
      V.Grp[[k]] <- V[((k-1)*M + 1):(k*M), ] # each Vk dim M * n
    }
    
    # 7.3.2. Update rho for next round
    if(prim.resid > mu * dual.resid){
      rho <- tau * rho
      rho.changed <- TRUE
    }else if(dual.resid > mu * prim.resid){
      rho <- rho / tau
      rho.changed <- TRUE
    }else{
      rho.changed <- FALSE
    }
    rho.path <- c(rho.path, rho)
    
    # 7.3.3. old <- new, end of one iterate
    P.old <- P.new
    Q.old <- Q.new
    U.old <- U.new
    
    obj <- tryCatch(
      objective.function.two.groups(A.X=A.X, A.Y=A.Y, d=d, B=P.new, lambda=lambda),
      error = function(e) {
        stop("Error in objective function calculation: Check objective.function.two.groups()")
      }
    )
    
    obj.func.path <- c(obj.func.path, obj)
  }
  
  if(iter == maxiter){
      cat("ADMM did not converge.")
  }
  
  
  result <- list(P=P.new, Q=Q.new, U=U.new,
                 prim.res=prim.resid.vec, dual.res=dual.resid.vec,
                 rho.path=rho.path,
                 tol.prim.path=tol.prim.path, tol.dual.path=tol.dual.path,
                 obj.func.path=obj.func.path)
  return(result)}

##############################################
##### 3.COMPUTATION
##########################################

len.t <- length(thres.ctrl)
n <- nrow(scores)
M <- n_basis
Mp <- ncol(scores)
p <- ceiling(ncol(scores)/M)
if (is.null(covariates)) {
  C <- data.frame(Zeros = rep(0, n))
  iU <-data.frame(rep(1, n))
  colnames(iU) <- "(Intercept)"
  q <- 0
} else {
  numeric_columns <- sapply(covariates, is.numeric)
  covariates[, numeric_columns] <- scale(covariates[, numeric_columns])
  C <- model.matrix(~ ., covariates)
  q <- ncol(C) - 1
  iU <- C
}

## DEFINITION OF THE DESIGN MATRIX

interM <- data.frame(Intercept = rep(1, n))
temp_groups <- c(0)

if (verbose) {
  cat("Comuputation of the design matrix \n ")
}
for(c in 1:(q + 1)){
  product <- scores * iU[, c]
  if (c != 1) {
    original_colnames <- colnames(product)
    iU_name <- colnames(iU)[c]
    new_colnames <- paste(iU_name, ":", original_colnames, sep = "")
    colnames(product) <- new_colnames
  }
  for (i in 1:p){temp_groups <- c(temp_groups, rep(i+(p*(c-1)),M))}
  interM <- cbind(interM, product)
}
interM <- interM[, -1]
temp_groups <- temp_groups[-1]

cat(paste("Processing node ", j,"\n"))
jth.range.y <- (j-1)*M+(1:M)

############
# INCLUDE THE PRE_SCREENING
###########
scr_index <- rep(FALSE, p)
if (pre_screen_corr) {
  if (verbose) { cat("Computing correlation matrix \n") }
  scoredata_matrix <- as.matrix(scores)
  Scor <- cor(scoredata_matrix)
  if (verbose) { cat("Done computing correlation matrix \n") }
  if (is.null(pre_screen_corr_threshold)) { pre_screen_corr_threshold <- quantile(abs(Scor), probs = 0.2) }
  for(it in 1:p){
      if(it != j && sum(abs(Scor[(j-1)*M + (1:M), (it-1)*M + (1:M)])>pre_screen_corr_threshold)>0){ scr_index[it] <- TRUE}
  }
} else if(pre_screen_dist){
  dist_j <- as.numeric(dist_matrix[j,])
  dist_j[j] <- Inf
  nearest_indices <- unique(dist_j)[order(unique(dist_j))][1:pre_screen_dist_neigh]
  scr_index <- dist_j %in% nearest_indices
} else if(pre_screen_random){
  set.seed(243+10*j)
  rand_neig <- sample(1:p,pre_screen_random_neigh)
  scr_index[rand_neig] <- TRUE
  scr_index[j] <- FALSE
} else {
  scr_index <- rep(TRUE, p)
  scr_index[j] <- FALSE
}

num_nodes_included_j <- sum(scr_index)
num_cov_j <- num_nodes_included_j * M * (q+1)

#########
#### CHECK
#######
# Q: A cosa servono d.array e  d ? Quali sono le dimensioni corrette?
# Per il momento li ho sostuiti com dei vettori di uni della dimensione coretta (?)
d.array <- matrix(1, nrow=p, ncol=num_cov_j)
d_out <- list()
norm.adj <- rep(NA,p)
for(k in 1:p){
  d_out[[k]] <- d.array[k,]
  norm.adj[k] <- norm(d_out[[k]],"2")
}

P.def <- matrix(0, num_cov_j, M)
Q.def <- matrix(0.1, num_cov_j, M)
U.def <- matrix(0.01, num_cov_j, M)

A.Y <- as.matrix(interM[, jth.range.y])
jth.range.x <- rep(scr_index, each = M)
jth.range.x <- rep(jth.range.x, times = (q+1))
A.X <- as.matrix(interM[, jth.range.x])
groups <- temp_groups[jth.range.x]
  
P <- P.def; Q <- Q.def; U <- U.def
  
SCV.mat <- matrix(NA, ceiling(L * p.rand.lam), ceiling(len.t * p.rand.thr))
#SCV.mat <- matrix(NA, L , len.t )
  
lambda.max <- lambda.sup(A.X, A.Y)
lambdas <- exp(seq(log(lambda.max), log(1e-4), length.out = L))

if(p.rand.lam == 1){
  random.sel.lambdas <- lambdas
} else{
set.seed(123+10*j)
random.sel.indexes <- sample(seq_along(lambdas), ceiling(L * p.rand.lam))
random.sel.lambdas <- lambdas[random.sel.indexes[order(random.sel.indexes)]]
}

L_random <- length(random.sel.lambdas)
  
for(l in 1:L_random){
    lambda <- random.sel.lambdas[l]
    # Step 1. Use the full dataset to estimate B.hat(lambda)
    if(l%%10 == 0){
      cat(paste("Lambda ",l," of node ",j,"\n",sep=""))}
    
    grp.lasso.result <- tryCatch({
      ADMM.grplasso.two.groups(A.X = A.X, A.Y = A.Y, d = d_out[[j]],
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
    for(k in 1:n_blocks){
      key <- unique(groups[(k-1)*M + (1:M)])
      if(length(key) > 1){
        message("Error in group definition for block ", k)
        next} else {
          P.frob[[key]] <- norm(P[(k-1)*M + (1:M), ], "F")}
    }
    if(p.rand.thr == 1){
      random.sel.thresholds <- thres.ctrl
    }else{
      set.seed(234+10*j)
      random.sel.thresholds <- thres.ctrl[sample(seq_along(thres.ctrl), ceiling(len.t * p.rand.thr))]
    }
    len.t.random <- length(random.sel.thresholds)
    
    for(ind.t in 1:len.t.random){
      threshold <- lambda * random.sel.thresholds[ind.t]
      # Step 2. finding N.hat.j according to each threhold defined by a combination of lambda  and thres.ctrl
      N.hat.jlt <- rep(FALSE, length(P.frob))
      for(n_block in 1:length(P.frob)){
        if(!(is.null(P.frob[[n_block]]))){
          if(P.frob[[n_block]] > threshold){
            N.hat.jlt[n_block] <- TRUE} 
        }
      }
      if(length(N.hat.jlt) < max(temp_groups)){N.hat.jlt <- c(N.hat.jlt,rep(FALSE,max(temp_groups)-length(N.hat.jlt)))}
      card.N.hat.jlt <- sum(N.hat.jlt)
      slice.pM <- rep(N.hat.jlt, each=M)
      # Check that takes the right columns
      A.X.eff <- as.matrix(interM[, slice.pM]) # select only the columns of A.X with selected features
      
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
lambda.optimal <- random.sel.lambdas[l.optimal]
t.optimal <- random.sel.thresholds[ind.t.optimal]
cat("Oprimal lambda: ", lambda.optimal, "\n")
cat("Oprimal thershold: ", t.optimal, "\n")
  
cat("Computing optimal neighboors \n")
grp.lasso.result <- tryCatch({
    ADMM.grplasso.two.groups(A.X = A.X, A.Y = A.Y, d = d_out[[j]],
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
    key <- unique(groups[(k-1)*M + (1:M)])
    if(length(key) > 1){
      message("Error in group definition for block ", k)
      next} else {
        P.frob[[key]] <- norm(P[(k-1)*M + (1:M), ], "F")
      }
  }
  
N.hat.optimal <- rep(FALSE, length(P.frob))
for(n_block in 1:length(P.frob)){
    if(!(is.null(P.frob[[n_block]]))){
      threshold <- t.optimal*lambda.optimal
      if(P.frob[[n_block]] > threshold){
        N.hat.optimal[n_block] <- TRUE} 
    }
  }
if(length(N.hat.optimal) < max(temp_groups)){N.hat.optimal <- c(N.hat.optimal,rep(FALSE,max(temp_groups)-length(N.hat.optimal)))}
cat("Optimal neighboors of node ", j, ":\n")
cat(N.hat.optimal)
cat("\n Computational time of: ")
cat( Sys.time() - time.start )
full_result_path = paste(output_path, name_output, "_", j, ".rda", sep ="")
save(N.hat.optimal, file =full_result_path)
cat("Output saved to: ", full_result_path,"\n" )

