setwd("/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies")
rm(list=ls(all=TRUE))

time.start.tot <- Sys.time()

data.path <- "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1"
func.path <- "/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Graph_estimation"
save.path <- "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1"
runtime.path <- "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1"

source(paste(func.path,"cFGGM_functions_two_groups.R", sep="/"))
cat("function source: ", paste(func.path,"cFGGM_functions_two_groups.R", sep="/"),"\n")
scores <- read.csv(paste(data.path, "freq_scores_reord.csv", sep ="/"))[, -1]
cat("score source: ", paste(data.path, "freq_scores_reord.csv", sep ="/"),"\n")
n_basis <- 8
n_nodes <- ncol(scores)/n_basis
n_samples <- nrow(scores)
names <- rep(NA, ncol(scores))
for(l in 1:ncol(scores)){
  names[l] <- paste("f",ceiling(l/n_basis) ,".",l%%n_basis, sep ="")
}
colnames(scores) <- names
covariates <- data.frame(group= c(rep(1,24),rep(0,26)))
covariates$group <- as.factor(covariates$group)
full_data <- cbind(covariates,scores)

L = 100
K = 5
thres.ctrl = c(0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0)
tol.abs =1e-4
tol.rel = 1e-4
eps = 1e-08
verbose = TRUE
cat("Parameters: \n")
cat("n_basis: ", n_basis ,"\n")
cat("L: ", L ,"\n")
cat("K : ", K  ,"\n")
cat("thres.ctrl : ", thres.ctrl  ,"\n")
cat("p_rand: 0.5 \n")

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
for(j in 1:(q + 1)){
  product <- scores * iU[, j]
  if (j != 1) {
      original_colnames <- colnames(product)
      iU_name <- colnames(iU)[j]
      new_colnames <- paste(iU_name, ":", original_colnames, sep = "")
      colnames(product) <- new_colnames
    }
    for (i in 1:p){temp_groups <- c(temp_groups, rep(i+(p*(j-1)),M))}
    interM <- cbind(interM, product)
  }
interM <- interM[, -1]
temp_groups <- temp_groups[-1]
  
# Q: A cosa servono d.array e  d ? Quali sono le dimensioni corrette?
# Per il momento li ho sostuiti com dei vettori di uni della dimensione coretta (?)
d.array <- matrix(1, nrow=p, ncol=(p-1)*M*(q+1))
d_out <- list()
norm.adj <- rep(NA,p)
for(k in 1:p){
    d_out[[k]] <- d.array[k,]
    norm.adj[k] <- norm(d_out[[k]],"2")
  }
  
P.def <- matrix(0, 2*(p-1)*M, M)
Q.def <- matrix(0.1, 2*(p-1)*M, M)
U.def <- matrix(0.01, 2*(p-1)*M, M)

G.mat <- matrix(NA, p, (q+1)*p)


for(j in 1:p){
time.start <- Sys.time()
cat(paste("Processing node ", j,"\n"))
jth.range.y <- (j-1)*M+(1:M)
A.Y <- as.matrix(interM[, jth.range.y])
jth.range.x <- c(jth.range.y,(j+p-1)*M+(1:M))
A.X <- as.matrix(interM[, -jth.range.x])
groups <- temp_groups[-jth.range.x]
                                                                             
P <- P.def; Q <- Q.def; U <- U.def
                                                                             
SCV.mat <- matrix(NA, ceiling(L * 0.5), ceiling(len.t * 0.5))
                                                                             
lambda.max <- lambda.sup(A.X, A.Y)
lambdas <- exp(seq(log(lambda.max), log(1e-4), length.out = L))
                                                                             
set.seed(123+10*j)
random.sel.indexes <- sample(seq_along(lambdas), ceiling(L * 0.5))
random.sel.lambdas <- lambdas[random.sel.indexes[order(random.sel.indexes)]]
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
        thresholds <- lambda * thres.ctrl
        set.seed(234+10*j)
        random.sel.thresholds <- thresholds[sample(seq_along(thresholds), ceiling(len.t * 0.5))]
        len.t.random <- length(random.sel.thresholds)
                                                                               
        for(ind.t in 1:len.t.random){
         threshold <- random.sel.thresholds[ind.t]
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
           SCV.single[k] <- norm(residual,"F")^2 + log(fold.k.size) * card.N.hat.jlt
                                                                                     
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
   t.optimal <- thres.ctrl[ind.t.optimal]
                                                                             
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
         if(P.frob[[n_block]] > t.optimal){
               N.hat.optimal[n_block] <- TRUE} 
            }
     }
      if(j == p){N.hat.optimal <- c(N.hat.optimal,FALSE)}
      cat("Optimal neighboors of node ", j, ":\n")
      cat(N.hat.optimal)
      cat("\n Computational time of: ")
      cat( Sys.time() - time.start )
      G.mat[j,] <- N.hat.optimal 
      cat("\n")
}

save(G.mat, file = paste(save.path,"res_1404_sequential.RData", sep="/"))
cat("Output saved to: ",paste(save.path,"res_1404_sequential.RData", sep="/"), "\n" )
Sys.time() - time.start.tot
