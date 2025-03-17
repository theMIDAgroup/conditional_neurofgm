setwd("/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Graph_estimation")
rm(list=ls(all=TRUE))
func.path <- "/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Previous_litt/Boxin_Zhao_FGM_Neighboorhood/Functions"
source(paste(func.path,"ADMM.new.R", sep="/"))      # ADMM Function
load("cFGGM_code_debug_data.RData")

covariates <- data.frame(group = full_data[, "group"])
scores <- full_data[, -ncol(full_data)]
n_basis <- M
scr = TRUE
gamma = NULL
verbose=TRUE

n <- nrow(scores)
M <- n_basis
Mp <- ncol(scores)
p <- ceiling(ncol(scores)/M)
Ip <- diag(rep(1, p))

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

if (n < (q+1)*Mp) {
  warning("The sample size is too small! Network estimate may be unreliable!")
}


## INITIAL SCREENING

if (scr) {
  if (verbose) {
    cat("Computing correlation matrix \n")
  }
  score_matrix <- as.matrix(scores)
  Scor <- cor(score_matrix)
  if (verbose) {
    cat("Done computing correlation matrix \n")
  }
  scr_index <- matrix(1, Mp, Mp)
  if (is.null(gamma)) {
    gamma <- quantile(abs(Scor), probs = 0.1)
  }
  for(k in 1:p){
    row_range <- ((M * (k - 1)) + 1):(M * k)
    scr_index[row_range, row_range] <- 0
    for (l in 1:p){
      col_range <- ((M * (l - 1)) + 1):(M * l)
      if (all(abs(Scor)[row_range, col_range] <= gamma))
        scr_index[row_range, col_range] <- 0}
  }
  colnames(scr_index) <- colnames(scores)
  rownames(scr_index) <- colnames(scores)
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

################################################
## PART 2. COMPUTE THE MULTIPLE REGRESSION
################################################

Delta_scores <- matrix(0, Mp, ncol(interM))
colnames(Delta_scores) <- colnames(interM)
rownames(Delta_scores) <- colnames(scores)
No_sim_Delta_hat <- matrix(0, Mp, ncol(interM))
colnames(No_sim_Delta_hat) <- colnames(interM)
rownames(No_sim_Delta_hat) <- colnames(scores)


time.start.gX <- proc.time()

G.list.gX <- list()
TPR.and.gX <- rep(0,L)
FPR.and.gX <- rep(0,L)
TPR.or.gX <- rep(0,L)
FPR.or.gX <- rep(0,L)

# Q: A cosa servono d.array e  d ? Quali sono le dimensioni corrette?
# Per il momento li ho sostuiti com dei vettori di uni della dimensione coretta (?)
d.array <- matrix(1, nrow=p, ncol=p*M*(q+1))

d_out <- list()
norm.adj <- rep(NA,p)
for(k in 1:p){
  d_out[[k]] <- d.array[k,]
  norm.adj[k] <- norm(d_out[[k]],"2")
}
grp.soft.thres <- function(V.Grp,M, d, lambda, rho){
  # V.Grp: a group of p matrices (Vk), each dim M * M
  # d: vector, length pM.
  # lambda: penalty param of group lasso
  # rho: penalty param of augmented Lagrangian
  
  p <- length(V.Grp)
  D.Grp <- list()
  P.Grp <- list() # output group of p matrices (Pk), each dim M * n
  for(k in 1:p){
    if (!is.null(V.Grp[[k]])){
      D.Grp[[k]] <- diag(d[(k-1)*M+(1:M)])
      left.plus <- matrix(0, M, M)
      norm.F <- norm(D.Grp[[k]] %*% V.Grp[[k]], "F")
      l=2
      for(l in 1:M){
        left.plus[l,l] <- max(0, 1 - lambda / rho * D.Grp[[k]][l,l]^2 / norm.F )
      }
      P.Grp[[k]] <- left.plus %*% V.Grp[[k]]
    } else {
      P.Grp[[k]] <- NULL
    }
    
  }
  
  return(P.Grp)
}




j=1
V.array <- foreach(j = 1:p, .combine="rbind") %dopar% {
  jth.range <- (((j-1)*M+1) : (j*M)) + 1
  Y <- as.matrix(interM[, jth.range])
  indexes_to_select <- rep(TRUE, ncol(interM ))
  if (scr) {
    zero_col_indices <- which(scr_index[(j*M), ] == 0)
    zero_col_names <- colnames(scr_index)[zero_col_indices]
    for (name in zero_col_names) {
      pattern <- paste0(":", name, "\\b")
      indexes_to_select[grepl(pattern, colnames(interM))] <- FALSE
      indexes_to_select[colnames(interM) == name] <- FALSE } 
  }else{
    protnames <- colnames(interM[jth.range])
    indexes_to_select[jth.range] <- FALSE
    for (name in protnames){
      pattern <- paste0(":", name, "\\b")
      indexes_to_select[grepl(pattern, colnames(interM))] <- FALSE }
  }
  indexes_to_select_reg <- indexes_to_select[-1]
  Xmat <- interM[,-1]
  Xmat[, !indexes_to_select_reg] <- 0
  Xmat <- as.matrix(Xmat)
  groups <- temp_groups[-1]
  groups[!indexes_to_select_reg] <- -1

  # default warm start PQU
  P.def <- matrix(0, ncol(Xmat), M)
  Q.def <- matrix(0.1, ncol(Xmat), M)
  # P e Q sono le matrici dei coefficienti delle regressioni, sulle righe ci sono le covaraite e le colonne sono legate a ogni 
  # componenete della Y
  U.def <- matrix(0.01, ncol(Xmat), M)
  # differenza tra P e Q
  
  P <- P.def; Q <- Q.def; U <- U.def
  
  V.j <- matrix(NA, L, p)
  l =99
  for(l in 1:L){
    lambda <- lambdas.gX[l]
    grp.lasso.result <- ADMM.grplasso(A.X = Xmat, A.Y = Y, groups=groups, d = d_out[[j]],
                                      lambda=lambda, rho.init=1,
                                      P.in=P, Q.in=Q, U.in=U,
                                      tol.abs = tol.abs, tol.rel = tol.rel,
                                      maxiter = 100)
    # Modifica della funzione ADMM.grplasso
    
    ADMM.grplasso <- function(A.X, A.Y,groups, d, lambda,
                              rho.init=1, maxiter=2000, mu=10, tau=2,
                              P.in, Q.in, U.in, tol.rel=1e-4, tol.abs=1e-4){
    n <- nrow(A.Y)
    M <- ncol(A.Y)
    p_original <- ncol(A.X) / M # corresponds to p-1 in theory
    p <- ncol(A.X)
    
    # 1.2. Recover matrices used in computation
    D <- diag(d)
    
    AX.D <- A.X %*% D #A.X
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
    
    
    V <- Q.old - U.old # pM * M
    V.Grp <- list()
    # NEW:DEFINE THE GROUPS
    for(k in unique(groups)){
      if(!k==-1){
        V.Grp[[k]] <- V[groups==k, ] # each Vk dim M * M
      }
    }
    
    # 1.5.2. indicator of rho's change
    rho <- rho.init
    iter =1
    for(iter in 1:maxiter){
      
      if(iter==1){
        rho.changed <- T # For iteration 1, we need it True for Q update
      }
      
      #############################
      ## 7.1. Updating Variables ##
      #############################
      # 7.1.1. Update P with group soft thresholding
      P.new <- matrix(NA, nrow=ncol(A.X), ncol=M)
      P.Grp <- grp.soft.thres(V.Grp,M, d, lambda, rho) # Group of p, each M * M
      for(k in unique(groups)){
        if(!k==-1){
          P.new[groups==k, ]<- P.Grp[[k]] # each Vk dim M * M
        } else {
          P.new[groups==k, ]<- 0 
        }
      }
      
      # 7.1.2. Update Q
      if(rho.changed){
        # if rho not changed, this term stays the same
        left.inv <- solve(DXXD.n + rho * diag(p))
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
      tol.prim <- tol.abs * sqrt( p_original*M*M) + tol.rel * max(norm(P.new, "F"), 
                                                        norm(Q.new, "F"))
      tol.dual <- tol.abs * sqrt( p_original*M*M) + tol.rel * norm(U.new, "F")
      # sqrt(pM^2) is because the Frob norms are in R^{pM * M}.
      
      # 7.2.3. Compare residuals with tolerance
      if(prim.resid < tol.prim & dual.resid < tol.dual){
        break
      }
      
      #####################################
      ## 7.3. Prepare for next iteration ##
      #####################################
      
      # 7.3.1. Update V
      V <- Q.new - U.new # pM * M
      for(k in unique(groups)){
        if(!k==-1){
          V.Grp[[k]] <- V[groups==k, ] # each Vk dim M * M
        }
      }
      
      # 7.3.2. Update rho for next round
      if(prim.resid > mu * dual.resid){
        rho <- tau * rho
        rho.changed <- T
      }else if(dual.resid > mu * prim.resid){
        rho <- rho / tau
        rho.changed <- T
      }else{
        rho.changed <- F
      }
      rho.path <- c(rho.path, rho)
      
      # 7.3.3. old <- new, end of one iterate
      P.old <- P.new
      Q.old <- Q.new
      U.old <- U.new
      
      
      # 7.4.1. Record objective function
      obj <- objective.function(A.X=A.X, A.Y=A.Y, d=d, B=P.new, lambda=lambda)
      obj.func.path <- c(obj.func.path, obj)
    }
    
    cat(paste("ADMM converges after ", iter," iterations.\n", sep=""))
    
    result <- list(P=P.new, Q=Q.new, U=U.new,
                   prim.res=prim.resid.vec, dual.res=dual.resid.vec,
                   rho.path=rho.path,
                   tol.prim.path=tol.prim.path, tol.dual.path=tol.dual.path,
                   obj.func.path=obj.func.path)
    return(result)
    }
    
    
    
    
    P <- grp.lasso.result$P
    Q <- grp.lasso.result$Q
    U <- grp.lasso.result$U
    
    # Process the estimated P into a neighborhood selection vector
    P.frob <- rep(0, p-1)
    for(k in 1:(p-1)){
      P.frob[k] <- norm(P[(k-1)*M + (1:M), ], "F")
    }
    
    threshold <- lambda * thres.ctrl
    
    # Neighbor recognization
    V.jl <- rep(0, p)
    for(juliet in 1:(p-1)){
      if(juliet < j){
        if(P.frob[juliet] > threshold)
          V.jl[juliet] <- 1
      }else{
        if(P.frob[juliet] > threshold)
          V.jl[juliet + 1] <- 1
      }
    }
    
    V.j[l,] <- V.jl
  }
  V.j
}


# From v.array to G.list[[1 to L]]
G.list.gX <- list()
for(l in 1:L){
  G.list.gX[[l]] <- matrix(NA, p, p)
  for(j in 1:p){
    G.list.gX[[l]][j,] <- V.array[(j-1)*L+l,]
  }
}


# Calculate TPR and FPR
for(l in 1:L){
  TPR.and.gX[l] <- prec.rec(G.true, G.list.gX[[l]], type="AND")$TPR
  FPR.and.gX[l] <- prec.rec(G.true, G.list.gX[[l]], type="AND")$FPR
  TPR.or.gX[l] <-  prec.rec(G.true, G.list.gX[[l]], type="OR" )$TPR
  FPR.or.gX[l] <-  prec.rec(G.true, G.list.gX[[l]], type="OR" )$FPR
}

time.end.gX <- proc.time()
runtime.gX <- (time.end.gX - time.start.gX)[3]


