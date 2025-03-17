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

# To be adapted to ours
CV_hyperparam <-
  function(x = NULL, # expression data, nrow: subjects; ncol: features.
           known_ppi = NULL, # previous known PPI
           covariates = NULL,       #covariates to regress on
           lambda = NULL,
           nlambda = 5,
           lambda.min.ratio = 0.1, # ratio for minial lambda in a grid search
           scr = TRUE, # optional screening to speed up
           gamma = NULL,
           rho = NULL,
           weight = 1.1,
           eta = 0,
           verbose = FALSE,
           eps = 1e-08,
           seed = NULL, # random seed for K fold sample split
           K = 5, # K fold to optimize lambda
           crit.cv = c("PBIC", "BIC", "loglik", "ploglik", "AIC",  "EBIC"),
           trace = c("progress", "print", "none"),
           ...) {
    n <- nrow(x)
    p <- ncol(x)
    

    # Initial checks on the covariates
    if (is.null(covariates)) {
      covariates = data.frame(Zeros = rep(0,n))
    }

    X = scale(x, center = TRUE, scale = TRUE)
    if(verbose){
      cat("Computing correlation matrix \n ")
    }
    S = cor(X)
    if(verbose){
      cat("Done with the matrix computation \n")
    }

    # Definition of the lambda sequence
    if (!is.null(lambda)){
      nlambda = length(lambda)}
    if (is.null(lambda)) {
      if (is.null(nlambda)) {
        nlambda = 15 }
      if (is.null(lambda.min.ratio)) {
        lambda.min.ratio = 0.1 }
      lambda.max = max(max(S - diag(p)), -min(S - diag(p)))
      lambda.min = lambda.min.ratio * lambda.max
      lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
      rm(lambda.max, lambda.min, lambda.min.ratio)
      gc()
    }
  
    
    #Match values
    crit.cv = match.arg(crit.cv)
    trace = match.arg(trace)
    lambda = sort(lambda)

    #Initialize the metrics
    CV_nLogLikelihood = CV_pNLogLikelihood = CV_Pbic.score = matrix(0, length(lambda), K)
    CV_ebic.score = CV_aic.score = CV_bic.score =matrix(0, length(lambda), K)
    #CV_density = CV_nLogLikelihood
    
    # set progress bar
    if (trace == "progress") {
      progress = txtProgressBar(max = K * length(lambda), style = 3)
    }
    # No need to create folds if K = 1
    if (K == 1) {
      # set sample size
      X.train = X
      empcov.train = S
      t = 1
      # loop over all tuning parameters
      for (loop.lambda in 1:length(lambda)) {
        # compute the penalized likelihood precision matrix
        # estimator
        fit = GGM_Estimation( x=X.train,
          known_ppi=known_ppi,
          covariates=covariates,
          scr=scr,
          gamma=gamma,
          lambda=lambda[loop.lambda],
          rho=rho,
          weight=weight,
          eta=eta,
          verbose=verbose,
          eps = eps
        )
        
        
        siginv = fit$invcov
        no.edge = (sum(abs(siginv) > eps)) / 2
        #CV_density[loop.lambda, t] = sparsePrec(fit$Adj)
        #Compute the observed negative loglikelihood
        CV_nLogLikelihood[loop.lambda, t] = tr(empcov.train %*% siginv) - 
                                  determinant(siginv, logarithm = T)$modulus
        
        #Compute the observed penalized negative loglikelihood
        CV_pNLogLikelihood[loop.lambda, t] = CV_nLogLikelihood[loop.lambda, t] + fit$lambda *
          sum(abs((1 - diag(p)) * siginv))
        
        #Compute the different metrices
        CV_bic.score[loop.lambda, t] = tr(empcov.train %*% siginv) -
          determinant(siginv, logarithm = T)$modulus + log(n) * no.edge / n
        CV_Pbic.score[loop.lambda, t] = CV_pNLogLikelihood[loop.lambda, t] + log(n) * no.edge / n
        CV_ebic.score[loop.lambda, t] = CV_bic.score[loop.lambda, t] + 4 * log(p) * no.edge / n
        CV_aic.score[loop.lambda, t] = CV_nLogLikelihood[loop.lambda, t] + 2 * no.edge / n
        
        # update progress bar
        if (trace == "progress") {
          setTxtProgressBar(progress, loop.lambda + (t - 1) *
                              length(lambda))
          
          # if not quiet, then print progress lambda
        } else if (trace == "print") {
          cat("\nFinished lambda = ", paste(loop.lambda, sep = ""))
        }
      }
      
      # if not quiet, then print progress kfold
      if (trace == "print") {
        cat("\nFinished fold ", paste(t, sep = ""))
      }
    }
    else {
      # designate folds and shuffle -- ensures randomized folds
      if (!is.null(seed)) {
        set.seed(seed)
      }
      ind = sample(n)
      
      # parse data into folds and perform CV
      for (t in 1:K) {
        # training set
        leave.out = ind[(1 + floor((t - 1) * n / K)):floor(t * n / K)]
        X.train = x[-leave.out, , drop = FALSE]
        covariates.train = covariates[-leave.out, , drop = FALSE]
        # # sample covariances
        X.train = scale(X.train, center = TRUE, scale = TRUE)
        empcov.train = cor(X.train)
        
        # loop over all tuning parameters
        for (loop.lambda in 1:length(lambda)) {
          # compute the penalized likelihood precision matrix
          # estimator
          fit = GGM_Estimation( x=X.train,
          known_ppi=known_ppi,
          covariates=covariates.train,
          scr=scr,
          gamma=gamma,
          lambda=lambda[loop.lambda],
          rho=rho,
          weight=weight,
          eta=eta,
          verbose=verbose,
          eps = eps
        )
        
          siginv = fit$invcov
          
          no.edge = (sum(abs(siginv) > eps)) / 2
          #CV_density[loop.lambda, t] = sparsePrec(fit$Adj)
          # compute the observed negative validation loglikelihood
          CV_nLogLikelihood[loop.lambda, t] = tr(empcov.train %*% siginv) - 
                                    determinant(siginv, logarithm = T)$modulus
          CV_pNLogLikelihood[loop.lambda, t] = CV_nLogLikelihood[loop.lambda, t] + 
                              fit$lambda *sum(abs((1 - diag(p)) * siginv))
          CV_bic.score[loop.lambda, t] = tr(empcov.train %*% siginv) - 
            determinant(siginv, logarithm = T)$modulus + log(n) * no.edge / n
          CV_Pbic.score[loop.lambda, t] = CV_pNLogLikelihood[loop.lambda, t] + 
                                          log(n) * no.edge / n
          CV_ebic.score[loop.lambda, t] = CV_bic.score[loop.lambda, t] + 4 *
                                          log(p) * no.edge / n
          CV_aic.score[loop.lambda, t] = CV_nLogLikelihood[loop.lambda, t] + 
                                          2 * no.edge / n
          # update progress bar
          if (trace == "progress") {
            setTxtProgressBar(progress, loop.lambda + (t - 1) * length(lambda))
            # if not quiet, then print progress lambda
          } else if (trace == "print") {
            cat("\nFinished lambda = ", paste(loop.lambda, sep = ""))
          }
        }
        # if not quiet, then print progress kfold
        if (trace == "print") {
          cat("\nFinished fold ", paste(t, sep = ""))
        }
      }
    }
    
    if (crit.cv == "loglik") {
      CV_errors = CV_nLogLikelihood
    }
    if (crit.cv == "ploglik") {
      CV_errors = CV_pNLogLikelihood
    }
    if (crit.cv == "PBIC") {
      CV_errors = CV_Pbic.score
    }
    if (crit.cv == "BIC") {
      CV_errors = CV_bic.score
    }
    if (crit.cv == "EBIC") {
      CV_errors = CV_ebic.score
    }
    if (crit.cv == "AIC") {
      CV_errors = CV_aic.score
    }
    
    allMetrics = cbind(
      lambda,
      #apply(CV_density, 1, mean),
      apply(CV_nLogLikelihood, 1, mean),
      apply(CV_pNLogLikelihood, 1, mean),
      apply(CV_bic.score, 1, mean),
      apply(CV_Pbic.score, 1, mean),
      apply(CV_ebic.score, 1, mean),
      apply(CV_aic.score, 1, mean)
    )
    colnames(allMetrics) <-
      c(
        "lambda",
        #"density",
        "nLogLikelihood",
        "pNLogLikelihood",
        "BIC",
        "PBIC",
        "EBIC",
        "AIC"
      )
    
    # determine optimal tuning parameters
    AVG = apply(CV_errors, 1, mean)
    STD = apply(CV_errors, 1, sd)
    if (K == 1) {
      bestLambda = lambda[which.min(AVG)]
      minerror = min(AVG)
    } else {
      minerror = min(AVG)
      minerror_sd = STD[which.min(AVG)]
      bestLambda = lambda[which.min(abs(AVG - (minerror + minerror_sd)))]
    }
    
    results = list(
      lambda = lambda,
      bestLambda = bestLambda,
      min.error = minerror,
      #CV_density = CV_density,
      avg.error = AVG,
      cv.error = CV_errors,
      CV_nLogLikelihood = CV_nLogLikelihood,
      CV_pNLogLikelihood = CV_pNLogLikelihood,
      CV_bic.score = CV_bic.score,
      CV_Pbic.score = CV_Pbic.score,
      CV_ebic.score = CV_ebic.score,
      CV_aic.score = CV_aic.score,
      allMetrics = allMetrics
    )
    class(results) = "CVahglasso"
    return(results)
    
  }

func.path <- "/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Previous_litt/Boxin_Zhao_FGM_Neighboorhood/Functions"
source(paste(func.path,"ADMM.new.R", sep="/"))

FGGReg_cov_estimation <- function(scores, # functional score on a defined basis, nrow: subjects; ncol: functions*n_basis.
                                  n_basis = 1, #Number of bases considered 
                                  covariates  = NULL, #Additional covariates to regress on
                                  scr = TRUE, # optional screening to speed up
                                  gamma = NULL, # Person correlation threshold for the screening
                                  lambda_prec = NULL, # Penalization term in the Lasso
                                  lambda_prec_type = "1se", # or "min"
                                  asparse = 0.75, # The relative weight to put on the `1-norm in sparse group lasso
                                  verbose = FALSE,
                                  eps = 1e-08) {
  
  ################################################
  ## PART 0. INITILA PREPROCESSING
  ################################################
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
  
  ################################################
  ## PART 1. INITIAL SCREENING
  ################################################
  
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
  
  ################################################
  ## PART 2. DEFINITION OF THE DESIGN MATRIX
  ################################################
  
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
  ## PART 3. COMPUTE THE MULTIPLE REGRESSION
  ################################################
  
  Delta_scores <- matrix(0, Mp, ncol(interM))
  colnames(Delta_scores) <- colnames(interM)
  rownames(Delta_scores) <- colnames(scores)
  No_sim_Delta_hat <- matrix(0, Mp, ncol(interM))
  colnames(No_sim_Delta_hat) <- colnames(interM)
  rownames(No_sim_Delta_hat) <- colnames(scores)
  

  ##################################################
  # DEVELOP THE SPARSE GROUP LASSO FOR MULTIPLE OUTCOME THEN MOVE ON
  ##################################################
  
  # Register parallel backend for the second loop
  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  if (verbose) {
        cat("Comuputation of adjacency matrix \n ")
      }

  results <- foreach(i = 1:p, .combine = 'cbind', .packages = c('glmnet', 'sparsegl')) %dopar% {
    Y <- matrix(Z0[, i], ncol = 1)
    protname <- colnames(Z0)[i]
    pattern <- paste0(":", protname, "\\b")
    indexes_to_select <- !grepl(pattern, colnames(interM))
    indexes_to_select[colnames(interM) == protname] <- FALSE
    indexes_to_select_penalty <- rep(TRUE, ncol(Z0))
    indexes_to_select_penalty[colnames(Z0) == protname] <- FALSE
    if (scr) {
      zero_col_indices <- which(scr_index[i, ] == 0)
      zero_col_names <- colnames(scr_index)[zero_col_indices]
      for (j in zero_col_names) {
        pattern <- paste0(":", j, "\\b")
        indexes_to_select[grepl(pattern, colnames(interM))] <- FALSE
        indexes_to_select[colnames(interM) == j] <- FALSE
        indexes_to_select_penalty[colnames(Z0) == j] <- FALSE
      }
    }
    indexes_to_select_reg <- indexes_to_select
    indexes_to_select_reg[1] <- FALSE
    Xmat <- interM[, indexes_to_select_reg]
    Xgroup <- groups[indexes_to_select_reg]
    
    w_sparsity <- rep(weight, ncol(Xmat))
    w_sparsity[Xgroup == 1] <- penaltyfactor[i, indexes_to_select_penalty]
    w_group[which(w_group==1)] <- sqrt(sum(Xgroup == 1))

    
    if (is.null(lambda_prec)) {
      set.seed(2024)
      mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup, pf_group = w_group,
                        pf_sparse = w_sparsity, asparse = asparse, intercept = FALSE)
      #mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup,
      #                  pf_sparse = w_sparsity, asparse = asparse, intercept = FALSE)
    } else {
      set.seed(2024)
      mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup, pf_group = w_group,
                        pf_sparse = w_sparsity, asparse = asparse, lambda = lambda_prec, intercept = FALSE)
      #mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup,
      #                  pf_sparse = w_sparsity, asparse = asparse, lambda = lambda_prec, intercept = FALSE)
    }
    if (lambda_prec_type == "min"){
    delta <- coef(mod, s = "lambda.min")[, 1]
    }else{
      delta <- coef(mod, s = "lambda.1se")[, 1]
    }
    
    res.model <- Y - cbind(1, as.matrix(Xmat)) %*% delta
    df_j <- sum(delta != 0)
    sigma <- sum(res.model^2) / (n - df_j)
    Delta_row <- rep(0, ncol(interM))
    Delta_row[indexes_to_select] <- as.vector(delta)
    No_sim_Delta_hat_row <- rep(0, ncol(interM))
    No_sim_Delta_hat_row[indexes_to_select] <- -sigma * as.vector(delta)
    list(Delta_row = Delta_row, No_sim_Delta_hat_row = No_sim_Delta_hat_row, sigma = sigma)
  }

  stopCluster(cl)

  # DEFINE THE ADIACENCY MATRIX OF THE FUNCTIONS INSTEAD OF THE SCORES
  Delta_fun <- matrix(0, p, ncol(interM))
  colnames(Delta_fun) <- colnames(interM)
  rownames(Delta_fun) <- colnames(scores)
  
  for (i in 1:p) {
    Delta[i, ] <- results[,i]$Delta_row
    No_sim_Delta_hat[i, ] <- results[,i]$No_sim_Delta_hat_row
    Sigma_hat[i] <- results[,i]$sigma
  }

  Dic_Delta_hat <- list()
  fact_groups <- as.factor(groups)
  if (!is.null(covariates)) {
    U_names <- colnames(C)[-1]
  } else {
    U_names <- c()
  }
  U_names <- c("Intercept", "Prot", U_names)

  for (k in 2:length(levels(fact_groups))) {
    id <- levels(fact_groups)[k]
    key <- U_names[k]
    value <- No_sim_Delta_hat[, fact_groups == id]
    Temp <- abs(value)
    Temp_transpose <- abs(t(value))
    max_matrix <- pmax(Temp, Temp_transpose)
    sim_value <- sign(value + t(value)) * max_matrix
    Dic_Delta_hat[[key]] <- sim_value
  }
  
  return(
    list(
      deltas = Delta,
      Sigma_hat = Sigma_hat,
      Dic_Delta_hat = Dic_Delta_hat
    )
  )
}


ADMM.grplasso <- function(A.X, A.Y, d, lambda,
                          rho.init=1, maxiter=2000, mu=10, tau=2,
                          P.in, Q.in, U.in, tol.rel=1e-4, tol.abs=1e-4){
  
  # Input:
  #   A.X, n*pM matrix. For convenience we denote dimension by p instead of p-1
  #   A.Y, n*M matrix
  #   d, a pM long vector
  
  # 1.1. Reading in the parameters
  n <- nrow(A.Y)
  M <- ncol(A.Y)
  p <- ncol(A.X) / M # corresponds to p-1 in theory
  
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
  for(k in 1:p){
    V.Grp[[k]] <- V[(k-1)*M + (1:M), ] # each Vk dim M * M
  }
  
  # 1.5.2. indicator of rho's change
  rho <- rho.init
  
  ##############################
  ## 2. Main part: Iterations ##
  ##############################
  for(iter in 1:maxiter){
    
    if(iter==1){
      rho.changed <- T # For iteration 1, we need it True for Q update
    }
    
    #############################
    ## 7.1. Updating Variables ##
    #############################
    # 7.1.1. Update P with group soft thresholding
    P.new <- matrix(NA, nrow=p*M, ncol=M)
    P.Grp <- grp.soft.thres(V.Grp, d, lambda, rho) # Group of p, each M * M
    for(k in 1:p){
      P.new[(k-1)*M + (1:M), ] <- P.Grp[[k]]
    }
    # 7.1.2. Update Q
    if(rho.changed){
      # if rho not changed, this term stays the same
      left.inv <- solve(DXXD.n + rho * diag(p*M))
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
    tol.prim <- tol.abs * sqrt(p*M*M) + tol.rel * max(norm(P.new, "F"), 
                                                      norm(Q.new, "F"))
    tol.dual <- tol.abs * sqrt(p*M*M) + tol.rel * norm(U.new, "F")
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
    V.Grp <- list()
    for(k in 1:p){
      V.Grp[[k]] <- V[(k-1)*M + (1:M), ] # each Vk dim M * n
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


############################# UPDATE grp.soft.thres
grp.soft.thres.multiple_groups <- function(V.Grp, M, d, lambda, rho){
  # V.Grp: a group of n_blocks matrices (Vk), each dim M * M
  # d: vector, length pM.
  # lambda: penalty param of group lasso
  # rho: penalty param of augmented Lagrangian
  
  n_blocks <- length(V.Grp)
  D.Grp <- list()
  P.Grp <- list() # output group of n_blocks matrices (Pk), each dim M * M
  for(k in 1:n_blocks){
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





#############################
### FUNCTIONS TO PERFROM TWO GROUPS COMPARISON
#############################

############################# UPDATE grp.soft.thres
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

############################# UPDATE objective.function
objective.function.two.groups <- function(A.X, A.Y, d, B, lambda){
  # A.X: matrix, n * 2(p-1)M
  # A.Y: matrix, n * M
  # d: vector, length 2(p-1)M
  # B: matrix, 2(p-1)M * M
  # Q: matrix, 2(p-1)M * M
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

############################# UPDATE ADMM.grplasso -> miss a parantesis or one more

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
  ##############################
  ## 2. Main part: Iterations ##
  ##############################
  for(iter in 1:maxiter){
    
    if(iter==1){
      rho.changed <- T # For iteration 1, we need it True for Q update
    }
    
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
      # if rho not changed, this term stays the same
      left.inv <- solve(DXXD.n + rho * diag(n_blocks*M))
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
    if(prim.resid < tol.prim & dual.resid < tol.dual){
      break
    }
    
    #####################################
    ## 7.3. Prepare for next iteration ##
    #####################################
    
    # 7.3.1. Update V
    V <- Q.new - U.new # pM * M
    V.Grp <- list()
    for(k in 1:n_blocks){
      V.Grp[[k]] <- V[(k-1)*M + (1:M), ] # each Vk dim M * n
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
    obj <- objective.function.two.groups(A.X=A.X, A.Y=A.Y, d=d, B=P.new, lambda=lambda)
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

############################# SAME AS THEIRS lambda.sup
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

############################# UPDATE lambda.sup.global
lambda.sup.global.two.groups <- function(design.matrix, M, p){
  
  n <- nrow(design.matrix)
  
  Y <- matrix(NA, n, M)
  X <- matrix(NA, n, 2*(p-1)*M)
  l.max <- rep(0,p)
  for(j in 1:p){
    jth.range.y <- (j-1)*M+(1:M)
    jth.range.x <- c(jth.range.y,(j+p-1)*M+(1:M))
    for(i in 1:n){
      Y[i,] <- as.matrix(design.matrix[i, jth.range.y])
      X[i,] <- as.matrix(design.matrix[i, -jth.range.x])
    }
    l.max[j] <- lambda.sup(X, Y)
  }
  lambda.max <- max(l.max)
  return(lambda.max)
}



FGGReg_diff_two_groups <- function(scores, # functional score on a defined basis, nrow: subjects; ncol: functions*n_basis.
                                  n_basis = 1, #Number of bases considered 
                                  covariates  = NULL, #Additional covariates to regress on
                                  L = 100, # How many penalization term in the Lasso to try
                                  thres.ctrl = 0, # recognition threshold epsilon_n = thres.ctrl * lambda_n,
                                  verbose = FALSE,
                                  tol.abs =1e-4 ,
                                  tol.rel = 1e-4,
                                  eps = 1e-08) {
  
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
  if (n < (q+1)*Mp) {
    warning("The sample size is too small! Network estimate may be unreliable!")
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
  
  ## COMPUTE THE MULTIPLE REGRESSION
  
  #The outcomes I'd like
  Delta_scores <- matrix(0, Mp, ncol(interM))
  colnames(Delta_scores) <- colnames(interM)
  rownames(Delta_scores) <- colnames(scores)
  No_sim_Delta_hat <- matrix(0, Mp, ncol(interM))
  colnames(No_sim_Delta_hat) <- colnames(interM)
  rownames(No_sim_Delta_hat) <- colnames(scores)
  
  #The outcomes they have
  time.start.gX <- proc.time()
  
  G.list.gX <- list()
  TPR.and.gX <- rep(0,L)
  FPR.and.gX <- rep(0,L)
  TPR.or.gX <- rep(0,L)
  FPR.or.gX <- rep(0,L)
  
  # Q: A cosa servono d.array e  d ? Quali sono le dimensioni corrette?
  # Per il momento li ho sostuiti com dei vettori di uni della dimensione coretta (?)
  d.array <- matrix(1, nrow=p, ncol=(p-1)*M*(q+1))
  d_out <- list()
  norm.adj <- rep(NA,p)
  for(k in 1:p){
    d_out[[k]] <- d.array[k,]
    norm.adj[k] <- norm(d_out[[k]],"2")
  }
  
  # Define the sequneces of lambda to try
  #TO DO: Add condition of doing so only if lambda is not provided
  lambda.max.gX <- lambda.sup.global.two.groups(interM,M, p)
  lambdas.gX <- seq(lambda.max.gX, 0, length.out=L)
  
  # default warm start PQU
  P.def <- matrix(0, 2*(p-1)*M, M)
  Q.def <- matrix(0.1, 2*(p-1)*M, M)
  U.def <- matrix(0.01, 2*(p-1)*M, M)
  
  # Repeat the following loop for each node
  V.array <- foreach(j = 1:p, .combine="rbind", .packages="fda") %dopar% {
    jth.range.y <- (j-1)*M+(1:M)
    A.Y <- as.matrix(interM[, jth.range.y])
    jth.range.x <- c(jth.range.y,(j+p-1)*M+(1:M))
    A.X <- as.matrix(interM[, -jth.range.x])
    groups <- temp_groups[-jth.range.x]
    
    P <- P.def; Q <- Q.def; U <- U.def
    
    
    V.j <- matrix(NA, L, 2*p)
    
    for(l in 1:L){
      lambda <- lambdas.gX[l]
      # Up to here
      grp.lasso.result <- ADMM.grplasso.two.groups(A.X = A.X, A.Y = A.Y, d = d_out[[j]],
                                        lambda=lambda, rho.init=1,
                                        P.in=P, Q.in=Q, U.in=U,
                                        tol.abs = tol.abs, tol.rel = tol.rel,
                                        maxiter = 400)
      P <- grp.lasso.result$P
      Q <- grp.lasso.result$Q
      U <- grp.lasso.result$U
      
      n_blocks <- nrow(P)/M
      # Process the estimated P into a neighborhood selection vector
      P.frob <- list()
      for(k in 1:n_blocks){
        key <- unique(groups[(k-1)*M + (1:M)])
        if(length(key) > 1){
            print("Error in group definition")
        } else {
          P.frob[[key]] <- norm(P[(k-1)*M + (1:M), ], "F")
        }
      }
      
      threshold <- lambda * thres.ctrl
      
      # Neighbor recognization
      V.jl <- rep(0, length(P.frob))
      for(n_block in 1:length(P.frob)){
          if(!(is.null(P.frob[[n_block]]))){
            if(P.frob[[n_block]] > threshold){
            V.jl[n_block] <- 1}
            }
      }
      
      V.j[l,] <- V.jl
    }
    V.j
  }
  
  time.end.gX <- proc.time()
  runtime.gX <- (time.end.gX - time.start.gX)[3]
  print("Computational time of: ")
  print(runtime.gX)
  return(V.array)
  }
  
 