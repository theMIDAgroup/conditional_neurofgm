setwd("/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies")
rm(list=ls(all=TRUE))

# List of required libraries
required_libraries <- c("parallel", "doParallel", "foreach", "rBayesianOptimization")

# Check and install missing packages
for (lib in required_libraries) {
  if (!require(lib, character.only = TRUE)) {
    message(paste("Installing missing package:", lib))
    install.packages(lib, dependencies = TRUE)
    library(lib, character.only = TRUE)
  } else {
    message(paste("Package", lib, "is already installed."))
  }
}

data.path <- "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1"
func.path <- "/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Graph_estimation"
save.path <- "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1"
runtime.path <- "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1"

source(paste(func.path,"cFGGM_functions_two_groups.R", sep="/"))

scores <- read.csv(paste(data.path, "freq_scores_reord.csv", sep ="/"))[, -1]
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


cl <- makeCluster(16)
registerDoParallel(cl)

L = 10
K = 5
thres.ctrl = c(0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0)

G.mat <- FGGReg_diff_two_groups_SCV_bayes_search(scores, # functional score on a defined basis, nrow: subjects; ncol: functions*n_basis.
                                         n_basis = n_basis, #Number of bases considered 
                                         covariates  = covariates, #Additional covariates to regress on
                                         L = L, # How many penalization term in the Lasso to try
                                         K = K,
                                         thres.ctrl = thres.ctrl, # recognition threshold epsilon_n = thres.ctrl * lambda_n,
                                         verbose = FALSE,
                                         tol.abs =1e-4 ,
                                         tol.rel = 1e-4,
                                         eps = 1e-08,
                                         name_log= paste(save.path,"/logs/try0104_VDM", sep=""))

save(G.mat, file = paste(save.path,"res_0104_Bayes_param.RData", spe="/"))
stopCluster(cl)



print( Sys.time() - start )
