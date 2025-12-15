rm(list=ls(all=TRUE))
packages <- c('yaml')
install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

#################################################
## USER DEFINED PARAMETERS (MODIFY THE PATH TO THE CORRECT YAML FILE)
#################################################
args <- commandArgs(trailingOnly = TRUE)
yaml_file_path = args[[1]]
config <- yaml.load_file(yaml_file_path)
time.start <- Sys.time()

#################################################
## 0. USER DEFINED PARAMETERS (MODIFY THIS PART)
#################################################
save_path = config$save_path
simulation_name = config$simulation_name
iteration = config$iteration
name_output = config$name_output
p <- config$p
n_basis <- config$M
output_path = paste0(save_path, simulation_name,"/", "seed_", iteration, "/", "results_lit_comparison/")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

#################################################
## 0.1 LOAD FUNCTIONS
#################################################
folder.dir <- "/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies/Pipeline_scripts/Simulation_paper/Step1/FuDGE"
func.path <- paste0(folder.dir, "/Function_Definitions")
model.path <- paste0(folder.dir, "/Example_Model2")


# Load the library and functions
source(paste(func.path,"ProxAlg_DFGM.R", sep="/"))
source(paste(func.path,"FFGL_ADMM.R", sep="/"))
source(paste(func.path,"FFGL2_ADMM.R", sep="/"))


suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(matrixcalc))
suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(quadprog))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(JGL))

score_path = paste0(save_path, simulation_name,"/", "seed_", iteration, "/", "fpc_scores_", name_output,".csv")
grouping_path = paste0(save_path, simulation_name,"/", "seed_", iteration, "/", "grouping_factor_", name_output,".csv")

scores <- read.csv(score_path)[, -1]
n_nodes <- ncol(scores)/n_basis
n_samples <- nrow(scores)
covariates <- data.frame( group = as.factor(read.csv(grouping_path)[,-1]))
full_data <- cbind(covariates,scores)
#str(full_data)


################### FuDGE #########################

#### Data Preparation
principle.score.X <- full_data[full_data$group == "1",][, -1]
principle.score.Y <- full_data[full_data$group == "2",][, -1]

principle.score.X.cen <- scale(principle.score.X, center=T, scale=F)
principle.score.Y.cen <- scale(principle.score.Y, center=T, scale=F)

estimated.cov.pc.X <- (t(principle.score.X.cen) %*% principle.score.X.cen) / (n_samples-1)
estimated.cov.pc.Y <- (t(principle.score.Y.cen) %*% principle.score.Y.cen) / (n_samples-1)

##### Setting hyperparams

lambdamax <- 2
lambdamin <- 0
steplength <- 0.05
lambda.v.fudge <- seq(lambdamax, lambdamin, -steplength)
Delta.initial <- matrix(0, nrow=(p*n_basis), ncol=(p*n_basis))

################### Estimate the Differential Graph #########################
AdjMats_lambda_FuDGE <- list()
for (lambda.ind in c(1:length(lambda.v.fudge)) ){
  
  ################### FuDGE #########################
  lambda.choose <- lambda.v.fudge[lambda.ind]
  DFGM_g <- ProxAlg_DFGM(estimated.cov.pc.X, estimated.cov.pc.Y, p, n_basis, 
                         lambda.choose, Intialization=Delta.initial)
  
  # Update initialization (warm start)
  Delta.initial <- DFGM_g$DelatMathat
  
  # Compute TPR and FPR
  FrobMat <- DFGM_g$blockFrob
  
  Edgeshat <- matrix(0, nrow=p, ncol=p)
  Edgeshat[which(FrobMat > 0)] <- 1
  AdjMats_lambda_FuDGE[[lambda.ind]]<- list(Diff=Edgeshat)
  
}

print(paste0("Computational time of: ",Sys.time() - time.start,"\n"))


full_result_path = paste(output_path, name_output, "_FUDGE_comparison", ".rda", sep ="")
save(AdjMats_lambda_FuDGE, file =full_result_path)
cat("\n Output saved to: ", full_result_path,"\n" )
