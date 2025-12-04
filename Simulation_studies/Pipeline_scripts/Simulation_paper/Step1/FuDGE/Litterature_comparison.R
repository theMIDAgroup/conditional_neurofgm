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


################### FGL #########################

#### Data Preparation
prin.score <- list(X=principle.score.X, Y=principle.score.Y)

##### Setting hyperparams
lambda1.choose.fgl <- 0.01
lambda2.v.fgl <- rev(seq(0, 1, length.out = length(lambda.v.fudge)))

################### FFGL/FFGL2 #########################

#### Data Preparation
cov.list <- list(X=estimated.cov.pc.X, Y=estimated.cov.pc.Y)
n.sam <- c(n_samples, n_samples)

##### Setting hyperparams
lambda2.v.ffgl <- rev(seq(0, 100, length.out = length(lambda.v.fudge)))
lambda1.choose.ffgl <- 0.1

################### Estimate the Differential Graph #########################
time.start <- Sys.time()



AdjMats_lambda_FuDGE <- list()
AdjMats_lambda_FGL <- list()
AdjMats_lambda_FFGL <- list()
AdjMats_lambda_FFGL2 <- list()
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
  
  ################### FGL #########################
  
  lambda2.choose <- lambda2.v.fgl[lambda.ind]
  JGL.result <- JGL(prin.score, penalty="fused", lambda1=lambda1.choose.fgl, 
                    lambda2=lambda2.choose, penalize.diagonal=FALSE, 
                    return.whole.theta=TRUE)
  
  ThetaX <- JGL.result$theta[[1]]
  ThetaY <- JGL.result$theta[[2]]
  Edgeshat <- matrix(0, nrow=p, ncol=p)
  Edgeshat.pop <- matrix(0, nrow=p, ncol=p)
  Edgeshat.group <- matrix(0, nrow=p, ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      diff <- frobenius.norm(ThetaX[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)] - 
                               ThetaY[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      pop <- frobenius.norm(ThetaX[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      group <- frobenius.norm(ThetaY[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      if (diff>0){
        Edgeshat[i,j] <- 1
      }
      if (pop>0){
         Edgeshat.pop[i,j] <- 1
      }
      if (group>0){
         Edgeshat.group[i,j] <- 1
      }

      
    }
  }
  AdjMats_lambda_FGL[[lambda.ind]]<- list(Pop=Edgeshat.pop,Group=Edgeshat.group, Diff=Edgeshat)
  
  ################### FFGL #########################
  
  lambda2.choose <- lambda2.v.ffgl[lambda.ind]
  FFGL.result <- FFGL_ADMM(cov.list, n.sam, p, n_basis, lambda1.choose.ffgl, lambda2.choose)
  ThetaX <- FFGL.result$est[[1]]
  ThetaY <- FFGL.result$est[[2]]
  Edgeshat <- matrix(0, nrow=p, ncol=p)
  Edgeshat.pop <- matrix(0, nrow=p, ncol=p)
  Edgeshat.group <- matrix(0, nrow=p, ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      diff <- frobenius.norm(ThetaX[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)] - 
                               ThetaY[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      pop <- frobenius.norm(ThetaX[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      group <- frobenius.norm(ThetaY[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      if (diff>0.02){ #I would change to zero
        Edgeshat[i,j] <- 1
      }
      if (pop>0.02){ #I would change to zero
         Edgeshat.pop[i,j] <- 1
      }
      if (group>0.02){ #I would change to zero
         Edgeshat.group[i,j] <- 1
      }

      
    }
  }
  AdjMats_lambda_FFGL[[lambda.ind]]<- list(Pop=Edgeshat.pop,Group=Edgeshat.group, Diff=Edgeshat)

  ################### FFGL2 #########################
  
  FFGL.result <- FFGL2_ADMM(cov.list, n.sam, p, n_basis, lambda1.choose.ffgl, lambda2.choose)
  hetaX <- FFGL.result$est[[1]]
  ThetaY <- FFGL.result$est[[2]]
  Edgeshat <- matrix(0, nrow=p, ncol=p)
  Edgeshat.pop <- matrix(0, nrow=p, ncol=p)
  Edgeshat.group <- matrix(0, nrow=p, ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      diff <- frobenius.norm(ThetaX[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)] - 
                               ThetaY[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      pop <- frobenius.norm(ThetaX[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      group <- frobenius.norm(ThetaY[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      if (diff>0.02){ #I would change to zero
        Edgeshat[i,j] <- 1
      }
      if (pop>0.02){ #I would change to zero
         Edgeshat.pop[i,j] <- 1
      }
      if (group>0.02){ #I would change to zero
         Edgeshat.group[i,j] <- 1
      }

      
    }
  }
  AdjMats_lambda_FFGL2[[lambda.ind]]<- list(Pop=Edgeshat.pop,Group=Edgeshat.group, Diff=Edgeshat)
  
}

print(paste0("Computational time of: ",Sys.time() - time.start,"\n"))


full_result_path = paste(output_path, name_output, "_litt_comparison", ".rda", sep ="")
save(AdjMats_lambda_FuDGE, AdjMats_lambda_FGL,AdjMats_lambda_FFGL, AdjMats_lambda_FFGL2, file =full_result_path)
cat("\n Output saved to: ", full_result_path,"\n" )
