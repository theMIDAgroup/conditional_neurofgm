rm(list=ls(all=TRUE))
packages <- c('yaml')
install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)
suppressPackageStartupMessages(library(yaml))

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
tot_iteration = config$tot_iteration
name_output = config$name_output
num_nodes <- config$p
n_groups <- 2
n_g1 <- config$n_g1
n_g2 <- config$n_g2

##############################################

#################################################
## Define the neede functions
prec.rec <- function(G.true, G.mat, type=c("AND","OR")){
  # a function to calculate TP, FP, TN, FN
  # and return precision, recall, TPR and FPR.
  # Input:
  #   G.true, the true p*p adjacency matrix
  #   G.mat, the estimated p*p adjacency matrix
  #   type,
  #     AND: when two nodes both recognize each other as neighbors
  #     OR: when either of them recognize each other as neighbors
  # Output:
  #   prec, precision
  #   rec, recall
  #   TPR and FPR, as their names
  
  
  p <- nrow(G.true)
  TP <- 0; TN <- 0; FP <- 0; FN <- 0
  
  if(type=="AND"){
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        if(G.mat[i,j] == 1 & G.mat[j,i] == 1){
          if(G.true[i,j] == 1)
            TP <- TP + 1
          else
            FP <- FP + 1
        }else{
          if(G.true[i,j] == 1)
            FN <- FN + 1
          else
            TN <- TN + 1
        }
      }
    }
  }else if(type=="OR"){
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        if(G.mat[i,j] == 1 | G.mat[j,i] == 1){
          if(G.true[i,j] == 1)
            TP <- TP + 1
          else
            FP <- FP + 1
        }else{
          if(G.true[i,j] == 1)
            FN <- FN + 1
          else
            TN <- TN + 1
        }
      }
    }
  }
  prec <- TP / (TP + FP)
  if(TP+FP==0) prec <- 0
  
  rec <- TP / (TP + FN)
  if(TP+FN==0) rec <- 0
  
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  F1 <- 2*TP/(2*TP+FP+FN)
  return(list(prec=prec, rec=rec, TPR=TPR, FPR=FPR, F1=F1))
}
#################################################

#################################################
##### Compute the metrices
results_metices <- data.frame(network=NA,	symm=NA,	prec=NA,	TPR=NA,	FPR=NA,	F1=NA,
                              Baseline=NA,	Differential=NA)
for(iteration in 1:tot_iteration){
  cat("Processing iteration ", iteration, "\n")
  output_path = paste0(save_path, simulation_name,"/", "seed_", iteration, "/", "results/")
  G.our <- matrix(NA, nrow = num_nodes, ncol = num_nodes*n_groups)
  for(j in 1: num_nodes){
    load(paste(output_path, name_output,"_" ,j, ".rda", sep=""))
    G.our[j, ] <- N.hat.optimal
  }
  foldname = paste0(save_path, simulation_name,"/", "seed_", iteration)
  load(paste(foldname,"/Ground_truth_",name_output,".rda", sep=""))
  G.pop <- (G.true.g2>0) & (G.true.g1>0)
  diag(G.pop) <- FALSE
  G.pop <-ifelse(G.pop, 1, 0)
  G.our.pop <- G.our[,1:num_nodes]
  G.our.pop <- ifelse(G.our.pop, 1, 0)
  
  res_and <- prec.rec(G.pop, G.our.pop, type="AND")
  results_metices <- rbind(results_metices,c("POP","AND",res_and$prec,
                                             res_and$TPR, res_and$FPR, res_and$F1,
                                             n_g1, n_g2))
  res_or <- prec.rec(G.pop, G.our.pop, type="OR")
  results_metices <- rbind(results_metices,c("POP","OR",res_or$prec,
                                             res_or$TPR, res_or$FPR, res_or$F1,
                                             n_g1, n_g2))
  G.diff <- (abs(G.true.g2 - G.true.g1)>0)
  diag(G.diff) <- FALSE
  G.diff <-ifelse(G.diff, 1, 0)
  
  G.our.diff <- G.our[,(num_nodes+1):(num_nodes*n_groups)]
  
  res_and <- prec.rec(G.diff, G.our.diff, type="AND")
  results_metices <- rbind(results_metices,c("DIFF","AND",res_and$prec,
                                             res_and$TPR, res_and$FPR, res_and$F1,
                                             n_g1, n_g2))
  res_or <- prec.rec(G.diff, G.our.diff, type="OR")
  results_metices <- rbind(results_metices,c("DIFF","OR",res_or$prec,
                                             res_or$TPR, res_or$FPR, res_or$F1,
                                             n_g1, n_g2))
  
}
results_metices <- results_metices[-1,]
write.csv(results_metices, paste(save_path, simulation_name, "/results_metices_",name_output, ".csv", sep=""))
cat("Results saved to ", paste(save_path, simulation_name, "/results_metices_",name_output, ".csv", sep=""), "\n")
#################################################
