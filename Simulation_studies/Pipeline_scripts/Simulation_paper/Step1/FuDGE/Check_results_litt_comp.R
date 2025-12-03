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
tot_iteration = config$tot_iteration
name_output = config$name_output
num_nodes <- config$p
n_comp <- config$M
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

  if(sum(G.mat) ==0 & sum(G.true)==0 ){
    TPR <- 0
    FPR <- 0
    F1 <- 1
  }
  return(list(prec=prec, rec=rec, TPR=TPR, FPR=FPR, F1=F1))
}
#################################################

#################################################
##### Compute the metrices
results_metices <- data.frame(network=NA,	method=NA, hyper=NA,	prec=NA,	TPR=NA,	FPR=NA,	F1=NA,
                              Baseline=NA,	Differential=NA)
for(iteration in 1:tot_iteration){
  cat("Processing iteration ", iteration, "\n")
  output_path = paste0(save_path, simulation_name,"/", "seed_", iteration, "/", "results_lit_comparison/")
  load(paste(output_path, name_output,"_litt_comparison", ".rda", sep=""))
  foldname = paste0(save_path, simulation_name,"/", "seed_", iteration)
  load(paste(foldname,"/Ground_truth_",name_output,".rda", sep=""))
  #G.pop <- (G.true.g2>0) & (G.true.g1>0)
  G.pop <- (G.true.g1>0)
  diag(G.pop) <- FALSE
  G.pop <-ifelse(G.pop, 1, 0)
  G.diff <- (abs(G.true.g2 - G.true.g1))>0
  diag(G.diff) <- FALSE
  G.diff <-ifelse(G.diff, 1, 0)
  G.group <- (G.true.g2>0)
  diag(G.group) <- FALSE
  G.group <-ifelse(G.group, 1, 0)
  
  G.FuDGE.diff.opt <- NULL
  G.FuDGE.diff.worse <- NULL
  F1.FuDGE.opt <- 0
  F1.FuDGE.worse <- 1
  
  G.FGL.diff.opt <- NULL
  G.FGL.diff.worse <- NULL
  F1.FGL.diff.opt <- 0
  F1.FGL.diff.worse <- 1
  
  G.FGL.pop.opt <- NULL
  G.FGL.pop.worse <- NULL
  F1.FGL.pop.opt <- 0
  F1.FGL.pop.worse <- 1

  G.FGL.group.opt <- NULL
  G.FGL.group.worse <- NULL
  F1.FGL.group.opt <- 0
  F1.FGL.group.worse <- 1


  G.FFGL.diff.opt <- NULL
  G.FFGL.diff.worse <- NULL
  F1.FFGL.diff.opt <- 0
  F1.FFGL.diff.worse <- 1
  
  G.FFGL.pop.opt <- NULL
  G.FFGL.pop.worse <- NULL
  F1.FFGL.pop.opt <- 0
  F1.FFGL.pop.worse <- 1

  G.FFGL.group.opt <- NULL
  G.FFGL.group.worse <- NULL
  F1.FFGL.group.opt <- 0
  F1.FFGL.group.worse <- 1

  G.FFGL2.diff.opt <- NULL
  G.FFGL2.diff.worse <- NULL
  F1.FFGL2.diff.opt <- 0
  F1.FFGL2.diff.worse <- 1
  
  G.FFGL2.pop.opt <- NULL
  G.FFGL2.pop.worse <- NULL
  F1.FFGL2.pop.opt <- 0
  F1.FFGL2.pop.worse <- 1

  G.FFGL2.group.opt <- NULL
  G.FFGL2.group.worse <- NULL
  F1.FFGL2.group.opt <- 0
  F1.FFGL2.group.worse <- 1


  for(lam in 1: length(AdjMats_lambda_FuDGE)){

    ################### FuDGE #########################
    G.FuDGE <- AdjMats_lambda_FuDGE[[lam]]
    G.FuDGE.diff <- G.FuDGE$Diff
    res_FuDGE <- prec.rec(G.diff, G.FuDGE.diff, type="AND")
    results_metices <- rbind(results_metices,c("DIFF","FuDGE", lam, res_FuDGE$prec,
                                               res_FuDGE$TPR, res_FuDGE$FPR, res_FuDGE$F1,
                                               n_g1, n_g2))
    if(res_FuDGE$F1 > F1.FuDGE.opt){
      F1.FuDGE.opt = res_FuDGE$F1
      G.FuDGE.diff.opt <- G.FuDGE.diff
    }
    
    if(res_FuDGE$F1 < F1.FuDGE.worse){
      F1.FuDGE.worse = res_FuDGE$F1
      G.FuDGE.diff.worse <- G.FuDGE.diff
    }
    
    ################### FGL #########################
    G.FGL <- AdjMats_lambda_FGL[[lam]]
    
    ## DIFERENTIAL
    G.FGL.diff <- G.FGL$Diff
    res_FGL <- prec.rec(G.diff, G.FGL.diff, type="AND")
    results_metices <- rbind(results_metices,c("DIFF","FGL", lam, res_FGL$prec,
                                               res_FGL$TPR, res_FGL$FPR, res_FGL$F1,
                                               n_g1, n_g2))
    
    if(res_FGL$F1 > F1.FGL.diff.opt){
      F1.FGL.diff.opt = res_FGL$F1
      G.FGL.diff.opt <- G.FGL.diff
    }
    
    if(res_FGL$F1 < F1.FGL.diff.worse){
      F1.FGL.diff.worse = res_FGL$F1
      G.FGL.diff.worse <- G.FGL.diff
    }
    
    ## POPULATION
    G.FGL.pop <- G.FGL$Pop
    res_FGL <- prec.rec(G.pop, G.FGL.pop, type="AND")
    results_metices <- rbind(results_metices,c("POP","FGL", lam, res_FGL$prec,
                                               res_FGL$TPR, res_FGL$FPR, res_FGL$F1,
                                               n_g1, n_g2))
    
    if(res_FGL$F1 > F1.FGL.pop.opt){
      F1.FGL.pop.opt = res_FGL$F1
      G.FGL.pop.opt <- G.FGL.pop
    }
    
    if(res_FGL$F1 < F1.FGL.pop.worse){
      F1.FGL.pop.worse = res_FGL$F1
      G.FGL.pop.worse <- G.FGL.pop
    }
    
    ## GROUP
    G.FGL.group <- G.FGL$Group
    res_FGL <- prec.rec(G.group, G.FGL.group, type="AND")
    results_metices <- rbind(results_metices,c("Group","FGL", lam, res_FGL$prec,
                                               res_FGL$TPR, res_FGL$FPR, res_FGL$F1,
                                               n_g1, n_g2))
    
    if(res_FGL$F1 > F1.FGL.group.opt){
      F1.FGL.group.opt = res_FGL$F1
      G.FGL.group.opt <- G.FGL.group
    }
    
    if(res_FGL$F1 < F1.FGL.group.worse){
      F1.FGL.group.worse = res_FGL$F1
      G.FGL.group.worse <- G.FGL.group
    }

    ################### FFGL #########################
    G.FFGL <- AdjMats_lambda_FFGL[[lam]]
    
    ## DIFERENTIAL
    G.FFGL.diff <- G.FFGL$Diff
    res_FFGL <- prec.rec(G.diff, G.FFGL.diff, type="AND")
    results_metices <- rbind(results_metices,c("DIFF","FFGL", lam, res_FFGL$prec,
                                               res_FFGL$TPR, res_FFGL$FPR, res_FFGL$F1,
                                               n_g1, n_g2))
    
    if(res_FFGL$F1 > F1.FFGL.diff.opt){
      F1.FFGL.diff.opt = res_FFGL$F1
      G.FFGL.diff.opt <- G.FFGL.diff
    }
    
    if(res_FFGL$F1 < F1.FFGL.diff.worse){
      F1.FFGL.diff.worse = res_FFGL$F1
      G.FFGL.diff.worse <- G.FFGL.diff
    }
    
    ## POPULATION
    G.FFGL.pop <- G.FFGL$Pop
    res_FFGL <- prec.rec(G.pop, G.FFGL.pop, type="AND")
    results_metices <- rbind(results_metices,c("POP","FFGL", lam, res_FFGL$prec,
                                               res_FFGL$TPR, res_FFGL$FPR, res_FFGL$F1,
                                               n_g1, n_g2))
    
    if(res_FFGL$F1 > F1.FFGL.pop.opt){
      F1.FFGL.pop.opt = res_FFGL$F1
      G.FFGL.pop.opt <- G.FFGL.pop
    }
    
    if(res_FFGL$F1 < F1.FFGL.pop.worse){
      F1.FFGL.pop.worse = res_FFGL$F1
      G.FFGL.pop.worse <- G.FFGL.pop
    }
    
    ## GROUP
    G.FFGL.group <- G.FFGL$Group
    res_FFGL <- prec.rec(G.group, G.FFGL.group, type="AND")
    results_metices <- rbind(results_metices,c("Group","FFGL", lam, res_FFGL$prec,
                                               res_FFGL$TPR, res_FFGL$FPR, res_FFGL$F1,
                                               n_g1, n_g2))
    
    if(res_FFGL$F1 > F1.FFGL.group.opt){
      F1.FFGL.group.opt = res_FFGL$F1
      G.FFGL.group.opt <- G.FFGL.group
    }
    
    if(res_FFGL$F1 < F1.FFGL.group.worse){
      F1.FFGL.group.worse = res_FFGL$F1
      G.FFGL.group.worse <- G.FFGL.group
    }
    
    ################### FFGL2 #########################
    G.FFGL2 <- AdjMats_lambda_FFGL2[[lam]]
    
    ## DIFERENTIAL
    G.FFGL2.diff <- G.FFGL2$Diff
    res_FFGL2 <- prec.rec(G.diff, G.FFGL2.diff, type="AND")
    results_metices <- rbind(results_metices,c("DIFF","FFGL2", lam, res_FFGL2$prec,
                                               res_FFGL2$TPR, res_FFGL2$FPR, res_FFGL2$F1,
                                               n_g1, n_g2))
    
    if(res_FFGL2$F1 > F1.FFGL2.diff.opt){
      F1.FFGL2.diff.opt = res_FFGL2$F1
      G.FFGL2.diff.opt <- G.FFGL2.diff
    }
    
    if(res_FFGL2$F1 < F1.FFGL2.diff.worse){
      F1.FFGL2.diff.worse = res_FFGL2$F1
      G.FFGL2.diff.worse <- G.FFGL2.diff
    }
    
    ## POPULATION
    G.FFGL2.pop <- G.FFGL2$Pop
    res_FFGL2 <- prec.rec(G.pop, G.FFGL2.pop, type="AND")
    results_metices <- rbind(results_metices,c("POP","FFGL2", lam, res_FFGL2$prec,
                                               res_FFGL2$TPR, res_FFGL2$FPR, res_FFGL2$F1,
                                               n_g1, n_g2))
    
    if(res_FFGL2$F1 > F1.FFGL2.pop.opt){
      F1.FFGL2.pop.opt = res_FFGL2$F1
      G.FFGL2.pop.opt <- G.FFGL2.pop
    }
    
    if(res_FFGL2$F1 < F1.FFGL2.pop.worse){
      F1.FFGL2.pop.worse = res_FFGL2$F1
      G.FFGL2.pop.worse <- G.FFGL2.pop
    }
    
    ## GROUP
    G.FFGL2.group <- G.FFGL2$Group
    res_FFGL2 <- prec.rec(G.group, G.FFGL2.group, type="AND")
    results_metices <- rbind(results_metices,c("Group","FFGL2", lam, res_FFGL2$prec,
                                               res_FFGL2$TPR, res_FFGL2$FPR, res_FFGL2$F1,
                                               n_g1, n_g2))
    
    if(res_FFGL2$F1 > F1.FFGL2.group.opt){
      F1.FFGL2.group.opt = res_FFGL2$F1
      G.FFGL2.group.opt <- G.FFGL2.group
    }
    
    if(res_FFGL2$F1 < F1.FFGL2.group.worse){
      F1.FFGL2.group.worse = res_FFGL2$F1
      G.FFGL2.group.worse <- G.FFGL2.group
    }

    
  }
  
  adj_df <- melt(G.FuDGE.diff.opt)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Diff Node adjacency Matrix - FuDGE (best) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Diff_Adjacency_matrix_FuDGE_best.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FuDGE.diff.worse)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Diff Node adjacency Matrix - FuDGE (worse) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Diff_Adjacency_matrix_FuDGE_worse.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FGL.diff.opt)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Diff Node adjacency Matrix - FGL (best) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Diff_Adjacency_matrix_FGL_best.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FGL.diff.worse)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Diff Node adjacency Matrix - FGL (worse) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Diff_Adjacency_matrix_FGL_worse.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FGL.pop.opt)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Pop Node adjacency Matrix - FGL (best) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Pop_Adjacency_matrix_FGL_best.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FGL.pop.worse)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Pop Node adjacency Matrix - FGL (worse) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Pop_Adjacency_matrix_FGL_worse.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FGL.group.opt)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Group Node adjacency Matrix - FGL (best) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Group_Adjacency_matrix_FGL_best.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FGL.group.worse)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Group Node adjacency Matrix - FGL (worse) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Group_Adjacency_matrix_FGL_worse.png", sep=""), width=8, height=8, dpi=300)


  adj_df <- melt(G.FFGL.diff.opt)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Diff Node adjacency Matrix - FFGL (best) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Diff_Adjacency_matrix_FFGL_best.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FFGL.diff.worse)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Diff Node adjacency Matrix - FFGL (worse) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Diff_Adjacency_matrix_FFGL_worse.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FFGL.pop.opt)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Pop Node adjacency Matrix - FFGL (best) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Pop_Adjacency_matrix_FFGL_best.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FFGL.pop.worse)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Pop Node adjacency Matrix - FFGL (worse) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Pop_Adjacency_matrix_FFGL_worse.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FFGL.group.opt)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Group Node adjacency Matrix - FFGL (best) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Group_Adjacency_matrix_FFGL_best.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FFGL.group.worse)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Group Node adjacency Matrix - FFGL (worse) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Group_Adjacency_matrix_FFGL_worse.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FFGL2.diff.opt)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Diff Node adjacency Matrix - FFGL2 (best) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Diff_Adjacency_matrix_FFGL2_best.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FFGL2.diff.worse)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Diff Node adjacency Matrix - FFGL2 (worse) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Diff_Adjacency_matrix_FFGL2_worse.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FFGL2.pop.opt)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Pop Node adjacency Matrix - FFGL2 (best) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Pop_Adjacency_matrix_FFGL2_best.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FFGL2.pop.worse)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Pop Node adjacency Matrix - FFGL2 (worse) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Pop_Adjacency_matrix_FFGL2_worse.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FFGL2.group.opt)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Group Node adjacency Matrix - FFGL2 (best) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Group_Adjacency_matrix_FFGL2_best.png", sep=""), width=8, height=8, dpi=300)

  adj_df <- melt(G.FFGL2.group.worse)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Group Node adjacency Matrix - FFGL2 (worse) ", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Group_Adjacency_matrix_FFGL2_worse.png", sep=""), width=8, height=8, dpi=300)

 
}
results_metices <- results_metices[-1,]
write.csv(results_metices, paste(save_path, simulation_name, "/litt_comp_test_results_metices_",name_output, ".csv", sep=""))
cat("Results saved to ", paste(save_path, simulation_name, "/litt_comp_test_results_metices_",name_output, ".csv", sep=""), "\n")
#################################################
