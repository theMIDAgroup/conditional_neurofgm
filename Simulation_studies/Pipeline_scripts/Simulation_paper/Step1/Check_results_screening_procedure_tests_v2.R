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
  #G.pop <- (G.true.g2>0) & (G.true.g1>0)
  G.pop <- (G.true.g1>0)
  diag(G.pop) <- FALSE
  G.pop <-ifelse(G.pop, 1, 0)
  G.our.pop <- G.our[,1:num_nodes]
  G.our.pop <- ifelse(G.our.pop, 1, 0)

  adj_df <- melt(G.our.pop)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Node adjacency Matrix - Baseline group", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Adjacency_matrix_node_pop.png", sep=""), width=8, height=8, dpi=300)

  
  res_and <- prec.rec(G.pop, G.our.pop, type="AND")
  results_metices <- rbind(results_metices,c("POP","AND",res_and$prec,
                                             res_and$TPR, res_and$FPR, res_and$F1,
                                             n_g1, n_g2))
  res_or <- prec.rec(G.pop, G.our.pop, type="OR")
  results_metices <- rbind(results_metices,c("POP","OR",res_or$prec,
                                             res_or$TPR, res_or$FPR, res_or$F1,
                                             n_g1, n_g2))
  G.diff <- (abs(G.true.g2 - G.true.g1))>0
  diag(G.diff) <- FALSE
  G.diff <-ifelse(G.diff, 1, 0)
  
  G.our.diff <- G.our[,(num_nodes+1):(num_nodes*n_groups)]

  adj_df <- melt(G.our.diff)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Node adjacency Matrix - Differential connection", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Adjacency_matrix_node_diff.png", sep=""), width=8, height=8, dpi=300)

  
  res_and <- prec.rec(G.diff, G.our.diff, type="AND")
  results_metices <- rbind(results_metices,c("DIFF","AND",res_and$prec,
                                             res_and$TPR, res_and$FPR, res_and$F1,
                                             n_g1, n_g2))
  res_or <- prec.rec(G.diff, G.our.diff, type="OR")
  results_metices <- rbind(results_metices,c("DIFF","OR",res_or$prec,
                                             res_or$TPR, res_or$FPR, res_or$F1,
                                             n_g1, n_g2))

  
  G.our.full.estimation <- matrix(0, nrow = num_nodes, ncol = num_nodes)

  for (i in 1:num_nodes) {
    for (j in 1:num_nodes) {
      if (G.our.pop[i, j] == 1 && G.our.diff[i, j] == 1) {
        G.our.full.estimation[i, j] <- -2  # reduced (present in both)
      } else if (G.our.pop[i, j] == 0 && G.our.diff[i, j] == 1) {
        G.our.full.estimation[i, j] <- 2   # additional (only in diff)
      } else if (G.our.pop[i, j] == 1 && G.our.diff[i, j] == 0) {
        G.our.full.estimation[i, j] <- 1   # common (only in pop)
      }
    }
  }

  # Save the matrix
  save(G.our.full.estimation, file = file.path(output_path, "G_our_full_estimation.rda"))

  # Prepare data for plotting
  adj_df <- melt(G.our.full.estimation)
  colnames(adj_df) <- c("Row", "Col", "Value")

 # Create plot
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = as.factor(Value))) +
    geom_tile() +
    scale_fill_manual(
      values = c("-2" = "blue", "1" = "black", "2" = "red", "0"="white"),
      labels = c("-2" = "Reduced edges", "1" = "Common edges", "2" = "Enlarge edges", "0"="Edges not present"),
      name = "Edge Type"
    ) +
    theme_minimal() +
    coord_fixed() +
    scale_y_reverse() +
    labs(title = "Full Estimated Adjacency Matrix - Differential contribution highlighted", x = "", y = "")

  # Save the plot
  ggsave(plot, filename = file.path(output_path, "Adjacency_matrix_node_full.png"), width = 8, height = 8, dpi = 300)

  adj_df <- adj_df %>%
    mutate(Value = ifelse(Value %in% c(1,2), 1, 0)) 

  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black") +
    theme_minimal() +
    coord_fixed() +
    scale_y_reverse() +  # to match matrix view
    labs(title = "Node adjacency Matrix - Differential group", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Adjacency_matrix_node_diff_group.png", sep=""), width=8, height=8, dpi=300)

  
  
  G.our.ratio.forb <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
  G.our.sum.forb <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
  G.our.sum.forb.logic <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
  G.our.ratio.forb.sum <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
  G.our.group <- G.our.pop
  for(j in 1: num_nodes){
    for(k in 1: num_nodes){
      if(G.our.diff[j,k]){
        load(paste(output_path, name_output,"_" ,j, "coeff.rda", sep=""))
        beta.g2.forbs <- norm(P.values[[k]] + P.values[[k+num_nodes]], "F")
        G.our.group[j,k] <- ifelse(beta.g2.forbs > threshold,1,0)
        if(P.frob[[k]] == 0){
          cat("Warning: Zero norm at denominator. Setting to the numerator value")
          G.our.ratio.forb[j,k] <- (P.frob[[k+num_nodes]])
          G.our.ratio.forb.sum[j,k] <- beta.g2.forbs
        }else{
          G.our.ratio.forb[j,k] <- (P.frob[[k+num_nodes]])/(P.frob[[k]])
          G.our.ratio.forb.sum[j,k] <- beta.g2.forbs/(P.frob[[k]])
        }
        G.our.sum.forb[j,k] <- beta.g2.forbs
        G.our.sum.forb.logic[j,k] <- beta.g2.forbs>threshold
      }
    }
  }
  
  G.group <- (G.true.g2>0)
  diag(G.group) <- FALSE
  G.group <-ifelse(G.group, 1, 0)
  
  adj_df <- melt(G.our.group)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black") +
    theme_minimal() +
    coord_fixed() +
    scale_y_reverse() +  # to match matrix view
    labs(title = "Node adjacency Matrix - Gorup connection", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Adjacency_matrix_node_diff_group_value_based.png", sep=""), width=8, height=8, dpi=300)
  
  
  res_and <- prec.rec(G.group, G.our.group, type="AND")
  results_metices <- rbind(results_metices,c("GROUP","AND",res_and$prec,
                                             res_and$TPR, res_and$FPR, res_and$F1,
                                             n_g1, n_g2))
  res_or <- prec.rec(G.group, G.our.group, type="OR")
  results_metices <- rbind(results_metices,c("GROUP","OR",res_or$prec,
                                             res_or$TPR, res_or$FPR, res_or$F1,
                                             n_g1, n_g2))
  
  adj_df <- melt(G.our.ratio.forb.sum, varnames = c("Row","Col"), value.name = "Value")
  
  # Option A — map directly with midpoint at 1
  p <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 1,               # fade to white near 1
      na.value = "white",
      name = "Value"
    ) +
    coord_fixed() +
    scale_y_reverse() +
    theme_minimal() +
    labs(title = "Changes in the estimated Adjacency Matrix",
         x = NULL, y = NULL)
  
  ggsave(p, filename=paste(output_path, "Adjacency_matrix_node_diff_weights.png", sep=""), width=8, height=8, dpi=300)
  
  adj_df <- adj_df %>%
    mutate(Value_logic = case_when(
      is.na(Value)      ~ "0",
      Value < 1         ~ "-2",
      Value > 1         ~ "2",
      TRUE              ~ "NA"
    ))
  
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = as.factor(Value_logic))) +
    geom_tile() +
    scale_fill_manual(
      values = c("-2" = "blue", "1" = "black", "2" = "red", "0"="white"),
      labels = c("-2" = "Reduced edges", "1" = "Common edges", "2" = "Enlarge edges", "0"="Edges not present"),
      name = "Edge Type"
    ) +
    theme_minimal() +
    coord_fixed() +
    scale_y_reverse() +
    labs(title = "Full Estimated Adjacency Matrix - Differential contribution highlighted", x = "", y = "")
  
  # Save the plot
  ggsave(plot, filename = file.path(output_path, "Adjacency_matrix_node_diff_weights_logic.png"), width = 8, height = 8, dpi = 300)
  
}
results_metices <- results_metices[-1,]
write.csv(results_metices, paste(save_path, simulation_name, "/test_results_metices_",name_output, ".csv", sep=""))
cat("Results saved to ", paste(save_path, simulation_name, "/test_results_metices_",name_output, ".csv", sep=""), "\n")
#################################################
