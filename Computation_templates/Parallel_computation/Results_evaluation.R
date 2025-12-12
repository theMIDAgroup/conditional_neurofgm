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
output_path = config$output_path
name_output = config$name_output
input_path = config$input_path
type = config$type

load(input_path)

if ((exists("covariates_df", envir = .GlobalEnv) && is.data.frame(covariates_df))){
  C <- model.matrix(~ ., covariates_df)
  n_groups <- ncol(C)
  cov_names <- colnames(C)
  cov_names[1] <- "population"
} else{ n_groups <- 1
cov_names <- "population" }

num_nodes = config$n_nodes
n_comp <- config$n_basis_for_dim_recustion

##############################################

G.our <- matrix(NA, nrow = num_nodes, ncol = num_nodes*n_groups)
for(j in 1: num_nodes){
    load(paste(output_path, name_output,"_" ,j, ".rda", sep=""))
    G.our[j, ] <- N.hat.optimal
}

G.our.symm = list()
G.our.symm.weighted = list()

for(i in 1:n_groups){
  G.cov <- G.our[,(num_nodes*(i-1) +1):(num_nodes*i)]
  G.cov.simm <- matrix(0, nrow = num_nodes, ncol = num_nodes)
  if(type=="AND"){
    for(j in 1:(num_nodes-1)){
      for(k in (i+1):num_nodes){
        if(G.cov[j,k] == 1 & G.cov[k,j] == 1){
          G.cov.simm[j,k] = 1
          G.cov.simm[k,j] = 1
        }
      }
    }
  }else if(type=="OR"){
    for(j in 1:(num_nodes-1)){
      for(k in (i+1):num_nodes){
        if(G.cov[j,k] == 1 | G.cov[k,j] == 1){
          G.cov.simm[j,k] = 1
          G.cov.simm[k,j] = 1
        }
      }
    }
  }else{
    G.cov.simm = G.cov
  }
  
  if(i==1){
    plot_title <- paste0("Node adjacency Matrix - ", cov_names[i])
  } else {
    plot_title <- paste0("Node adjacency Matrix - Differential contribution ", cov_names[i])
  }
  

  G.our.symm[[cov_names[i]]] <- G.cov.simm

  adj_df <- melt(G.cov.simm)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = plot_title, x = "", y = "")
  ggsave(plot, filename=paste0(output_path, name_output, "_Adjacency_matrix_node_",cov_names[i], ".png"), width=8, height=8, dpi=300)


  if(i !=1){
    G.cov.simm.weighted <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
    for(j in 1: num_nodes){
      for(k in 1: num_nodes){
        if(G.cov[j,k]){
          load(paste(output_path, name_output,"_" ,j, "coeff.rda", sep=""))
          beta.g2.forbs <- norm(P.values[[k]] + P.values[[k+num_nodes]], "F")
          if(P.frob[[k]] == 0){
            cat("Warning: Zero norm at denominator. Setting to the numerator value")
            G.cov.simm.weighted[j,k] <- beta.g2.forbs
          }else{
            G.cov.simm.weighted[j,k] <- beta.g2.forbs/(P.frob[[k]])
          }
        }
      }
    }


  
  G.our.symm.weighted[[cov_names[i]]] <- G.cov.simm.weighted

  adj_df <- melt(G.cov.simm.weighted, varnames = c("Row","Col"), value.name = "Value")
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
    labs(title = paste0("Full Estimated Adjacency Matrix - Differential contribution ", cov_names[i]), x = "", y = "")

  ggsave(plot, filename = paste0(output_path, name_output, "_Adjacency_matrix_node_diff_weights_",cov_names[i],".png"), width = 8, height = 8, dpi = 300)
}
  
}

# Save the matrix
save(G.our.symm,G.our.symm.weighted, file = paste0(output_path, name_output, "_Adj_estimation.rda"))
cat("Results saved to ", paste(output_path, name_output, "_Adj_estimation.rda", sep=""), "\n")
#################################################
