load("/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1/res_0304_sequential.RData")
setwd("/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1")
G.pop <- read.table("adj_pop.txt")
View(G.pop)
G.group <- read.table("adj_group.txt")

  
num_nodes <- 64

# Create an empty adjacency matrix (64x64) initialized to 0
adj_matrix_pop <- matrix(0, nrow = num_nodes, ncol = num_nodes)

# Fill the adjacency matrix based on the edges
for (i in 1:nrow(G.pop)) {
  node1 <- G.pop$V1[i]
  node2 <- G.pop$V2[i]
  
  # Set the matrix entries to 1 (assuming an unweighted graph)
  adj_matrix_pop[node1, node2] <- 1
  adj_matrix_pop[node2, node1] <- 1  # If the graph is undirected
}

# Print the adjacency matrix
print(adj_matrix_pop)

# Create an empty adjacency matrix (64x64) initialized to 0
adj_matrix_g <- matrix(0, nrow = num_nodes, ncol = num_nodes)

# Fill the adjacency matrix based on the edges
for (i in 1:nrow(G.group)) {
  node1 <- G.group$V1[i]
  node2 <- G.group$V2[i]
  
  # Set the matrix entries to 1 (assuming an unweighted graph)
  adj_matrix_g[node1, node2] <- 1
  adj_matrix_g[node2, node1] <- 1  # If the graph is undirected
}

# Print the adjacency matrix
print(adj_matrix_g)

G.mat <-  ifelse(G.mat, 1, 0)
G.mat 

G.true <- cbind(adj_matrix_pop, adj_matrix_g)

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
  if(TP+FP==0) prec <- 1
  
  rec <- TP / (TP + FN)
  if(TP+FN==0) rec <- 0
  
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  F1 <- 2*TP/(2*TP+FP+FN)
  return(list(prec=prec, rec=rec, TPR=TPR, FPR=FPR, F1=F1))
}


prec.rec(G.true, G.mat, type="AND")
prec.rec(G.true, G.mat, type="OR")

G.mat.pop <-  ifelse(G.mat[,1:64], 1, 0)
G.mat.pop 

prec.rec(adj_matrix_pop, G.mat.pop, type="AND")
prec.rec(adj_matrix_pop, G.mat.pop, type="OR")




G.mat.g <-  ifelse(G.mat[,65:128], 1, 0)
G.mat.g

prec.rec(adj_matrix_g, G.mat.g, type="AND")
prec.rec(adj_matrix_g, G.mat.g, type="OR")


#### USING OR
G.mat.pop.sim <- G.mat.pop + t(G.mat.pop)
G.mat.pop.sim <- ifelse(G.mat.pop.sim>0, 1, 0)
table(adj_matrix_pop,G.mat.pop.sim)
##### USING AND
G.mat.pop.sim <- G.mat.pop + t(G.mat.pop)
G.mat.pop.sim <- ifelse(G.mat.pop.sim>1, 1, 0)
table(adj_matrix_pop,G.mat.pop.sim)

sum(G.mat.g)
table(adj_matrix_g, G.mat.g)

#### USING OR
G.mat.pop.sim <- G.mat.pop + t(G.mat.pop)
G.mat.pop.sim <- ifelse(G.mat.pop.sim>0, 1, 0)
table(adj_matrix_g,G.mat.pop.sim)

##### USING AND
G.mat.pop.sim <- G.mat.pop + t(G.mat.pop)
G.mat.pop.sim <- ifelse(G.mat.pop.sim>1, 1, 0)
table(adj_matrix_g,G.mat.pop.sim)

