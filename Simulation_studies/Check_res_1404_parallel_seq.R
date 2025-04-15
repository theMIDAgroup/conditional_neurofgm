rm(list=ls(all=TRUE))
setwd("/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1")

#result_folder <- "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1/results/fourier_scores/"
#result_folder <- "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1/results/pca_scores/"
result_folder <- "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1/results/original_scores/"
result_name <- "rand_hyper_search_"

num_nodes <- 64
n_groups <- 2
G.mat <- matrix(NA, nrow = num_nodes, ncol = num_nodes*n_groups)
for(j in 1: num_nodes){
  load(paste(result_folder, result_name, j, ".rda", sep=""))
  G.mat[j, ] <- N.hat.optimal
}

G.pop <- read.table("adj_pop.txt")
View(G.pop)
G.group <- read.table("adj_group.txt")

adj_matrix_pop <- matrix(0, nrow = num_nodes, ncol = num_nodes)
for (i in 1:nrow(G.pop)) {
  node1 <- G.pop$V1[i]
  node2 <- G.pop$V2[i]
  
  adj_matrix_pop[(node1+1), (node2+1)] <- 1
  adj_matrix_pop[(node2+1), (node1+1)] <- 1
}

# Create an empty adjacency matrix (64x64) initialized to 0
adj_matrix_g <- matrix(0, nrow = num_nodes, ncol = num_nodes)

for (i in 1:nrow(G.group)) {
  node1 <- G.group$V1[i]
  node2 <- G.group$V2[i]
  
  # Set the matrix entries to 1 (assuming an unweighted graph)
  adj_matrix_g[(node1+1), (node2+1)] <- 1
  adj_matrix_g[(node2+1), (node1+1)] <- 1  # If the graph is undirected
}


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
sum(G.mat[,65:128])

prec.rec(adj_matrix_g, G.mat.g, type="AND")
prec.rec(adj_matrix_g, G.mat.g, type="OR")


#### USING OR
G.mat.pop.sim <- G.mat.pop + t(G.mat.pop)
G.mat.pop.sim <- ifelse(G.mat.pop.sim>0, 1, 0)
table(adj_matrix_pop,G.mat.pop.sim)
table(adj_matrix_g,G.mat.pop.sim)
##### USING AND
G.mat.pop.sim <- G.mat.pop + t(G.mat.pop)
G.mat.pop.sim <- ifelse(G.mat.pop.sim>1, 1, 0)
table(adj_matrix_pop,G.mat.pop.sim)

sum(G.mat.g)

#### USING OR
G.mat.g.sim <- G.mat.g+ t(G.mat.g)
G.mat.g.sim <- ifelse(G.mat.g.sim>0, 1, 0)
table(adj_matrix_g,G.mat.g.sim)

##### USING AND
G.mat.g.sim <- G.mat.g + t(G.mat.g)
G.mat.g.sim <- ifelse(G.mat.g.sim>1, 1, 0)
table(adj_matrix_g,G.mat.g.sim)


############################################
## CHECK PARAL and SEQUENTIAL give the same results
############################################

G.mat.parall <- G.mat
load("/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1/res_1404_sequential.RData")

dim(G.mat) #8192 entries
sum(G.mat == G.mat.parall) #Computation is the same

############################################
## CHECK ZHAO
############################################

rm(list=ls(all=TRUE))
setwd("/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1")

save.path <- "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1/results/zhao_1g"
load(paste(save.path,"/SCV.A.050.RunInd",1234,".Rdata",sep=""))

prec.rec(adj_matrix_pop, G.mat.gX, type="AND")
prec.rec(adj_matrix_pop, G.mat.gX, type="OR")

