# Required packages
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


time.start <- Sys.time()
print("STARTING: Simulation of the data")
#############################################
###### PART 0: Upload all the data and libraries needed
#############################################
func.path <- config$func_path
save.path <- config$save_path
simualtion.name <- config$simulation_name
iteration <- config$iteration
seed_g1 <- config$seed_g1
seed_g2 <- config$seed_g2
name_output <- config$name_output

p <- config$p
n_g1 <- config$n_g1
n_g2 <- config$n_g2
mu <- config$mu
tau <- config$tau 
model_g1 <- config$model_g1
model_g2 <- config$model_g2
red.number <- config$red_number 

rec_basis_type <- config$rec_basis_type
rec_basis_number <- config$rec_basis_number
M <- config$M

suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(matrixcalc))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(mvtnorm))
suppressPackageStartupMessages(library(abind))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(pheatmap))
source(paste(func.path,"FPCA.score.R", sep="/"))    # Computing FPCA scores
source(paste(func.path,"bases.func.R", sep="/"))    # Basis functions for FPCA
dir.create(save.path)
dir.create(paste0(save.path, simualtion.name))

####################################
#     PART 1: DATA GENERATION      
####################################
foldname = paste0(save.path, simualtion.name,"/", "seed_", iteration, "/")
dir.create(foldname)

# Generate Random Functions and Observation h_ijk       
# h_ijk = Fourier basis func(t_k) %*% delta_ij + iid error

# Generate precision matrix and real adjacency matrix for G1
source(paste(func.path,model_g1, sep="/")) # For Model A generation
set.seed(seed_g1*iteration)
if(model_g1 == "A.prelim.func.R"){
  Theta.g1 <- cov.mat.model.A(p, mu) # p*mu by p*mu large square matrix
}else if (model_g1== "A.prelim.func.red.R"){
  Theta.g1 <- cov.mat.model.A.red(p, mu, red.number) # p*mu by p*mu large square matrix
}


adj_df <- melt(Theta.g1)
colnames(adj_df) <- c("Row", "Col", "Value")

# Plot
plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Score adjacency Matrix", x = "", y = "")
ggsave(plot, filename=paste(foldname, "/Adjacency_matrix_scores_g1_", name_output, ".png", sep=""), width=8, height=8, dpi=300)


G.true.g1 <- matrix(0, p, p) # p by p adjacency matrix
for(i in 1:p){
  for(j in 1:p){
    if(sum(abs(Theta.g1[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])) > 0)
      G.true.g1[i,j] <- frobenius.norm(Theta.g1[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])
  }
}
adj_df <- melt(G.true.g1)
colnames(adj_df) <- c("Row", "Col", "Value")
plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Node adjacency Matrix", x = "", y = "")
ggsave(plot, filename=paste(foldname, "/Adjacency_matrix_node_g1_",name_output , ".png", sep=""), width=8, height=8, dpi=300)

# Generating delta
set.seed(seed_g1*iteration)
delta <- rmvnorm(n_g1, sigma = solve(Theta.g1))

# Observation time
obs.time <- seq(1/tau, 1, 1/tau) # vector of observation time points of delta

# Fourier basis function for data generation
b.mat.list <- list()
for(j in 1:p){
  b.mat.list[[j]] <- fda.fourier.mat(obs.time, mu)
}

# 4. Observations h_ijk
h <- array(0, c(n_g1, p, tau))
for(i in 1:n_g1){
  for(j in 1:p){
    set.seed(seed_g1*iteration + j)
    h[i,j,] <- b.mat.list[[j]] %*% 
      matrix(delta[i, ((j-1)*mu+1) : (j*mu)], ncol=1) + rnorm(tau, 0, 0.5)
  }
}

h_g1 <- h
png(paste(foldname, "/Data_simulated_g1_",name_output , ".png", sep=""), width = 1200, height = 800, res = 150)
matplot(t(h_g1[,1,]), type = "l", lty = 1, col = 1:6)
dev.off()

# Generate precision matrix and real adjacency matrix for G2
source(paste(func.path,model_g2, sep="/")) # For Model A generation
set.seed(seed_g2*iteration)
if(model_g2 == "A.prelim.func.R"){
  Theta.g2 <- cov.mat.model.A(p, mu) # p*mu by p*mu large square matrix
} else if (model_g2== "A.prelim.func.red.R"){
  Theta.g2 <- cov.mat.model.A.red(p, mu, red.number ) # p*mu by p*mu large square matrix
}

adj_df <- melt(Theta.g2)
colnames(adj_df) <- c("Row", "Col", "Value")

# Plot
plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Score adjacency Matrix", x = "", y = "")
ggsave(plot, filename=paste(foldname, "/Adjacency_matrix_scores_g2_", name_output, ".png", sep=""), width=8, height=8, dpi=300)

G.true.g2 <- matrix(0, p, p) # p by p adjacency matrix
for(i in 1:p){
  for(j in 1:p){
    if(sum(abs(Theta.g2[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])) > 0)
      G.true.g2[i,j] <- frobenius.norm(Theta.g2[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])
  }
}
adj_df <- melt(G.true.g2)
colnames(adj_df) <- c("Row", "Col", "Value")
plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  coord_fixed() +
  scale_y_reverse() +  # to match matrix view
  labs(title = "Node adjacency Matrix", x = "", y = "")
ggsave(plot, filename=paste(foldname, "/Adjacency_matrix_node_g2_", name_output, ".png", sep=""), width=8, height=8, dpi=300)

# 1. Generating delta
set.seed(seed_g2*iteration)
delta <- rmvnorm(n_g2, sigma = solve(Theta.g2))
# 2. Observation time
obs.time <- seq(1/tau, 1, 1/tau) # vector of observation time points of delta
# 3. Fourier basis function for data generation
b.mat.list <- list()
for(j in 1:p){
  b.mat.list[[j]] <- fda.fourier.mat(obs.time, mu)
}

# 4. Observations h_ijk
h <- array(0, c(n_g2, p, tau))
for(i in 1:n_g2){
  for(j in 1:p){
    set.seed(seed_g2*iteration + j)
    h[i,j,] <- b.mat.list[[j]] %*% 
      matrix(delta[i, ((j-1)*mu+1) : (j*mu)], ncol=1) + rnorm(tau, 0, 0.5)
  }
}

h_g2 <- h
png(paste(foldname, "/Data_simulated_g2_",name_output , ".png", sep=""), width = 1200, height = 800, res = 150)
matplot(t(h_g2[,1,]), type = "l", lty = 1, col = 1:6)
dev.off()
h <- abind(h_g1, h_g2, along = 1)
group=c(rep(1,n_g1), rep(2,n_g2))
save(h,group, file=paste(foldname, "/Original_data_generated_", name_output, ".RData", sep=""))

save(G.true.g1, G.true.g2, file=paste(foldname,"/Ground_truth_", name_output, ".rda", sep=""))


####################################
#     PART 2: GAIN FPC SCORE       #
####################################

if (rec_basis_type == "bsplines"){
  basis <- create.bspline.basis(rangeval=c(0, 1), nbasis=rec_basis_number)
} else{
  basis <- create.fourier.basis(rangeval=c(0,1),nbasis=rec_basis_number)
}

fpc.score <- numeric(0)
for(j in 1:p){
  obs.val.matrix <- matrix(0, nrow=tau, ncol=dim(h)[1])
  for (i in c(1:dim(h)[1])){
    obs.val.vec <- as.vector(h[i, j, ])
    obs.val.matrix[, i] <- obs.val.vec
  }
  fd.object.array <- Data2fd(argvals=obs.time, y=obs.val.matrix, basisobj=basis)
  fpc.score <- cbind(fpc.score, pca.fd(fd.object.array, nharm=M)$scores)
}

n_nodes <- ncol(fpc.score)/M
n_samples <- nrow(fpc.score)
names <- rep(NA, ncol(fpc.score))
for(l in 1:ncol(fpc.score)){
  names[l] <- paste("f",ceiling(l/M) ,".",l%%M, sep ="")
}
colnames(fpc.score) <- names

write.csv(fpc.score, paste(foldname, "fpc_scores_", name_output, ".csv", sep=""))
write.csv(data.frame(group=group),paste(foldname, "grouping_factor_", name_output, ".csv", sep=""))

fpc.score <- as.data.frame(fpc.score)
scoredata_matrix <- as.matrix(fpc.score)
Scor <- as.data.frame(cor(scoredata_matrix))
colnames(Scor) <- colnames(fpc.score)
rownames(Scor) <- colnames(fpc.score)

save(Scor, file=paste(foldname,"Score_corr_matrix_", name_output, ".rda", sep=""))

node_Scor <- matrix(NA, nrow = p, ncol = p)
for (i in 1:p) {
  for (j in 1:p) {
    row_idx <- ((i - 1) * M + 1):(i * M)
    col_idx <- ((j - 1) * M + 1):(j * M)
    
    block <- abs(Scor[row_idx, col_idx])
    node_Scor[i, j] <- max(block, na.rm = TRUE)
  }
}

names_node_corr <- rep(NA, p)
for(l in 1:p){
  names_node_corr[l] <- paste("f",l, sep ="")
}
colnames(node_Scor) <- names_node_corr
rownames(node_Scor) <- names_node_corr

save(node_Scor, file=paste(foldname,"Node_corr_matrix_", name_output, ".rda", sep=""))
print("END: Simulation of the data completed in: ")
Sys.time()-time.start
output_path = paste0(foldname, "results/")
dir.create(output_path)
