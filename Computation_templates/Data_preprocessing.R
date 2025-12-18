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
print("STARTING: Proprocessing of the data")
#############################################
###### PART 0: Upload all the data and libraries needed
#############################################
output_path <- config$output_path
name_output <- config$name_output

observed_functional_data_path <- config$observed_functional_data_path

p <- config$n_nodes
tau <- config$time_points 
time_range <- config$time_range


rec_basis_type <- config$rec_basis_type
rec_basis_number <- config$rec_basis_number
M <- config$dim_red_basis_number

suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(matrixcalc))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(mvtnorm))
suppressPackageStartupMessages(library(abind))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(pheatmap))
dir.create(output_path)

h <- get(load(observed_functional_data_path))
####################################
#     PART 2: GAIN FPC SCORE       #
####################################
if (rec_basis_type == "bsplines"){
  basis <- create.bspline.basis(rangeval=c(time_range[1], time_range[2]), nbasis=rec_basis_number)
} else{
  basis <- create.fourier.basis(rangeval=c(time_range[1], time_range[2]),nbasis=rec_basis_number)
}

obs.time <- seq(time_range[1] + (time_range[2]-time_range[1])/tau, time_range[2], (time_range[2]-time_range[1])/tau)

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

write.csv(fpc.score, paste(output_path, "fpc_scores_", name_output, ".csv", sep=""))
scores_df <- fpc.score
save(scores_df, file=paste(output_path,"scores_df_", name_output, ".RData", sep=""))
print(paste0("END: Preprocessing of the data completed in: ", Sys.time()-time.start))

