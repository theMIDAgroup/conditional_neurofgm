rm(list = ls())

library("igraph")

# Step 1. Load data
folder = ""
path_sensors = paste0(folder, "EEG_dataset\\", "Position_list.Rdata")
path_adj = paste0(folder, "EEG_dataset\\seed_1\\results\\", "Test1_Adj_estimation.rda")

load(path_sensors)  # load pos.list
load(path_adj)  # load AdjMat

# Step 2. Define sensor names and positions
p <- 64
node.names <- numeric(p)
for (i in 1:p){
  node.names[i] <- pos.list[[i]]
}
position.List <- list("FPZ"=c(0, 0.8), "AFZ"=c(0, 0.6), "FZ"=c(0, 0.4), FCZ=c(0, 0.2), 
                      "CZ"=c(0, 0), "CPZ"=c(0, -0.2), "PZ"=c(0, -0.4), "POZ"=c(0, -0.6), 
                      "OZ"=c(0, -0.8), "nd"=c(0,-1), "C2"=c(0.2, 0), "C4"=c(0.4, 0), 
                      "C6"=c(0.6, 0), "T8"=c(0.8, 0), "Y"=c(1, 0), "C1"=c(-0.2, 0), 
                      "C3"=c(-0.4, 0), "C5"=c(-0.6, 0), "T7"=c(-0.8, 0), "X"=c(-1, 0), 
                      "FP2"=c(0.2, 0.76), "FP1"=c(-0.2, 0.76), "AF2"=c(0.26, 0.62), 
                      "AF1"=c(-0.26, 0.62), "AF8"=c(0.45, 0.68), "AF7"=c(-0.45, 0.68), 
                      "F2"=c(0.21, 0.41), "F1"=c(-0.21, 0.41), "F4"=c(0.4, 0.45), 
                      "F3"=c(-0.4, 0.45), "F6"=c(0.55, 0.5), "F5"=c(-0.55, 0.5), 
                      "F8"=c(0.65, 0.55), "F7"=c(-0.65, 0.55), "FC2"=c(0.25, 0.21), 
                      "FC1"=c(-0.25, 0.21), "FC4"=c(0.5, 0.22), "FC3"=c(-0.5, 0.22), 
                      "FC6"=c(0.7, 0.26), "FC5"=c(-0.7, 0.26), "FT8"=c(0.9, 0.31), 
                      "FT7"=c(-0.9, 0.31), "CP2"=c(0.25, -0.21), "CP4"=c(0.5, -0.22), 
                      "CP6"=c(0.7, -0.26), "TP8"=c(0.9, -0.31), "CP1"=c(-0.25, -0.21), 
                      "CP3"=c(-0.5, -0.22), "CP5"=c(-0.7, -0.26), "TP7"=c(-0.9, -0.31), 
                      "P2"=c(0.21, -0.41), "P4"=c(0.4, -0.45), "P6"=c(0.55, -0.5), 
                      "P8"=c(0.65, -0.55), "P1"=c(-0.21, -0.41), "P3"=c(-0.4, -0.45), 
                      "P5"=c(-0.55, -0.5), "P7"=c(-0.65, -0.55), "PO2"=c(0.2, -0.62), 
                      "PO8"=c(0.45, -0.68), "PO1"=c(-0.2, -0.62), "PO7"=c(-0.45, -0.68), 
                      "O2"=c(0.2, -0.85), "O1"=c(-0.2, -0.85))
layMat <- matrix(NA, nrow=64, ncol=2)
for (i in 1:length(node.names)){
  x <- node.names[i]
  layMat[i, 1] <- unlist(position.List)[paste(x, 1, sep="")]
  layMat[i, 2] <- unlist(position.List)[paste(x, 2, sep="")]
}

# Step 3. Prepare adjancecy matrix
adj_group = G.our.symm.weighted$group
adj_group = log10(adj_group)        # Logscale
adj_group[is.na(adj_group)] = 0
diag(adj_group) <- 0

vals <- adj_group[adj_group != 0 & !is.na(adj_group)]
greates_vals <- vals[abs(vals)>2]
perc_greatest <- length(greates_vals)/length(vals)

abs_thr <- quantile(abs(vals), probs =0.0001, na.rm = TRUE)
adj_extreme <- adj_group
adj_extreme[abs(adj_group) < abs_thr] <- 0

# Step 4. Plot
#pdf("diff_network_AUD.pdf", width = 11, height = 7)
adj_plot = adj_extreme
colnames(adj_plot) <- node.names
row.names(adj_plot) <- node.names
net <- graph_from_adjacency_matrix(adj_plot, mode = "undirected",
            weighted = TRUE,diag = FALSE)
w <- E(net)$weight
max_range = max(abs(w))
ncol <- 1000
col_fun <- colorRampPalette(c("darkblue", "white", "darkred"))
cols <- col_fun(ncol)
E(net)$color <- cols[
  as.numeric(cut(w, breaks = seq(-max_range, max_range, length.out = ncol+1),
                 include.lowest = TRUE))
]
layout(matrix(c(1, 2, 3), nrow = 1), widths = c(4, 0.5, 0.5))
par(mar = c(0, 0, 0, 0))
ord <- order(abs(E(net)$weight))
# Rebuild graph with reordered edges
net2 <- graph_from_data_frame(
  d = as_data_frame(net, what = "edges")[ord, ],
  directed = FALSE,
  vertices = as_data_frame(net, what = "vertices")
)
#net <- net %>% igraph::permute.edges(ord)
#alpha_vals <- scales::rescale(abs(E(net2)$weight), to = c(0.3, 1))
E(net2)$color <- adjustcolor(E(net2)$color, alpha.f = 0.7)
plot(net2, layout = layMat,
     vertex.size = 3, edge.width=2, vertex.label.color = "black", vertex.label.font = 2,
     vertex.label.dist=1, vertex.label.cex = 2, vertex.color = "orange", 
     edge.curved = 0.2, margin =0)
# ---- ColorbarS ----
par(fig = c(0.8, 0.83, 0.2, 0.8), new = TRUE, mar = c(2,0,2,0))
image(
  z = t(matrix(seq(-max_range, max_range, length.out = ncol+1), ncol = 1)),
  col = cols,
  xaxt='n',
  yaxt='n'
)
axis(4, at=seq(0,1,length.out=5), cex.axis = 2,
     labels=round(seq(-max_range, max_range, length.out=5),2))

mtext(
  expression('Edge weight (log'[10]*' scale)'),
  side = 4,
  line = 4,
  cex = 1.7
)

box()
#dev.off()



