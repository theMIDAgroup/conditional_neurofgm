rm(list = ls())
# setwd("C:\\Users\\saras\\OneDrive\ -\ unige.it\\Documenti\\network_differenziali")

library(igraph)
library(ggraph)
library(extrafont)
library(patchwork)
library (pracma)
folder = "EEG_dataset/"
folder_res = "EEG_dataset/seed_1/results/"


plot_single_graph <- function(adjm, title){
  
  graph <- graph.adjacency(adjm, weighted = TRUE, mode = "undirected")
  
  node_names <- read.csv(
    paste0(folder, "position_list.txt"),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  V(graph)$name <- unlist(node_names, use.names = FALSE)
  
  gp <- ggraph(graph, layout = "linear", circular = TRUE) +
    
    geom_edge_arc(
      aes(edge_colour = ifelse(weight == 0, NA, weight)),
      show.legend = TRUE
    ) +
    
    geom_node_point(shape = 21, size = 3, fill = "darkgrey") +
    
    theme_graph(base_family = "Calibri") +
    
    scale_edge_colour_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      na.value = "transparent",
      oob = scales::squish_infinite,
      breaks = scales::pretty_breaks(n = 5),
      name = "Edge weight"
    ) +
    
    guides(
      fill = guide_colorbar(
        title.position = "top",
        barheight = unit(4, "cm"),
        barwidth  = unit(0.4, "cm")
      )
    ) +
    
    coord_fixed(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4)) +
    
    geom_node_text(
      aes(
        label = name,
        x = x * 1.2,
        y = y * 1.2,
        hjust = ifelse(
          atan2(y, x) > -pi/2 & atan2(y, x) < pi/2, 0, 1
        )
      ),
      size = 2
    ) +
    
    ggtitle(title)
  
  return(gp)
}



plot_diff_graph <- function(adjm0, adjm1,title,  leg=0){
  
  
  # --- Coding scheme
  coddiff <- data.frame(
    rel   = c("Equal", "Add", "Rem", "Incr", "Red"),
    value = c(0.1, 1, -1, 0.5, -0.5),
    color = c("black", "palegreen2", "red", "palegreen2", "red"),
    style = c("solid", "solid", "solid", "dashed", "dashed"),
    stringsAsFactors = FALSE
  )
  
  diff <- matrix(0, nrow = nrow(adjm0), ncol = ncol(adjm1))
  diff[adjm0 & adjm1 & adjm0 == adjm1] <- coddiff$value[coddiff$rel == "Equal"]
  diff[adjm0 & !adjm1]                 <- coddiff$value[coddiff$rel == "Rem"]
  diff[!adjm0 & adjm1]                 <- coddiff$value[coddiff$rel == "Add"]
  diff[adjm0 & adjm1 & adjm0 < adjm1] <- coddiff$value[coddiff$rel == "Incr"]
  diff[adjm0 & adjm1 & adjm0 > adjm1] <- coddiff$value[coddiff$rel == "Red"]
  
  
  graph <- graph.adjacency(diff, weighted = TRUE, mode = "undirected")
  
  node_names <- read.csv(paste0(folder, "position_list.txt"), header = FALSE, stringsAsFactors = FALSE) 
  node_names <- unlist(node_names, use.names = FALSE) # add names to the graph 
  V(graph)$name <- node_names
  
  values_col <- setNames(coddiff$color, coddiff$value)
  values_style <- setNames(coddiff$style, coddiff$value)
  labels_legend <- setNames(c("In both groups", "Only in G1", "Only in G0", "Reduced in G1", "Enhanced in G1"), coddiff$value)
  E(graph)$weight <- factor(E(graph)$weight, levels = coddiff$value)   
  # --- Prepare weight as factor with *all* possible levels
  
  # --- Plot
  gp <- ggraph(graph, layout = "linear", circular = TRUE) + 
    geom_edge_arc(aes(color = weight, linetype=weight)) + 
    geom_node_point(shape = 21, size = 3, fill = "darkgrey") +
    theme_graph(base_family = "Calibri") +
    scale_edge_color_manual(values = values_col, labels = labels_legend,
                            drop = FALSE) +  # <- keeps all legend entries
    scale_edge_linetype_manual(values = values_style, labels = labels_legend,
                               drop = FALSE) + 
    coord_fixed(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4)) +
    guides(fill = guide_none()) + ggtitle(title) +
    geom_node_text(
      aes(
        label = name,
        x = x * 1.2,
        y = y * 1.2,
        hjust = ifelse(atan2(y, x) > -pi/2 & atan2(y, x) < pi/2, 0, 1)
      ),
      size = 2
    )
  
  if(!leg){
    gp <- gp + theme(legend.position = "none", plot.margin = margin(1, 1, 1, 1))
  }else{
    gp <- gp + theme(legend.text = element_text(size = 20),
                     legend.title = element_blank(),
                     plot.margin = margin(1, 1, 1, 1))
  }
  return(gp)
}
thre = 2

load(paste0(folder, "Seed_1/results/Test1_Adj_estimation.rda"))
adj_group = G.our.symm.weighted$group
# adj_group[adj_group < thre] = 0
adj_group = log(adj_group)
adj_group[is.na(adj_group)] = 0
adj_group[adj_group> -6 & adj_group < 1] = 0
# adj_group[adj_group <0 & adj_group > -1  ] = 0
plot_single_graph(adj_group, "Differential network")
