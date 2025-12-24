rm(list = ls())
main.path <- "~/Real data analysis"
setwd(main.path)
library(igraph)
load("Separate.RData")

Omega.list <- array(list(),c(3,2))

for(i in 1:3){
  for(j in 1:2){
    Omega.list[i,j] <- as.matrix(x[i,j])
  }
}

channel <- read.csv("~/Real data analysis/channel.txt")
coord <- read.delim("~/Real data analysis/standard_1005_2D.tsv")
coord$label <- toupper(coord$label)
Coord <- coord[which(coord$label %in% channel$channel),] 
Coord <- data.frame(label = Coord$label, X = Coord$x,
                    Y = Coord$y)
addition <- data.frame(label = rbind("X","Y","nd"),
                       X = rbind(0.07635,0.07635,0.07635),
                       Y = rbind(0.85,-0.0808,-0.85))
Coord <- rbind(Coord,addition)
Coord1 <- Coord[match(channel$channel,Coord$label),]
Coord1 <- data.frame(Coord1,location = channel$location)
label <- Coord1$label;Coordinate <- Coord1[,c(2,3)];location <- Coord1$location

layout(matrix(c(1,2,3,4,5,6), nrow = 3, byrow = T), heights = c(1,1,1))
eeg.group <- c("Non Alcoholic","Alcoholic")
method <- c("Separate", "U-MANGO", "H-MANGO")
crit <- 1e-3
for(i in 1:3){
  for(j in 1:2){
    if(i == 1){
      par(mar=c(0,0,6.5,0))
    }else{
      par(mar = c(0,0,0,0))
    }
    connect.map <- array(0,dim(Omega.list[i,j][[1]]))
    p <- length(diag(connect.map))
    colors = c("#3399FF","#33FF33","#FF66CC")
    colnames(connect.map) = rownames(connect.map) <- label
    for(k in 1:p){
      for(l in 1:p){
        connect.map[k,l] <- ifelse(abs(Omega.list[i,j][[1]][k,l])!=0,1,0)
      }
    }
    connect.map <- (t(connect.map) + connect.map)/2
    nonzero <- sum(connect.map - diag(diag(connect.map))!=0)
    
    R <- graph_from_adjacency_matrix(connect.map,mode = "undirected", diag = FALSE,vertex)
    l <- as.matrix(Coordinate)
    if(!is.null(location)){
      V(R)$color <- colors[location]
    }
    V(R)$shape <- 'circle';V(R)$label.cex = 5;V(R)$size = 12
    edge_list <- as_edgelist(R)
    plot(R, layout = l, edge.color = 'black',edge.width = 5,
         rescale = T, main = NULL, margin = 0)
    if(i == 1){
      title(main = eeg.group[j],cex.main = 10)
    }
  }
}
mtext(text = "Separate", side = 3, outer = T, line = -30, cex = 15)
mtext(text = "U-MANGO", side = 3, outer = T, cex = 15, line = -200)
mtext(text = "H-MANGO", side = 3, outer = T, cex = 15, line = -370)