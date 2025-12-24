rm(list = ls())
main.path <- "~/Real data analysis"
data.path <- "/EEG datasets"
setwd(paste0(main.path,data.path))
library(igraph)
load("EEG_dataset.RData")

x <- aperm(eegdata$X.mat,c(2,1,3))
x <- x[,31:256,]
y <- eegdata$y.vec
dimen <- dim(x)[1:2] # dimension of x
n <- dim(x)[3] # sample size
K <- 2 # order of tensor

x1 <- x[,,which(y==0)]
x2 <- x[,,which(y==1)]
x1 <- x1 - replicate(dim(x1)[K+1],apply(x1,c(1,2),mean))
x2 <- x2 - replicate(dim(x2)[K+1],apply(x2,c(1,2),mean))

source("~/Separate.R")
source("~/Basic Setting.R")

con1 <- Matrix.Surrogate(x1,type = "Constant");con2 <- Matrix.Surrogate(x2, type = "Constant")
w1 <- sweep(x1,3,con1,"/");w2 <- sweep(x2,3,con2,"/")

# Eventually, I will conduct the precision matrix under H-MANGO.
mat1 <- Matrix.Surrogate(x1,type = "Matrix");mat2 <- Matrix.Surrogate(x2,type = "Matrix")
w.mat1 <- array(0,dim(x1));w.mat2 <- array(0,dim(x2))
for(i in 1:dim(x1)[K+1]){
  eval(parse(text = paste0("w.mat1[",paste0(rep(",",K),collapse = ""),"i] = solve(diag(mat1[,i])) %*% x1[",paste0(rep(",",K),collapse = ""),'i]')))
}
for(i in 1:dim(x2)[K+1]){
  eval(parse(text = paste0("w.mat2[",paste0(rep(",",K),collapse = ""),"i] = solve(diag(mat2[,i])) %*% x2[",paste0(rep(",",K),collapse = ""),'i]')))
}

fit1 <- cv.Separate(x = w.mat1)
fit2 <- cv.Separate(x = w.mat2)
fit.adapt1.matrix <- cv.Separate.penalty(x = w.mat1, Omega = fit1$Omegahat, penalty = "adapt")
fit.adapt2.matrix <- cv.Separate.penalty(x = w.mat2, Omega = fit2$Omegahat, penalty = "adapt")

fit1 <- cv.Separate(x = w1)
fit2 <- cv.Separate(x = w2)
fit.adapt1.constant <- cv.Separate.penalty(x = w1, Omega = fit1$Omegahat, penalty = "adapt")
fit.adapt2.constant <- cv.Separate.penalty(x = w2, Omega = fit2$Omegahat, penalty = "adapt")

fit1 <- cv.Separate(x = x1)
fit2 <- cv.Separate(x = x2)
fit.adapt1.original <- cv.Separate.penalty(x = x1, Omega = fit1$Omegahat, penalty = "adapt")
fit.adapt2.original <- cv.Separate.penalty(x = x2, Omega = fit2$Omegahat, penalty = "adapt")

library(igraph)
Omega.list <- array(list(),c(3,2))
Omega.list[1,1][[1]] <- fit.adapt1.original$Omegahat[[1]];Omega.list[1,2][[1]] <- fit.adapt2.original$Omegahat[[1]]
Omega.list[2,1][[1]] <- fit.adapt1.constant$Omegahat[[1]];Omega.list[2,2][[1]] <- fit.adapt2.constant$Omegahat[[1]]
Omega.list[3,1][[1]] <- fit.adapt1.matrix$Omegahat[[1]];Omega.list[3,2][[1]] <- fit.adapt2.matrix$Omegahat[[1]]

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