library(igraph)
source("/Basic Setting.R")

# Edge weight assigning function
edge.weight <- function(edge_list,Omega,index,crit = 2.5){
  idx.change <- change.index(edge_list = edge_list, index = index)
  weight <- order(Omega[idx.change],decreasing = F) * (crit / nrow(idx.change))
  return(weight)
}

EEG.Graph.plot <- function(Omega.list, label = NULL,Coordinate = NULL,
                           colors = c("#3399FF","#33FF33","#FF66CC"),location = NULL,
                           shape = 'circle', label.cex = 1, margin = -0.1,
                           size = 12, edge_weight = T, rescale = T,
                           set.dir = NULL, plot.feature = NULL, adjust.layout = 1,
                           type = "precision", width.png = 1426,
                           height.png = 704, scale = 5, edge.color = 'black',
                           plot.title = NULL,thr = 1e-3,nonzero.thres = 1e3,
                           plot.main.title = NULL,crit.edge = 5, cex.main = 3){
  idx <- data.frame(label = label, index = seq(1:length(label)))
  q <- length(Omega.list)
  if(!is.list(Omega.list)){
    q <- 1
    if(ncol(Omega.list) != nrow(Omega.list)){
      stop("Omega needs to be squared matrix.")
    }
    p <- length(diag(Omega.list))
    if(type == "precision"){
      connect.map <- array(0,c(p,p));crit <- 0
      if(sum(Omega.list[[i]] - diag(diag(Omega.list[[i]]))!=0) >=nonzero.thres){
        crit <- sort(abs(Omega.list[[i]] - diag(diag(Omega.list[[i]]))),decreasing = T)[nonzero.thres]
      }else{
        crit <- sort(abs(Omega.list[[i]] - diag(diag(Omega.list[[i]]))),decreasing = T)[sum(Omega.list[[i]] - diag(diag(Omega.list[[i]]))!=0)]
      }
      connect.map <- array(0,c(p,p))
      if(is.null(label)){
        colnames(connect.map) = rownames(connect.map) = seq(1:p)
      }else{
        colnames(connect.map) = rownames(connect.map) = label
      }
      for(i in 1:(p-1)){
        for(j in (i+1):p){
          connect.map[i,j] <- ifelse(abs(Omega.list[i,j])>=crit,1,0)
        }
      }
      connect.map <- t(connect.map) + connect.map + diag(p)
      nonzero <- sum(connect.map - diag(diag(connect.map))!=0)
      R <- graph_from_adjacency_matrix(connect.map,mode = "undirected", diag = FALSE, vertex)
      l <- as.matrix(Coordinate)
      if(!is.null(location)){
        V(R)$color <- colors[location]
      }
      V(R)$shape <- shape;V(R)$label.cex = label.cex;V(R)$size = size
      if(!is.null(set.dir)){
        setwd(set.dir)
        if(grepl(".png",plot.title)){
          plot.title1 <- gsub(".png",paste0("_",type,".png"),plot.title)
          png(plot.title1,width = width.png,height = height.png)
          par(mfrow = c(1,q))
          if(edge_weight){
          edge_list <- as_edgelist(R)
          weight <- edge.weight(edge_list = edge_list,Omega.list,index = idx,crit = crit.edge)
          plot(R, layout = l * adjust.layout, edge.color = edge.color,
               rescale = rescale, edge.width = weight,
               main = NULL, margin = margin)
          title(main = plot.main.title, cex.main = cex.main)
          }else{
            plot(R, layout = l * adjust.layout, edge.color = edge.color,
                 rescale = rescale, main = NULL, margin = margin)
            title(main = plot.main.title, cex.main = cex.main)
          }
          dev.off()
        }
      }else{
        par(mfrow = c(1,q))
        if(edge_weight){
        edge_list <- as_edgelist(R)
        weight <- edge.weight(edge_list = edge_list,Omega.list,index = idx,crit = crit.edge)
        plot(R, layout = l * adjust.layout, edge.color = edge.color,
             rescale = rescale, edge.width = weight,
             main = NULL, margin = margin)
        title(main = plot.main.title, cex.main = cex.main)
        }else{
          plot(R, layout = l * adjust.layout, edge.color = edge.color,
               rescale = rescale, main = NULL, margin = margin)
          title(main = plot.main.title, cex.main = cex.main)
        }
      }
    }else{
      A <- cov2cor(solve(Omega.list))
      if(sum(A[[i]] - diag(diag(A[[i]]))!=0) >=nonzero.thres){
        crit[i] <- sort(abs(A[[i]] - diag(diag(A[[i]]))),decreasing = T)[nonzero.thres]
      }else{
        crit[i] <- sort(abs(A[[i]] - diag(diag(A[[i]]))),decreasing = T)[sum(A[[i]] - diag(diag(A[[i]]))!=0)]
      }
      connect.map <- array(0,c(p,p))
      if(is.null(label)){
        colnames(connect.map) = rownames(connect.map) = seq(1:p)
      }else{
        colnames(connect.map) = rownames(connect.map) = label
      }
      for(i in 1:(p-1)){
        for(j in (i+1):p){
          connect.map[i,j] <- ifelse(abs(A[i,j])>=crit,1,0)
        }
      }
      connect.map <- t(connect.map) + connect.map  + diag(p)
      nonzero <- sum(connect.map - diag(diag(connect.map))!=0)
      R <- graph_from_adjacency_matrix(connect.map,mode = "undirected", diag = FALSE, vertex)
      l <- as.matrix(Coordinate)
      if(!is.null(location)){
        V(R)$color <- colors[location]
      }
      V(R)$shape <- shape;V(R)$label.cex = label.cex;V(R)$size = size
      if(!is.null(set.dir)){
        setwd(set.dir)
        if(grepl(".png",plot.title)){
          plot.title1 <- gsub(".png",paste0("_",type,".png"),plot.title)
          png(plot.title1,width = width.png,height = height.png)
          par(mfrow = c(1,q))
          if(edge_weight){
          edge_list <- as_edgelist(R)
          weight <- edge.weight(edge_list = edge_list,A,index = idx,crit = crit.edge)
          plot(R, layout = l * adjust.layout, edge.color = edge.color,
               rescale = rescale, edge.width = weight,
               main = NULL, margin = margin)
          title(main = plot.main.title, cex.main = cex.main)
          dev.off()
          }else{
            plot(R, layout = l * adjust.layout, edge.color = edge.color,
                 rescale = rescale, main = NULL, margin = margin)
            dev.off()
          }
        }
      }else{
        par(mfrow = c(1,q), oma=c(0,0,0,0), mar=c(0,0,0.5,0))
        if(edge_weight){
        edge_list <- as_edgelist(R)
        weight <- edge.weight(edge_list = edge_list,Omega.list,index = idx,crit = crit.edge)
        plot(R, layout = l * adjust.layout, edge.color = edge.color,
             rescale = rescale, edge.width = weight,
             main = NULL, margin = margin)
        title(main = plot.main.title, cex.main = cex.main)
        }else{
          plot(R, layout = l * adjust.layout, edge.color = edge.color,
               rescale = rescale,
               main = NULL, margin = margin)
          title(main = plot.main.title, cex.main = cex.main)
        }
      }
    }
  }else{
    for(i in 1:q){
      if(ncol(Omega.list[[i]])!= nrow(Omega.list[[i]])){
        stop("The element of Omega list needs to be squared matrix.")
      }
    }
    p <- length(diag(Omega.list[[1]]))
    if(type == "precision"){
      connect.map <- array(0,c(p,p,q));crit <- array(NA,q)
      for(i in 1:q){
        if(sum(Omega.list[[i]] - diag(diag(Omega.list[[i]]))!=0) >=nonzero.thres){
          crit[i] <- sort(abs(Omega.list[[i]] - diag(diag(Omega.list[[i]]))),decreasing = T)[nonzero.thres]
        }else{
          crit[i] <- sort(abs(Omega.list[[i]] - diag(diag(Omega.list[[i]]))),decreasing = T)[sum(Omega.list[[i]] - diag(diag(Omega.list[[i]]))!=0)]
        }
        if(is.null(label)){
          if(!is.null(plot.feature)){
            dimnames(connect.map) <- list(seq(1:p),seq(1:p),plot.feature)
          }else{
            dimnames(connect.map) <- list(seq(1:p),seq(1:p),seq(1:q))
          }
        }else{
          if(!is.null(plot.feature)){
            dimnames(connect.map) <- list(label,label,plot.feature)
          }else{
            dimnames(connect.map) <- list(label,label,seq(1:q))
          }
        }
      }
      nonzero <- array(NA,q)
      for(k in 1:q){
        for(i in 1:(p-1)){
          for(j in (i+1):p){
            connect.map[i,j,k] <- ifelse(abs(Omega.list[[k]][i,j])>=crit[k],1,0)
          }
        }
        connect.map[,,k] <- t(connect.map[,,k]) + connect.map[,,k] + diag(p)
        nonzero[k] <- sum(connect.map[,,k] - diag(diag(connect.map[,,k]))!=0)
      }
      if(!is.null(set.dir)){
        setwd(set.dir)
        if(grepl(".png",plot.title)){
          plot.title1 <- gsub(".png",paste0("_",type,".png"),plot.title)
          png(plot.title1,width = width.png,height = height.png)
          par(mfrow = c(1,q))
          for(k in 1:q){
            R <- graph_from_adjacency_matrix(connect.map[,,k],mode = "undirected", diag = FALSE, vertex)
            l <- as.matrix(Coordinate)
            if(!is.null(location)){
              V(R)$color <- colors[location]
            }
            V(R)$shape <- shape;V(R)$label.cex = label.cex;V(R)$size = size
            if(edge_weight){
            edge_list <- as_edgelist(R)
            weight <- edge.weight(edge_list = edge_list,Omega.list[[k]],index = idx,crit = crit.edge)
            plot(R, layout = l * adjust.layout, edge.color = edge.color,
                 rescale = rescale, edge.width = weight,
                 main = NULL, margin = margin)
            title(main = plot.main.title[k], cex.main = cex.main)
            
            }else{
              plot(R, layout = l * adjust.layout, edge.color = edge.color,
                   rescale = rescale, 
                   main = NULL, margin = margin)
              title(main = plot.main.title[k], cex.main = cex.main)
            }
          }
          dev.off()
        }
      }else{
        par(mfrow = c(1,q),  oma=c(0,0,0,0), mar=c(0,0,0.5,0))
        for(k in 1:q){
          R <- graph_from_adjacency_matrix(connect.map[,,k],mode = "undirected", diag = FALSE, vertex)
          l <- as.matrix(Coordinate)
          if(!is.null(location)){
            V(R)$color <- colors[location]
          }
          V(R)$shape <- shape;V(R)$label.cex = label.cex;V(R)$size = size
          edge_list <- as_edgelist(R)
          if(edge_weight){
          weight <- edge.weight(edge_list = edge_list,Omega.list[[k]],index = idx,crit = crit.edge)
          plot(R, layout = l * adjust.layout, edge.color = edge.color,
               rescale = rescale, edge.width = weight,
               main = NULL, margin = margin)
          title(main = plot.main.title[k], cex.main = cex.main)
          }else{
            plot(R, layout = l * adjust.layout, edge.color = edge.color,
                 rescale = rescale, 
                 main = NULL, margin = margin)
            title(main = plot.main.title[k], cex.main = cex.main)
          }
        }
      }
    }else{
      connect.map <- array(0,c(p,p,q));crit <- array(NA,q)
      A <- array(list(),q)
      
      for(i in 1:q){
        A[[i]] <- cov2cor(solve(Omega.list[[i]]))
        if(sum(A[[i]] - diag(diag(A[[i]]))!=0) >=nonzero.thres){
          crit[i] <- sort(abs(A[[i]] - diag(diag(A[[i]]))),decreasing = T)[nonzero.thres]
        }else{
          crit[i] <- sort(abs(A[[i]] - diag(diag(A[[i]]))),decreasing = T)[sum(A[[i]] - diag(diag(A[[i]]))!=0)]
        }
        if(is.null(label)){
          if(!is.null(plot.feature)){
            dimnames(connect.map) <- list(seq(1:p),seq(1:p),plot.feature)
          }else{
            dimnames(connect.map) <- list(seq(1:p),seq(1:p),seq(1:q))
          }
        }else{
          if(!is.null(plot.feature)){
            dimnames(connect.map) <- list(label,label,plot.feature)
          }else{
            dimnames(connect.map) <- list(label,label,seq(1:q))
          }
        }
      }
      nonzero <- array(NA,q)
      for(k in 1:q){
        for(i in 1:(p-1)){
          for(j in (i+1):p){
            connect.map[i,j,k] <- ifelse(abs(A[[k]][i,j])>=crit[k],1,0)
          }
        }
        connect.map[,,k] <- t(connect.map[,,k]) + connect.map[,,k] + diag(p)
        nonzero[k] <- sum(connect.map[,,k] - diag(diag(connect.map[,,k]))!=0)
      }
      if(!is.null(set.dir)){
        setwd(set.dir)
        if(grepl(".png",plot.title)){
          plot.title1 <- gsub(".png",paste0("_",type,".png"),plot.title)
          png(plot.title1,width = width.png,height = height.png)
          par(mfrow = c(1,q))
          for(k in 1:q){
            R <- graph_from_adjacency_matrix(connect.map[,,k],mode = "undirected", diag = FALSE, vertex)
            l <- as.matrix(Coordinate)
            if(!is.null(location)){
              V(R)$color <- colors[location]
            }
            V(R)$shape <- shape;V(R)$label.cex = label.cex;V(R)$size = size
            if(edge_weight){
            edge_list <- as_edgelist(R)
            weight <- edge.weight(edge_list = edge_list,A[[k]],index = idx,crit = crit.edge)
            plot(R, layout = l * adjust.layout, edge.color = edge.color,
                 rescale = rescale, edge.width = weight,
                 main = NULL, margin = margin)
            title(main = plot.main.title[k], cex.main = cex.main)
            }else{
              plot(R, layout = l * adjust.layout, edge.color = edge.color,
                   rescale = rescale,
                   main = NULL, margin = margin)
              title(main = plot.main.title[k], cex.main = cex.main)
            }
          }
          dev.off()
        }
      }else{
        par(mfrow = c(1,q),  oma=c(0,0,0,0), mar=c(0,0,0.5,0))
        for(k in 1:q){
          R <- graph_from_adjacency_matrix(connect.map[,,k],mode = "undirected", diag = FALSE, vertex)
          l <- as.matrix(Coordinate)
          if(!is.null(location)){
            V(R)$color <- colors[location]
          }
          V(R)$shape <- shape;V(R)$label.cex = label.cex;V(R)$size = size
          if(edge_weight){
          edge_list <- as_edgelist(R)
          weight <- edge.weight(edge_list = edge_list,A[[k]],index = idx,crit = crit.edge)
          plot(R, layout = l * adjust.layout, edge.color = edge.color,
               rescale = rescale, edge.width = weight,
               main = NULL, margin = margin)
          title(main = plot.main.title[k], cex.main = cex.main)
          }else{
            plot(R, layout = l * adjust.layout, edge.color = edge.color,
                 rescale = rescale, main = NULL, margin = margin)
            title(main = plot.main.title[k], cex.main = cex.main)
          }
        }
      }
    }
  }
  if(edge_weight){
    return(weight)
  }
}

nonzero.edge.count <- function(Omega){
  p <- length(diag(Omega))
  Omega1 <- Omega - diag(diag(Omega))
  connect.map <- array(0,c(p,p))
  for(i in 1:p){
    for(j in 1:p){
      connect.map[i,j] <- ifelse(Omega1[i,j]!=0,1,0)
    }
  }
  if(norm(connect.map - t(connect.map),"M")==0){
    count <- sum(connect.map!=0)/2
  }else{
    count1 <-length(intersect(which(connect.map == t(connect.map)) ,which(connect.map !=0)))/2
    count2 <- nrow(index(which(connect.map != t(connect.map)),dim = c(p,p))) / 2
    count <- count1 + count2
  }
  return(count)
}

