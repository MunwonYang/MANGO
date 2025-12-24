source("~/Heavy_tail/Basic Setting.R")
Tlasso.edit <- function(x,val,T = 5, nlambda = c(50,20), 
                        lambda.min.ratio = 1e-4,
                        norm.type = 2, thres = 1e-5,
                        scale = NULL, thr = 1e-3,
                        maxit = 1e3){
  if (!is.array(x)){
    stop('argument data should be an array')
  } 
  
  dim.vec = dim(x) # m1 , m2 , ... , mK , n
  K = length(dim.vec) - 1
  if(is.null(scale)){
    scale <- rep(1,K)
  }
  
  if (  length(dim.vec)<2   ) {
    stop('argument data should be at least a 2-dimension array')
  }
  
  n = as.integer(dim.vec[K+1])
  m.vec = dim.vec[1:K]
  
  n.val <- as.integer(dim(val)[K+1])
  for(i in 1:n.val){
    d=0
    eval(parse(text=paste('d=val[',paste(rep(',',K),collapse=''),'i]')))
    if(sum(dim(d)!=m.vec)!=0){
      stop("Dimension of Validation data is different with Training data.")
    }
  }
  
  if (!(norm.type==1 | norm.type==2)){
    stop('argument norm.type should be 1 or 2')
  } else if (!( (T==round(T)) & ((T>1)|(T==1) )) ) {
    stop('argument T should be a positive integer')
  } 
  Omega.list=list()
  for (k in 1:K) {
    Omega.list[[k]] = diag(m.vec[k])
  }
  Omega.list.sqrt=Omega.list
  Lambda.list <- array(list(),K)
  for(k in 1:K){
    S.array = array(0,c(m.vec[k],m.vec[k],n))
    for(i in 1:n){
      # set k-th matrix into identity
      Omega.list.sqrt[[k]]=diag(m.vec[k])
      # the ith observation
      
      d=0
      eval(parse(text=paste('d=x[',paste(rep(',',K),collapse=''),'i]')))
      
      Vi = k_unfold( as.tensor(ttl( as.tensor(d) , Omega.list.sqrt , ms=1:K)@data) ,m=k)@data
      S.array[,,i] = Vi %*% t(Vi)
    }
    S.mat = apply(S.array,c(1,2),mean) * m.vec[k] / prod(m.vec)
    lambda.max <- max(max(S.mat - diag(m.vec[k])),-min(S.mat - diag(m.vec[k])))
    lambda.min <- lambda.min.ratio* lambda.max
    Lambda.list[[k]] <- exp(seq(log(lambda.min),log(lambda.max), length.out = nlambda[1]))
  }
  count <- 0
  for (iter in 1:T) {
    # update from 1st to Kth mode
    Omega.list.old=Omega.list # store an old list for termination after current iteration
    Loglik <- array(list(),K)
    Lambda <- array(0,K)
    check <- array(NA,K)
    for(k in 1:K){
      # compute the estimation of covariance matrix to construct likelihood
      loglik <- array(0,length(Lambda.list[[k]]))
      S.array = array(0,c(m.vec[k],m.vec[k],n))
      for(i in 1:n){
        # set k-th matrix into identity
        Omega.list.sqrt[[k]]=diag(m.vec[k])
        # the ith observation
        
        d=0
        eval(parse(text=paste('d=x[',paste(rep(',',K),collapse=''),'i]')))
        
        Vi = k_unfold( as.tensor(ttl( as.tensor(d) , Omega.list.sqrt , ms=1:K)@data) ,m=k)@data
        S.array[,,i] = Vi %*% t(Vi)
      }
      S.mat = apply(S.array,c(1,2),mean) * m.vec[k] / prod(m.vec)
      S.array = array(0,c(m.vec[k],m.vec[k],n.val))
      for(i in 1:n.val){
        # set k-th matrix into identity
        Omega.list.sqrt[[k]]=diag(m.vec[k])
        # the ith observation
        
        d=0
        eval(parse(text=paste('d=val[',paste(rep(',',K),collapse=''),'i]')))
        
        Vi = k_unfold( as.tensor(ttl( as.tensor(d) , Omega.list.sqrt , ms=1:K)@data) ,m=k)@data
        S.array[,,i] = Vi %*% t(Vi)
      }
      S.val = apply(S.array,c(1,2),mean) * m.vec[k] / prod(m.vec)
      for(i in 1:nlambda[1]){
        Out1 <- glasso(S.mat,rho = Lambda.list[[k]][i],
                       penalize.diagonal = FALSE, thr = thr,
                       maxit = maxit)$wi
        calc <- log(eigen(Out1)$values)
        loglik[i] <- -sum(S.val * Out1) + sum(calc)
      }
      ind <- which.max(loglik)
      if(ind != 1 | ind != nlambda[1]){
        check[k] <- TRUE
      }
      Loglik[[k]] <- loglik
      Lambda[k] <- Lambda.list[[k]][ind]
      Out.result <- glasso(S.mat, rho = Lambda[k],
                           penalize.diagonal = FALSE, thr = thr, 
                           maxit = maxit)$wi
      # normalize matrix and stored into Omega.list
      if(norm.type == 2){
        Omega.list[[k]] = as.matrix(Out.result) / norm(as.matrix(Out.result),type='F')
      }else if(norm.type == 1){
        Omega.list[[k]] = as.matrix(Out.result) / as.matrix(Out.result)[1,1]
      }
      
      # store the square root of kth precision matrix for likelihood in next update
      Omega.list.sqrt[[k]]=sqrtm(Omega.list[[k]])
    }
    diff=0
    for (i in 1:K){
      diff=diff+norm(Omega.list.old[[i]]-Omega.list[[i]],type='F')
    }
    if ( diff < thres ) {
      print('Output changes less than arguement thres after certain iteration. Terminate the algorithm.')
    }
    count <- count+1
  }
  
  Omega.list=list()
  for (k in 1:K) {
    Omega.list[[k]] = diag(m.vec[k])
  }
  Omega.list.sqrt=Omega.list
  Lambda.list1 <- optimal.lambda.list(Lambda.list,Loglik,nlambda = nlambda[2])
  count <- 0
  for (iter in 1:T) {
    # update from 1st to Kth mode
    Omega.list.old=Omega.list # store an old list for termination after current iteration
    Loglik1 <- array(list(),K)
    Lambda <- array(0,K)
    check <- array(NA,K)
    for(k in 1:K){
      # compute the estimation of covariance matrix to construct likelihood
      loglik <- array(0,nlambda[2])
      S.array = array(0,c(m.vec[k],m.vec[k],n))
      for(i in 1:n){
        # set k-th matrix into identity
        Omega.list.sqrt[[k]]=diag(m.vec[k])
        # the ith observation
        
        d=0
        eval(parse(text=paste('d=x[',paste(rep(',',K),collapse=''),'i]')))
        
        Vi = k_unfold( as.tensor(ttl( as.tensor(d) , Omega.list.sqrt , ms=1:K)@data) ,m=k)@data
        S.array[,,i] = Vi %*% t(Vi)
      }
      S.mat = apply(S.array,c(1,2),mean) * m.vec[k] / prod(m.vec)
      S.array = array(0,c(m.vec[k],m.vec[k],n.val))
      for(i in 1:n.val){
        # set k-th matrix into identity
        Omega.list.sqrt[[k]]=diag(m.vec[k])
        # the ith observation
        
        d=0
        eval(parse(text=paste('d=val[',paste(rep(',',K),collapse=''),'i]')))
        
        Vi = k_unfold( as.tensor(ttl( as.tensor(d) , Omega.list.sqrt , ms=1:K)@data) ,m=k)@data
        S.array[,,i] = Vi %*% t(Vi)
      }
      S.val = apply(S.array,c(1,2),mean) * m.vec[k] / prod(m.vec)
      for(i in 1:nlambda[2]){
        Out1 <- glasso(S.mat,rho = Lambda.list1[[k]][i],
                       penalize.diagonal = FALSE, thr = thr,
                       maxit = maxit)$wi
        calc <- log(eigen(Out1)$values)
        loglik[i] <- -sum(S.val * Out1) + sum(calc)
      }
      ind <- which.max(loglik)
      if(ind != 1 | ind != nlambda[2]){
        check[k] <- TRUE
      }
      Loglik1[[k]] <- loglik
      Lambda[k] <- Lambda.list1[[k]][ind]
      Out.result <- glasso(S.mat, rho = Lambda[k],
                           penalize.diagonal = FALSE, thr = thr, 
                           maxit = maxit)$wi
      # normalize matrix and stored into Omega.list
      if(norm.type == 2){
        Omega.list[[k]] = as.matrix(Out.result) / norm(as.matrix(Out.result),type='F')
      }else if(norm.type == 1){
        Omega.list[[k]] = as.matrix(Out.result) / as.matrix(Out.result)[1,1]
      }
      
      # store the square root of kth precision matrix for likelihood in next update
      Omega.list.sqrt[[k]]=sqrtm(Omega.list[[k]])
    }
    diff=0
    for (i in 1:K){
      diff=diff+norm(Omega.list.old[[i]]-Omega.list[[i]],type='F')
    }
    if ( diff < thres ) {
      print('Output changes less than arguement thres after certain iteration. Terminate the algorithm.')
      result <- list()
      result$Omegahat<- Omega.list
      result$iteration <- count
      result$explore.loglik <- Loglik
      result$final.loglik <- Loglik1
      result$explore.lambda <- Lambda.list
      result$final.lambda <- Lambda.list1
      result$Lambda <- Lambda
      return(result)
    }
    count <- count+1
  }  
  Result <- list()
  Result$Omegahat <- Omega.list
  Result$iteration <- count
  Result$explore.loglik <- Loglik
  Result$final.loglik <- Loglik1
  Result$explore.lambda <- Lambda.list
  Result$final.lambda <- Lambda.list1
  Result$lambda <- Lambda
  return(Result)
}

Tlasso.adapt <- function(x,val,Omega,gam = 1, norm.type = 2,
                         rate = 1e-5, T = 5,
                         nlambda = c(20,10),penalty = "adapt",
                         thr = 1e-3,thres = 1e-5,maxit = 1e3){
  if (!is.array(x)){
    stop('argument data should be an array')
  } 
  
  dim.vec = dim(x) # m1 , m2 , ... , mK , n
  K = length(dim.vec) - 1
  if(is.null(scale)){
    scale <- rep(1,K)
  }
  
  if (  length(dim.vec)<2   ) {
    stop('argument data should be at least a 2-dimension array')
  }
  
  n = as.integer(dim.vec[K+1])
  m.vec = dim.vec[1:K]
  
  n.val <- as.integer(dim(val)[K+1])
  for(i in 1:n.val){
    d=0
    eval(parse(text=paste('d=val[',paste(rep(',',K),collapse=''),'i]')))
    if(sum(dim(d)!=m.vec)!=0){
      stop("Dimension of Validation data is different with Training data.")
    }
  }
  
  if (!(norm.type==1 | norm.type==2)){
    stop('argument norm.type should be 1 or 2')
  } else if (!( (T==round(T)) & ((T>1)|(T==1) )) ) {
    stop('argument T should be a positive integer')
  } 
  weight.matrix <- array(list(),K)
  Lambda.list <- array(list(),K)
  Omega.list=Omega
  Omega.list.sqrt <- array(list(),K)
  for (k in 1:K) {
    Omega.list.sqrt[[k]] = sqrtm(Omega.list[[k]])
  }
  for(k in 1:K){
    Theta <- Omega[[k]]
    weight.matrix[[k]] <- eval(parse(
      text =  paste0(
        penalty,
        "_deriv(Theta = Theta, gam = gam)"
      )
    ))
    S.array = array(0,c(m.vec[k],m.vec[k],n))
    for(i in 1:n){
      # set k-th matrix into identity
      Omega.list.sqrt[[k]]=diag(m.vec[k])
      # the ith observation
      d=0
      eval(parse(text=paste('d=x[',paste(rep(',',K),collapse=''),'i]')))
      
      Vi = k_unfold( as.tensor(ttl( as.tensor(d) , Omega.list.sqrt , ms=1:K)@data) ,m=k)@data
      S.array[,,i] = Vi %*% t(Vi)
    }
    S.mat = apply(S.array,c(1,2),mean) * m.vec[k] / prod(m.vec)
    lambda.max <- max(max(S.mat - diag(m.vec[k])),-min(S.mat - diag(m.vec[k])))
    lambda.min <- rate * lambda.max
    Lambda.list[[k]] <- exp(seq(log(lambda.min),log(lambda.max),length.out = nlambda[1]))
    }
    for(iter in 1:T){
      Omega.list.old=Omega.list # store an old list for termination after current iteration
      Loglik1 <- array(list(),K)
      Lambda <- array(0,K)
      check <- array(NA,K)
      for(k in 1:K){
        # compute the estimation of covariance matrix to construct likelihood
        loglik <- array(0,length(Lambda.list[[k]]))
        S.array = array(0,c(m.vec[k],m.vec[k],n))
        for(i in 1:n){
          # set k-th matrix into identity
          Omega.list.sqrt[[k]]=diag(m.vec[k])
          # the ith observation
          
          d=0
          eval(parse(text=paste('d=x[',paste(rep(',',K),collapse=''),'i]')))
          
          Vi = k_unfold( as.tensor(ttl( as.tensor(d) , Omega.list.sqrt , ms=1:K)@data) ,m=k)@data
          S.array[,,i] = Vi %*% t(Vi)
        }
        S.mat = apply(S.array,c(1,2),mean) * m.vec[k] / prod(m.vec)
        S.array = array(0,c(m.vec[k],m.vec[k],n.val))
        for(i in 1:n.val){
          # set k-th matrix into identity
          Omega.list.sqrt[[k]]=diag(m.vec[k])
          # the ith observation
          
          d=0
          eval(parse(text=paste('d=val[',paste(rep(',',K),collapse=''),'i]')))
          
          Vi = k_unfold( as.tensor(ttl( as.tensor(d) , Omega.list.sqrt , ms=1:K)@data) ,m=k)@data
          S.array[,,i] = Vi %*% t(Vi)
        }
        S.val = apply(S.array,c(1,2),mean) * m.vec[k] / prod(m.vec)
        for(i in 1:nlambda[1]){
          Out1 <- glasso(S.mat,rho = Lambda.list[[k]][i] * weight.matrix[[k]],
                         penalize.diagonal = FALSE, thr = thr,
                         maxit = maxit)$wi
          calc <- log(eigen(Out1)$values)
          loglik[i] <- -sum(S.val * Out1) + sum(calc)
        }
        ind <- which.max(loglik)
        if(ind != 1 | ind != nlambda[1]){
          check[k] <- TRUE
        }
        Loglik1[[k]] <- loglik
        Lambda[k] <- Lambda.list[[k]][ind]
        Out.result <- glasso(S.mat, rho = Lambda[k],
                             penalize.diagonal = FALSE, thr = thr, 
                             maxit = maxit)$wi
        # normalize matrix and stored into Omega.list
        if(norm.type == 2){
          Omega.list[[k]] = as.matrix(Out.result) / norm(as.matrix(Out.result),type='F')
        }else if(norm.type == 1){
          Omega.list[[k]] = as.matrix(Out.result) / as.matrix(Out.result)[1,1]
        }
        
        # store the square root of kth precision matrix for likelihood in next update
        Omega.list.sqrt[[k]]=sqrtm(Omega.list[[k]])
      }
      diff=0
      for (i in 1:K){
        diff=diff+norm(Omega.list.old[[i]]-Omega.list[[i]],type='F')
      }
      if ( diff < thres ) {
        print('Output changes less than arguement thres after certain iteration. Terminate the algorithm.')
      }
    }
    Lambda.list1 <- optimal.lambda.list(lambda = Lambda.list,
                                        loglik = Loglik1,
                                        nlambda = nlambda[2])
    Omega.list=Omega
    Omega.list.sqrt <- array(list(),K)
    for (k in 1:K) {
      Omega.list.sqrt[[k]] = sqrtm(Omega.list[[k]])
    }
    count <- 0
    for (iter in 1:T) {
      # update from 1st to Kth mode
      Omega.list.old=Omega.list # store an old list for termination after current iteration
      Loglik2 <- array(list(),K)
      Lambda <- array(0,K)
      check <- array(NA,K)
      for(k in 1:K){
        # compute the estimation of covariance matrix to construct likelihood
        loglik <- array(0,length(Lambda.list1[[k]]))
        S.array = array(0,c(m.vec[k],m.vec[k],n))
        for(i in 1:n){
          # set k-th matrix into identity
          Omega.list.sqrt[[k]]=diag(m.vec[k])
          # the ith observation
          
          d=0
          eval(parse(text=paste('d=x[',paste(rep(',',K),collapse=''),'i]')))
          
          Vi = k_unfold( as.tensor(ttl( as.tensor(d) , Omega.list.sqrt , ms=1:K)@data) ,m=k)@data
          S.array[,,i] = Vi %*% t(Vi)
        }
        S.mat = apply(S.array,c(1,2),mean) * m.vec[k] / prod(m.vec)
        S.array = array(0,c(m.vec[k],m.vec[k],n.val))
        for(i in 1:n.val){
          # set k-th matrix into identity
          Omega.list.sqrt[[k]]=diag(m.vec[k])
          # the ith observation
          
          d=0
          eval(parse(text=paste('d=val[',paste(rep(',',K),collapse=''),'i]')))
          
          Vi = k_unfold( as.tensor(ttl( as.tensor(d) , Omega.list.sqrt , ms=1:K)@data) ,m=k)@data
          S.array[,,i] = Vi %*% t(Vi)
        }
        S.val = apply(S.array,c(1,2),mean) * m.vec[k] / prod(m.vec)
        for(i in 1:nlambda[2]){
          Out1 <- glasso(S.mat,rho = Lambda.list1[[k]][i] * weight.matrix[[k]],
                         penalize.diagonal = FALSE, thr = thr,
                         maxit = maxit)$wi
          calc <- log(eigen(Out1)$values)
          loglik[i] <- -sum(S.val * Out1) + sum(calc)
        }
        ind <- which.max(loglik)
        if(ind != 1 | ind != nlambda[2]){
          check[k] <- TRUE
        }
        Loglik2[[k]] <- loglik
        Lambda[k] <- Lambda.list1[[k]][ind]
        Out.result <- glasso(S.mat, rho = Lambda[k] * weight.matrix[[k]],
                             penalize.diagonal = FALSE, thr = thr, 
                             maxit = maxit)$wi
        # normalize matrix and stored into Omega.list
        if(norm.type == 2){
          Omega.list[[k]] = as.matrix(Out.result) / norm(as.matrix(Out.result),type='F')
        }else if(norm.type == 1){
          Omega.list[[k]] = as.matrix(Out.result) / as.matrix(Out.result)[1,1]
        }
        
        # store the square root of kth precision matrix for likelihood in next update
        Omega.list.sqrt[[k]]=sqrtm(Omega.list[[k]])
      }
      diff=0
      for (i in 1:K){
        diff=diff+norm(Omega.list.old[[i]]-Omega.list[[i]],type='F')
      }
      if ( diff < thres ) {
        print('Output changes less than arguement thres after certain iteration. Terminate the algorithm.')
        result <- list()
        result$Omegahat <- Omega.list
        result$iteration <- count
        result$explore.loglik <- Loglik1
        result$final.loglik <- Loglik2
        result$explore.lambda <- Lambda.list
        result$final.lambda <- Lambda.list1
        result$Lambda <- Lambda
        result$weight <- weight.matrix
        return(result)
      }
      count <- count+1
    }  
    Result <- list()
    Result$Omegahat <- Omega.list
    Result$iteration <- count
    Result$explore.loglik <- Loglik1
    Result$final.loglik <- Loglik2
    Result$explore.lambda <- Lambda.list
    Result$final.lambda <- Lambda.list1
    Result$Lambda <- Lambda
    Result$weight <- weight.matrix
    return(Result)
}
