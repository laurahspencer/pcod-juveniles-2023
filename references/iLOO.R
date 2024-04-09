iLOO <- function(sub) {
  
  #estimate sequencing depth and compute cutoff
  sd <- mean(apply(sub,2,sum,na.rm=T))    
  sdcut <- 1/sd
  
  #implement iterative scheme
  outl <- apply(sub, 1, function(z) {
    x <- z
    nbp <- 0
    out <- rep(0,0)
    track <- c(1:length(x))
    
    #iterative scheme
    while((min(nbp,na.rm=T) < sdcut) & (length(x)>2)) {
      
      #build matrix with rows representing leave-one out observation
      mat <- matrix(rep(x,length(x)),ncol=length(x),byrow=T)
      diag(mat) <- NA
      tmp <- t(mat)
      mat <- t(matrix(tmp[!is.na(tmp)],nrow=(length(x)-1),ncol=(length(x))))
      
      #fit negative binomial or Poisson distribution
      nbfit <- apply(mat,1,function(y) {
        if(length(y)>1) {
          v <- var(y)
          m <- mean(y) 
          if(all(y==0)=="TRUE") {
            output <- NA
          } else if (v>m) {
            p <- mean(y)/var(y)
            r <- mean(y)^2/(var(y)-mean(y))
            output <- c(p,r)
          } else {
            lamb <- mean(y)
            output <- c(lamb)
          }
        } else output <- NA
        
        list(output)
      })
      
      nbfit <- lapply(nbfit, "[[", 1)
      #compute probabilities for leave-one out observation  
      nbp <- rep(0,0)
      for (i in 1:length(nbfit)) {
        if(length(nbfit[[i]])==2) {
          nbp <- c(nbp,dnbinom(x[i],prob=nbfit[[i]][1],size=nbfit[[i]][2]))
        } else {
          nbp <- c(nbp,dpois(x[i],lambda=nbfit[[i]][1]))
        }
      }
      
      #compare probabilities to cutoff
      sel <- which(nbp < sdcut)
      if(length(sel)>0) x <- x[-sel]
      if(length(out)==0) {
        out <- c(out,track[sel])
      } else {
        out <- c(out,track[-out][sel])
      }
      
      fout <- rep(NA,length(z))
      if(length(out)>0) {
        fout[out] <- z[out]
      }
    }
    
    list(fout)
  })    
  
  #new data matrix with outliers (all other data is NA’ed)
  identout <- matrix(unlist(outl),nrow=nrow(sub),ncol=ncol(sub),byrow=T)
  colnames(identout) <- colnames(sub)
  rownames(identout) <- rownames(sub)
  
  return(identout)
  
} #end of iLOO function
