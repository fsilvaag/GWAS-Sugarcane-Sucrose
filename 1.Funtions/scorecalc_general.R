.score.calc_general <- function(geno,y,Z,X,K,Hinv,ploidy,model,min.MAF,max.geno.freq) {
  m <- ncol(geno)
  n <- nrow(Z)
  scores <- numeric(m)*NA
  beta.out <- scores
  general <- length(grep("general",model,fixed=T))>0 
  P3D <- !is.null(Hinv)
  ef = list()
  for (i in 1:m) {
    S <- .design.score(geno[,i],model,ploidy,min.MAF,max.geno.freq)
    if (!is.null(S)) {
      v1 <- ncol(S)
      X2 <- cbind(X,Z%*%S)
      p <- ncol(X2)
      v2 <- n - p                 
      if (!P3D) {			
        out <- try(mixed.solve(y=y,X=X2,Z=Z,K=K,return.Hinv=TRUE),silent=TRUE)
        if (class(out)!="try-error") { 
          Hinv <- out$Hinv 
        } else {
          Hinv <- NULL
        }
      } 
      
      if (!is.null(Hinv)) {
        W <- crossprod(X2, Hinv %*% X2) 
        Winv <- try(solve(W),silent=TRUE)
        if (!inherits(Winv,what="try-error")) {
          beta <- Winv %*% crossprod(X2, Hinv %*% y)
          resid <- y - X2 %*% beta
          s2 <- as.double(crossprod(resid, Hinv %*% resid))/v2
          Q <- s2 * Winv[(p-v1+1):p,(p-v1+1):p]
          Tt <- solve(Q, silent= TRUE)
          if (!inherits(Tt,what="try-error")) {
            V <- beta[(p+1-v1):p]
            Fstat <- crossprod(V,Tt%*%V)/v1
            x <- v2/(v2+v1*Fstat)
            scores[i] <- -log10(pbeta(x, v2/2, v1/2)) 
            if (!general) {beta.out[i] <- beta[p]} else {
              asd = c(beta[1,],beta[rownames(beta) %in% colnames(S),])
              names(asd)[1] <- paste0("x", min(unique(geno[,i])))
              ef[[colnames(geno)[i]]] <-  asd
              rm(asd);rm(X2)
            }                    
          }
        }
      }
    }
  }
  return(list(score=scores,beta=beta.out, gener = ef))			
}
