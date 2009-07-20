`plsr1` <-
function(x, y, nc=2, scaled=TRUE)
{
    # ============ checking arguments ============
    X <- as.matrix(x)
    Y <- as.matrix(y)
    n <- nrow(X)
    p <- ncol(X)
    # ============ setting inputs ==============
    if (scaled) Xx<-scale(X) else Xx<-scale(X,scale=F)
    if (scaled) Yy<-scale(Y) else Yy<-scale(Y,scale=F)
    X.old <- Xx
    Y.old <- Yy
    Th <- matrix(NA, n, nc)# matrix of X-scores
    Ph <- matrix(NA, p, nc)# matrix of X-loadings
    Wh <- matrix(NA, p, nc)# matrix of raw-weights
    Uh <- matrix(NA, n, nc)# matrix of Y-scores
    ch <- rep(NA, nc)# vector of y-loadings
    # ============ pls regression algorithm ==============
    for (h in 1:nc)
    {
        w.old <- t(X.old) %*% Y.old / sum(Y.old^2)
        w.new <- w.old / sqrt(sum(w.old^2)) # normalization
        t.new <- X.old %*% w.new
        p.new <- t(X.old) %*% t.new / sum(t.new^2) 
        c.new <- t(Y.old) %*% t.new / sum(t.new^2)
        u.new <- Y.old / as.vector(c.new)
        Y.old <- Y.old - t.new%*%c.new# deflate y.old
        X.old <- X.old - (t.new %*% t(p.new))# deflate X.old
        Th[,h] <- round(t.new, 4)
        Ph[,h] <- round(p.new, 4)
        Wh[,h] <- round(w.new, 4)
        Uh[,h] <- round(u.new, 4)
        ch[h] <- round(c.new, 4)        
    }
    Ws <- round(Wh %*% solve(t(Ph)%*%Wh), 4)# modified weights
    Bs <- round(as.vector(Ws %*% ch), 4) # std beta coeffs    
    Br <- round(Bs * (rep(sd(Y),p)/apply(X,2,sd)), 4)   # beta coeffs
    cte <- as.vector(round(mean(y) - Br%*%apply(X,2,mean), 4))# intercept
    y.hat <- round(X%*%Br+cte, 4)# y predicted
    resid <- round(as.vector(Y - y.hat), 4)# residuals
    R2 <- round(as.vector(cor(Th, Yy))^2, 4)  # R2 coefficients    
    names(Br) <- colnames(X)
    names(resid) <- rownames(Y)
    names(y.hat) <- rownames(Y)
    names(R2) <- paste(rep("t",nc),1:nc,sep="")
    res <- list(coeffs=Br, cte=cte, R2=R2[1:nc], resid=resid, y.pred=y.hat)    
    return(res)
}

