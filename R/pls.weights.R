`pls.weights` <-
function(X, IDM, blocks, modes, scheme)
{
    lvs <- nrow(IDM)
    mvs <- ncol(X)
    lvs.names <- rownames(IDM)
    mvs.names <- colnames(X)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    # initialize arbitrary outer weights (value of '1' for all weights)
    out.ws <- as.list(1:lvs)
    for (j in 1:lvs)
        out.ws[[j]] <- rep(1,blocks[j])
    # outer design matrix 'ODM' and matrix of outer weights 'W'
    ODM <- matrix(0, mvs, lvs)
    for (j in 1:lvs)
        ODM[blocklist==j,j] <- rep(1,blocks[j])
    dimnames(ODM) <- list(mvs.names, lvs.names)
    W <- ODM %*% diag(1/sd(X %*% ODM),lvs,lvs)
    w.old <- rowSums(W)    
    w.dif <- 1
    itermax <- 1
    repeat 
    {            
        Y <- X %*% W  # external estimation of LVs 'Y'
        # matrix of inner weights 'e' 
        E <- switch(scheme, 
               "centroid" = sign(cor(Y) * (IDM + t(IDM))),
               "factor" = cor(Y) * (IDM + t(IDM)))
        if (is.null(E)) {   # path weighting scheme
            E <- IDM
            for (j in 1:lvs) {
                if (length(which(IDM[j,]==1)) > 0)
                    E[which(IDM[j,]==1),j] <- lm(Y[,j]~Y[,which(IDM[j,]==1)]-1)$coef
                if (length(which(IDM[,j]==1)) > 0)
                    E[which(IDM[,j]==1),j] <- cor(Y[,j], Y[,which(IDM[,j]==1)])
            } 
        }            
        Z <- Y %*% E  # internal estimation of LVs 'Z'
        Z <- scale(Z)
        # computing outer weights 'w'
        for (j in 1:lvs)
        {
            X.blok = X[,which(blocklist==j)] 
            if (modes[j]=="A")# reflective way
                ODM[which(blocklist==j),j] = cov(Z[,j], X.blok)
            if (modes[j]=="B")# formative way
                ODM[which(blocklist==j),j] = solve.qr(qr(X.blok),Z[,j])
        }
        W <- ODM %*% diag(1/sd(X %*% ODM),lvs,lvs)
        w.new = rowSums(W)                
        w.dif <- sum((w.old - w.new)^2)  # difference of out.weights 
        w.old <- w.new
        if (sum(w.dif^2)<1e-05 || itermax==300) break
        itermax <- itermax + 1
    } # end repeat       
    res.ws <- list(w.new, W)
    if (itermax==300) res.ws=NULL
    return(res.ws)
}

