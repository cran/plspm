`pls.weights` <-
function(X, IDM, sets, modes, scheme, blocklist)
{
    lvs <- nrow(IDM)
    mvs <- ncol(X)
    lvs.names <- rownames(IDM)
    mvs.names <- colnames(X)
    blocks <- unlist(lapply(sets, length))
    # initialize arbitrary outer weights (value of '1' for all weights)
    out.ws <- sets
    for (k in 1:lvs)
        out.ws[[k]] <- rep(1,length(sets[[k]]))    
    # outer design matrix 'ODM' and matrix of outer weights 'W'
    ODM <- matrix(0, mvs, lvs)
    for (k in 1:lvs)
        ODM[which(blocklist==k),k] <- rep(1,blocks[k])
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
            for (k in 1:lvs) {
                if (length(which(IDM[k,]==1)) > 0)
                    E[which(IDM[k,]==1),k] <- lm(Y[,k]~Y[,which(IDM[k,]==1)]-1)$coef
                if (length(which(IDM[,k]==1)) > 0)
                    E[which(IDM[,k]==1),k] <- cor(Y[,k], Y[,which(IDM[,k]==1)])
            } 
        }            
        Z <- Y %*% E  # internal estimation of LVs 'Z'
        Z <- scale(Z)
        # computing outer weights 'w'
        for (k in 1:lvs)
        {
            X.blok = X[,which(blocklist==k)] 
            if (modes[k]=="A")# reflective way
                ODM[which(blocklist==k),k] = cov(Z[,k], X.blok)
            if (modes[k]=="B")# formative way
                ODM[which(blocklist==k),k] = solve.qr(qr(X.blok),Z[,k])
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
} # end function pls.weights

