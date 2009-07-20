`pls.boot` <-
function(DM, IDM, blocks, modes, scheme, scaled, br, plsr)
{
    lvs <- nrow(IDM)
    lvs.names <- rownames(IDM)
    mvs <- ncol(DM)
    mvs.names <- colnames(DM)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1    
    bootnum <- br
    if(scaled) X=scale(DM) else X=scale(DM, scale=FALSE)
    # =============== computation of the original plspm model ================
    out.ws <- pls.weights(X, IDM, blocks, modes, scheme)
    wgs.orig <- out.ws[[1]]
    cor.XY <- cor(X, X%*%out.ws[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
         w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    Y.lvs <- X %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs)
    pathmod <- pls.paths(IDM, Y.lvs, plsr)
    Path <- pathmod[[2]]
    path.orig <- as.vector(Path[which(IDM==1)])
    r2.orig <- pathmod[[3]][which(endo==1)]
    Path.efs <- pls.efects(Path)
    loadcomu <- pls.loads(X, Y.lvs, blocks)    
    loads.orig <- loadcomu[[1]]
    # ========================= Bootstrap Validation =========================
    path.labs <- NULL
    efs.labs <- NULL
    for (j in 1:lvs)
        for (i in j:lvs)
             if (IDM[i,j]==1) 
                 path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))    
    WEIGS <- matrix(NA, bootnum, mvs)
    LOADS <- matrix(NA, bootnum, mvs)
    PATHS <- matrix(NA, bootnum, sum(IDM))
    TOEFS <- matrix(NA, bootnum, nrow(Path.efs))
    RSQRS <- matrix(NA, bootnum, sum(endo))
    for (i in 1:bootnum)
    {
        boot.obs <- sample(1:nrow(X), nrow(X), replace=TRUE)
        X.boot <- X[boot.obs,]
        w.boot <- pls.weights(X.boot, IDM, blocks, modes, scheme)
        if (is.null(w.boot)) stop("Bootstrapping failed") 
        WEIGS[i,] <- w.boot[[1]]
        Y.boot <- X.boot %*% w.boot[[2]]
        pathmod <- pls.paths(IDM, Y.boot, plsr)
        P.boot <- pathmod[[2]]
        Toef.boot <- pls.efects(P.boot)
        PATHS[i,] <- as.vector(P.boot[which(IDM==1)])
        TOEFS[i,] <- Toef.boot[,4]
        RSQRS[i,] <- pathmod[[3]][which(endo==1)]
        l.boot <- pls.loads(X.boot, Y.boot, blocks)    
        LOADS[i,] <- l.boot[[1]]
    }
    # Outer weights
    colnames(WEIGS) <- mvs.names
    t.wb <- matrix(NA, mvs, 2)
    for (j in 1:ncol(WEIGS))
        t.wb[j,] <- round(unlist(t.test(WEIGS[,j],mu=wgs.orig[j])[c(1,3)]), 4)
    WB <- data.frame(Original=round(wgs.orig,4), Mean.Boot=round(apply(WEIGS,2,mean), 4), 
                Std.Err=round(apply(WEIGS,2,sd),4), t.statis=t.wb[,1], p.value=t.wb[,2],
                perc.025=round(apply(WEIGS,2,mean) - qt(.975,nrow(X))*apply(WEIGS,2,sd), 4), 
                perc.975=round(apply(WEIGS,2,mean) + qt(.975,nrow(X))*apply(WEIGS,2,sd), 4))
    # Loadings
    colnames(LOADS) <- mvs.names
    t.lb <- matrix(NA, mvs, 2)
    for (j in 1:ncol(LOADS)) {
        if (loads.orig[j]!=1)
            t.lb[j,] <- round(unlist(t.test(LOADS[,j],mu=loads.orig[j])[c(1,3)]), 4)
        else
            t.lb[j,] <- c(0,1)
    }
    LB <- data.frame(Original=round(loads.orig,4), Mean.Boot=round(apply(LOADS,2,mean), 4), 
                Std.Err=round(apply(LOADS,2,sd),4), t.statis=t.lb[,1], p.value=t.lb[,2],
                perc.025=round(apply(LOADS,2,mean) - qt(.975,nrow(X))*apply(LOADS,2,sd), 4), 
                perc.975=round(apply(LOADS,2,mean) + qt(.975,nrow(X))*apply(LOADS,2,sd), 4))
    # Path Coefficients
    colnames(PATHS) <- path.labs
    t.pb <- matrix(NA, sum(IDM), 2)
    for (j in 1:ncol(PATHS))
        t.pb[j,] <- round(unlist(t.test(PATHS[,j],mu=path.orig[j])[c(1,3)]), 4)
    PB <- data.frame(Original=round(path.orig,4), Mean.Boot=round(apply(PATHS,2,mean), 4), 
                Std.Err=round(apply(PATHS,2,sd),4), t.statis=t.pb[,1], p.value=t.pb[,2],
                perc.025=round(apply(PATHS,2,mean) - qt(.975,nrow(X))*apply(PATHS,2,sd), 4), 
                perc.975=round(apply(PATHS,2,mean) + qt(.975,nrow(X))*apply(PATHS,2,sd), 4))
    # Total Effects
    colnames(TOEFS) <- Path.efs[,1]
    t.teb <- matrix(NA, nrow(Path.efs), 2)
    for (j in 1:ncol(TOEFS))
        t.teb[j,] <- round(unlist(t.test(TOEFS[,j],mu=Path.efs[j,4])[c(1,3)]), 4)
    TE <- data.frame(Original=round(Path.efs[,4],4), Mean.Boot=round(apply(TOEFS,2,mean), 4), 
                Std.Err=round(apply(TOEFS,2,sd), 4), t.statis=t.teb[,1], p.value=t.teb[,2],
                perc.025=round(apply(TOEFS,2,mean) - qt(.975,nrow(X))*apply(TOEFS,2,sd), 4), 
                perc.975=round(apply(TOEFS,2,mean) + qt(.975,nrow(X))*apply(TOEFS,2,sd), 4))
    # R Squared
    colnames(RSQRS) <- lvs.names[endo==1]
    t.rb <- matrix(NA, sum(endo), 2)
    for (j in 1:ncol(RSQRS))
        t.rb[j,] <- round(unlist(t.test(RSQRS[,j],mu=r2.orig[j])[c(1,3)]), 4)
    RB <- data.frame(Original=round(r2.orig,4), Mean.Boot=round(apply(RSQRS, 2, mean), 3), 
                Std.Err=round(apply(RSQRS,2,sd),3), t.statis=t.rb[,1], p.value=t.rb[,2],
                perc.025=round(apply(RSQRS,2,mean) - qt(.975,nrow(X))*apply(RSQRS,2,sd), 3), 
                perc.975=round(apply(RSQRS,2,mean) + qt(.975,nrow(X))*apply(RSQRS,2,sd), 3))
    # Bootstrap Results
    res.boot <- list(weights=WB, loadings=LB, paths=PB, rsq=RB, total.efs=TE)
    return(res.boot)
}

