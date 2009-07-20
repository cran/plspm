`locals.test` <-
function(X, pls, g)
{
    # ======================== local.test function ======================
    # Internal function related with "rebus.test" to calculate the
    # the permutation test (for paths, loadings, and gof) among 2 groups
    # =========================== arguments ===============================
    # X: data matrix related with g
    # pls: an object of class "plspm"
    # g: a factor with 2 levels indicating the groups to be compared

    # ========================== INPUTS SETTING ==========================
    IDM <- pls$model[[1]]# Inner Design Matrix
    blocks <- pls$model[[2]]# cardinality of blocks
    scheme <- pls$model[[3]]# inner weighting scheme
    modes <- pls$model[[4]]# measurement modes
    scaled <- pls$model[[5]]# type of scaling
    plsr <- pls$model[[7]]# pls-regression
    lvs <- nrow(IDM)
    lvs.names <- rownames(IDM)
    mvs <- sum(blocks)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    reps <- 100
    path.labs <- NULL
    efs.labs <- NULL
    for (j in 1:lvs)
        for (i in j:lvs)
             if (IDM[i,j]==1) 
                 path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))

    # ====================== Group1 model estimation =====================
    g1.lab <- levels(g)[1]
    group1 <- which(g==levels(g)[1])
    if(scaled) X.g1=scale(X[group1,]) else X.g1=scale(X[group1,], scale=FALSE)
    wgs.g1 <- pls.weights(X.g1, IDM, blocks, modes, scheme)
    if (is.null(wgs.g1)) stop("The algorithm is non convergent") 
    cor.XY <- cor(X.g1, X.g1%*%wgs.g1[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
         w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    Y1.lvs <- round(X.g1 %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs), 4)
    dimnames(Y1.lvs) <- list(rownames(X.g1), lvs.names)
    # Path coefficients 
    pathmod.g1 <- pls.paths(IDM, Y1.lvs, plsr)
    innmod.g1 <- pathmod.g1[[1]]
    Path.g1 <- pathmod.g1[[2]]
    R2.g1 <- pathmod.g1[[3]]    
    path.g1 <- as.vector(Path.g1[which(IDM==1)])
    names(path.g1) <- path.labs
    # calculating loadings and communalities for each class
    loadcomu <- pls.loads(X.g1, Y1.lvs, blocks)   
    load.g1 <- loadcomu[[1]]
    # gof
    gof.g1 <- pls.gof(load.g1^2, R2.g1, blocks, IDM)

    # ====================== Group2 model estimation =====================
    g2.lab <- levels(g)[2]
    group2 <- which(g==levels(g)[2])
    if(scaled) X.g2=scale(X[group2,]) else X.g2=scale(X[group2,], scale=FALSE)
    wgs.g2 <- pls.weights(X.g2, IDM, blocks, modes, scheme)
    if (is.null(wgs.g2)) stop("The algorithm is non convergent") 
    cor.XY <- cor(X[group2,], X.g2%*%wgs.g2[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
         w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    Y2.lvs <- round(X.g2 %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs), 4)
    dimnames(Y2.lvs) <- list(rownames(X.g2), lvs.names)
    # Path coefficients 
    pathmod.g2 <- pls.paths(IDM, Y2.lvs, plsr)
    innmod.g2 <- pathmod.g2[[1]]
    Path.g2 <- pathmod.g2[[2]]
    R2.g2 <- pathmod.g2[[3]]    
    path.g2 <- as.vector(Path.g2[which(IDM==1)])
    names(path.g2) <- path.labs
    # calculating loadings and communalities for each class
    loadcomu <- pls.loads(X.g2, Y2.lvs, blocks)   
    load.g2 <- loadcomu[[1]]
    # gof
    gof.g2 <- pls.gof(load.g2^2, R2.g2, blocks, IDM)

    # ====================== Group Comparison =====================
    difpath.orig <- abs(path.g1-path.g2)
    difload.orig <- abs(load.g1-load.g2)
    difgof.orig <- abs(gof.g1 - gof.g2)
    group1 <- which(g==levels(g)[1])
    group2 <- which(g==levels(g)[2])
    ng1 <- length(group1)
    ng2 <- length(group2)
    difpath.perm <- matrix(0, reps, sum(IDM))
    difload.perm <- matrix(0, reps, mvs)
    difgof.perm <- rep(0, reps)
    for (i in 1:reps)# multigroup permutation
    {
        permu <- sample(1:(ng1+ng2), ng1+ng2)
        samg1 <- permu[1:ng1]
        samg2 <- permu[(ng1+1):(ng1+ng2)]
        if(scaled) X.g1=scale(X[samg1,]) else X.g1=scale(X[samg1,], scale=FALSE)
        if(scaled) X.g2=scale(X[samg2,]) else X.g2=scale(X[samg2,], scale=FALSE)
        wgs.g1 <- pls.weights(X.g1, IDM, blocks, modes, scheme)
        wgs.g2 <- pls.weights(X.g2, IDM, blocks, modes, scheme)
        if (is.null(wgs.g1)) stop("Non convergence in bootstrap samples") 
        if (is.null(wgs.g2)) stop("Non convergence in bootstrap samples") 
        cor.XY <- cor(X.g1, X.g1%*%wgs.g1[[2]])
        w.sig <- rep(NA,lvs)
        for (j in 1:lvs) 
             w.sig[j] <- ifelse(sum(sign(cor.XY[blocklist==j,j]))<=0,-1,1)
        Y1.lvs <- round(X.g1 %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs), 4)
        cor.XY <- cor(X.g2, X.g2%*%wgs.g2[[2]])
        w.sig <- rep(NA,lvs)
        for (j in 1:lvs) 
             w.sig[j] <- ifelse(sum(sign(cor.XY[blocklist==j,j]))<=0,-1,1)
        Y2.lvs <- round(X.g2 %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs), 4)
        pathmod.g1 <- pls.paths(IDM, Y1.lvs, plsr)
        paths.g1 <- pathmod.g1[[2]]    
        pathmod.g2 <- pls.paths(IDM, Y2.lvs, plsr)
        paths.g2 <- pathmod.g2[[2]]
        loadcomu <- pls.loads(X.g1, Y1.lvs, blocks)   
        loads.g1 <- loadcomu[[1]]
        loadcomu <- pls.loads(X.g2, Y2.lvs, blocks)   
        loads.g2 <- loadcomu[[1]]  
        gofs.g1 <- pls.gof(loads.g1^2, R2.g1, blocks, IDM)  
        gofs.g2 <- pls.gof(loads.g2^2, R2.g2, blocks, IDM)
        # difference between groups
        pp1 <- as.vector(paths.g1[which(IDM==1)])
        pp2 <- as.vector(paths.g2[which(IDM==1)])
        difpath.perm[i,] <- abs(pp1 - pp2)
        difload.perm[i,] <- abs(loads.g1 - loads.g2)
        difgof.perm[i] <- abs(gofs.g1 - gofs.g2)
    }   
    # p-value for path coefficients
    path.perm <- difpath.orig 
    for (j in 1:sum(IDM))         
        path.perm[j] <- length(which(difpath.orig[j]<difpath.perm[,j])) + 1
    path.val <- (1/(reps+1))*path.perm 
    signi.path <- rep("no",length(path.val))
    signi.path[path.val<0.05] <- "yes"
    res.path <- round(cbind(path.g1, path.g2, difpath.orig, path.val), 4)
    res1 <- data.frame(res.path, signi.path)
    colnames(res1) <- c(paste(rep("Class",2),levels(g),sep="."), 
                         "diff.abs", "p.value", "sig.05")  
    # p-values for loadings
    load.perm <- difload.orig 
    for (j in 1:mvs)
        load.perm[j] <- length(which(difload.orig[j]<difload.perm[,j])) + 1
    load.val <- (1/(reps+1))*load.perm 
    signi.load <- rep("no",length(load.val))
    signi.load[load.val<0.05] <- "yes"
    res.load <- round(cbind(load.g1, load.g2, difload.orig, load.val), 4)
    res2 <- data.frame(res.load, signi.load)
    colnames(res2) <- c(paste(rep("Class",2),levels(g),sep="."), 
                         "diff.abs", "p.value", "sig.05")  
    # p-values for gof
    gof.perm <- length(which(difgof.orig<difgof.perm)) + 1
    gof.val <- (1/(reps+1))*gof.perm 
    signi.gof <- rep("no",length(gof.val))
    signi.gof[gof.val<0.05] <- "yes"
    res3 <- data.frame(round(gof.g1,4), round(gof.g2,4), 
               round(difgof.orig,4), round(gof.val,4), signi.gof)
    names(res3) <- c(paste(rep("Class",2),levels(g),sep="."), 
                         "diff.abs", "p.value", "sig.05")  
    # list with results 
    resul <- list(paths=res1, loadings=res2, gof=res3)
    return(resul)
}

