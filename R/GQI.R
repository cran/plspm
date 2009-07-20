`GQI` <-
function(pls, part)
{
    # ========================== GQI function ==========================
    # Function to calculate Group Quality Index (GQI)        
    # =========================== arguments ==============================
    # pls: object of class "plspm"
    # part: vector with units memberships / or categorical variable
    
    IDM <- pls$model[[1]]# Inner Design Matrix
    blocks <- pls$model[[2]]# cardinality of blocks
    scheme <- pls$model[[3]]# inner weighting scheme
    modes <- pls$model[[4]]# measurement modes
    scaled <- pls$model[[5]]# type of scaling
    plsr <- FALSE 
    DM <- pls$data
    lvs <- nrow(IDM)
    lvs.names <- rownames(IDM)
    mvs <- sum(blocks)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1  
    # data scaling (standardized data)
    X <- scale(DM)

    clas.part <- part 
    # number of clusters
    nclus <- nlevels(factor(clas.part))
    w.locals <- as.list(1:nclus)# outer.weights
    LV.locals <- as.list(1:nclus)# std latent variables 
    loads.locals <- as.list(1:nclus)# loadings
    path.locals <- as.list(1:nclus)# path coefficients
    R2.locals <- as.list(1:nclus)# R2
    comu.locals <- as.list(1:nclus)# mvs communalities
    outres.locals <- as.list(1:nclus)# communality residuals
    innres.locals <- as.list(1:nclus)# structural residuals
    out.term <- as.list(1:nclus)# outer term for GQI
    inn.term <- as.list(1:nclus)# inner term for GQI
    gqi.locals <- rep(0,nclus)# pseudo-gqi for each class

    # define MV matrix for each initial class
    split.DM <- as.list(1:nclus)
    split.X <- as.list(1:nclus)
    for (k in 1:nclus)
        split.DM[[k]] <- DM[clas.part==k,]            

    # local models computation
    for (k in 1:nclus)
    {   
        mean.k <- apply(split.DM[[k]],2,mean)# local mean
        sd.k <- apply(split.DM[[k]],2,sd)# local std.dev
        # spliting data matrix for each class
        split.X[[k]] <- scale(split.DM[[k]], center=mean.k, scale=sd.k)
        # calculating outer weights for each class
        out.ws  <- pls.weights(split.X[[k]], IDM, blocks, modes, scheme)
        w.locals[[k]] <- out.ws[[2]]
        # calculating LV scores for each class
        LV.locals[[k]] <- split.X[[k]] %*% out.ws[[2]]
        # calculating path coefficients for each class
        pathmod <- pls.paths(IDM, LV.locals[[k]], plsr)
        path.locals[[k]] <- pathmod[[2]]
        R2.locals[[k]] <- pathmod[[3]][endo==1]
        # calculating loadings and communalities for each class
        loadcomu <- pls.loads(split.X[[k]], LV.locals[[k]], blocks)    
        loads.locals[[k]] <- loadcomu[[1]]
        comu.locals[[k]] <- loadcomu[[2]]
        # computation of communality residuals (squared)
        out.res <- split.X[[k]]
        for (j in 1:lvs)
        {
            q <- which(blocklist==j) 
            X.hat <- LV.locals[[k]][,j] %*% t(loads.locals[[k]][q])
            out.res[,q] <- (split.X[[k]][,q] - X.hat)^2# outer residuals
        }
        outres.locals[[k]] <- out.res
        # computation of inner residuals (squared)
        if (sum(endo)!=1)
            Y.hat <- LV.locals[[k]] %*% t(path.locals[[k]][endo==1,])   
        if (sum(endo)==1)
            Y.hat <- LV.locals[[k]] %*% path.locals[[k]][endo==1,]        
        innres.locals[[k]] <- (LV.locals[[k]][,endo==1] - Y.hat)^2
        # outer and inner terms of GQI formula
        out.term[[k]] <- mean((1-(colSums(outres.locals[[k]])/colSums(split.X[[k]]^2))))
        if (sum(endo)==1)
            inn.term[[k]] <- mean((1-(colSums(innres.locals[[k]])/sum(LV.locals[[k]][,endo==1]^2))))
        if (sum(endo)!=1)
            inn.term[[k]] <- mean((1-(colSums(innres.locals[[k]])/colSums(LV.locals[[k]][,endo==1]^2))))
        gqi.locals[k] <- out.term[[k]] * inn.term[[k]]
    }
    unit.prop <- unlist(lapply(split.X, nrow))/nrow(DM)# proportion of units in each class
    GQI <- sqrt(sum(gqi.locals * unit.prop))
    return(GQI)
}

