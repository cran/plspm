plspm.fit <-
function(x, inner.mat, sets, modes=NULL, scheme="centroid", scaled=TRUE)
{
    # =========================== ARGUMENTS ====================================
    # x: a numeric matrix or data.frame containing the manifest variables
    # inner.mat: a square boolean matrix indicating the path relations
    # sets: a list of vectors of indices indicating the manifest variables 
    #       for each block
    # modes: a character vector indicating the measurement type for each 
    #        latent variable. "A" reflective, "B" formative
    # scheme: a character string indicating the inner weighting scheme 
    #         to be used: "factor", "centroid", or "path"
    # scaled: a logical value indicating whether scale data is performed
    # ===========================================================================

    # ==================== Checking function arguments ===================
    if (!is.matrix(x) && !is.data.frame(x))
        stop("Invalid object 'x'. Must be a numeric matrix or data frame.")
    if (is.null(rownames(x))) 
        rownames(x) <- 1:nrow(x)
    if (is.null(colnames(x))) 
        colnames(x) <- paste("MV", 1:ncol(x), sep="")
    if (!is.matrix(inner.mat))
        stop("Invalid argument 'inner.mat'. Must be a matrix.")
    if (nrow(inner.mat)!=ncol(inner.mat))
        stop("Invalid argument 'inner.mat'. Must be a square matrix.")
    for (j in 1:ncol(inner.mat))
        for (i in 1:nrow(inner.mat)) {
            if (i<=j)
                if (inner.mat[i,j]!=0) 
                    stop("argument 'inner.mat' must be a lower triangular matrix")
            if (length(intersect(inner.mat[i,j], c(1,0)))==0)
                stop("elements in 'inner.mat' must be '1' or '0'")
        }
    if (is.null(dimnames(inner.mat)))
        lvs.names <- paste("LV", 1:ncol(inner.mat), sep="")
    if (!is.null(rownames(inner.mat)))
        lvs.names <- rownames(inner.mat)
    if (!is.null(colnames(inner.mat)))
        lvs.names <- colnames(inner.mat)
    if (!is.list(sets))
        stop("Invalid argument 'sets'. Must be a list.")
    if (length(sets) != nrow(inner.mat))
        stop("Number of rows of 'inner.mat' does not coincide with length of 'sets'.")
    if (is.null(modes)) {
        modes <- rep("A",length(sets))
        warning("Argument 'modes' missing. Default reflective 'modes' is used.")
    }
    if (length(sets) != length(modes)) {
        warning("Invalid length of 'modes'. Default reflective 'modes' is used.")
        modes <- rep("A", length(sets))
    }
    for (i in 1:length(modes))
        if (modes[i]!="A" && modes[i]!="B") modes[i]<-"A"
    if (!is.na(pmatch(scheme, "centroid"))) 
        scheme <- "centroid"
    SCHEMES <- c("centroid", "factor")
    scheme <- pmatch(scheme, SCHEMES)
    if (is.na(scheme)) {
        warning("Invalid argument 'scheme'. Default 'scheme=centroid' is used.")   
        scheme <- "centroid"
    }
    if (!is.logical(scaled)) {
        warning("Invalid argument 'scaled'. Default 'scaled=TRUE' is used.")
        scaled <- TRUE
    }
    plsr <- FALSE

    # ========================== INPUTS SETTING ==========================
    IDM <- inner.mat
    dimnames(IDM) <- list(lvs.names, lvs.names)
    lvs <- nrow(IDM)
    lvs.names <- rownames(IDM)
    mvs <- sum(unlist(lapply(sets, length)))
    blocks <- unlist(lapply(sets, length))
    names(blocks) <- lvs.names
    blocklist <- sets
    for (k in 1:length(sets))
         blocklist[[k]] <- rep(k,blocks[k])
    blocklist <- unlist(blocklist)
    Mode <- modes
    Mode[modes=="A"] <- "Reflective"
    Mode[modes=="B"] <- "Formative"   
    # building data matrix 'DM'
    DM <- matrix(NA, nrow(x), mvs)
    mvs.names <- rep(NA, mvs)
    for (k in 1:lvs)
    {        
        DM[,which(blocklist==k)] <- as.matrix(x[,sets[[k]]])
        mvs.names[which(blocklist==k)] <- colnames(x)[sets[[k]]]
    }
    dimnames(DM) <- list(rownames(x), mvs.names)
    # apply the selected scaling
    one.vec <- rep(1,nrow(DM))
    center <- diag(1,nrow(DM),nrow(DM)) - one.vec%*%t(one.vec)/nrow(DM)
    X <- center %*% DM
    if (scaled)   # standard data (var=1)
    {
        stdev.X <- sd(DM) * sqrt((nrow(DM)-1)/nrow(DM)) 
        X <- X %*% diag(1/stdev.X, ncol(DM), ncol(DM))
    }
    dimnames(X) <- list(rownames(x), mvs.names)

    # ==================== Stage 1: Iterative procedure ==================
    out.ws <- .pls.weights(X, IDM, blocks, modes, scheme)
    if (is.null(out.ws)) stop("The pls algorithm is non convergent") 
    out.weights <- round(out.ws[[1]], 4)
    cor.XY <- cor(X, X%*%out.ws[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
         w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    Z.lvs <- X %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs)
    Y.lvs <- Z.lvs
    if (!scaled) 
        Y.lvs <- DM %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs)
    dimnames(Y.lvs) <- list(rownames(X), lvs.names)
    dimnames(Z.lvs) <- list(rownames(X), lvs.names)
    # ============ Stage 2: Path coefficients and total effects ==========
    pathmod <- .pls.paths(IDM, Y.lvs, plsr)
    innmod <- pathmod[[1]]
    Path <- round(pathmod[[2]], 4)
    R2 <- round(pathmod[[3]], 4)
    # ========== Stage 3: Measurement loadings and communalities =========
    loadcomu <- .pls.loads(X, Y.lvs, blocks)    
    loads <- round(loadcomu[[1]], 4)
    comu <- loadcomu[[2]]
    # ========================= Measurement model ========================
    outmod <- as.list(1:lvs)
    for (j in 1:lvs)
    {
        aux <- which(blocklist==j)
        outmod[[j]] <- round(cbind(weights=out.weights[aux], std.loads=loads[aux], 
                             communal=comu[aux]), 4)
    }
    names(outmod) <- lvs.names
    # =========================== Basic Results ==========================
    skem <- switch(scheme, "centroid"="centroid", "factor"="factor", "path"="path")
    model <- list(IDM=IDM, blocks=blocks, scheme=skem, modes=modes, scaled=scaled, obs=nrow(X))
    res <- list(outer.mod=outmod, inner.mod=innmod, latents=Z.lvs, scores=Y.lvs,
               out.weights=out.weights, loadings=loads, path.coefs=Path, r.sqr=R2, 
               model=model)
    class(res) <- c("plspm.fit", "plspm")
    return(res)
}

