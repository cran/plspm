`plspm` <-
function(x, inner.mat, sets, modes=NULL, scheme="centroid", 
                  scaled=TRUE, boot.val=FALSE, br=NULL, plsr=FALSE)
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
    # boot.val:a logical value indicating whether bootstrap validation is done 
    # br: an integer indicating the number of bootstraps resamples, used 
    #     only when boot.val=TRUE, (50 <= br <= 500)
    # plsr: a logical value for calculating path coeffs by pls-regression
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
    SCHEMES <- c("centroid", "factor", "path")
    scheme <- pmatch(scheme, SCHEMES)
    if (is.na(scheme)) {
        warning("Invalid argument 'scheme'. Default 'scheme=centroid' is used.")   
        scheme <- "centroid"
    }
    if (!is.logical(scaled)) {
        warning("Invalid argument 'scaled'. Default 'scaled=TRUE' is used.")
        scaled <- TRUE
    }
    if (!is.logical(boot.val)) {
        warning("Invalid argument 'boot.val'. No bootstrap validation is done.")
        boot.val <- FALSE
    }   
    if (boot.val) {
        if (!is.null(br)) {        
            if (mode(br)!="numeric" || length(br)!=1 || (br%%1)!=0 ||
                br<50 || br>500) {
                warning("Invalid argument 'br'. Default 'br=100' is used.")   
                br <- 100
            } 
        } else
            br <- 100
    }
    if (!is.logical(plsr)) plsr<-FALSE

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
    if(scaled) X=scale(DM) else X=scale(DM, scale=FALSE)
    # ==================== Stage 1: Iterative procedure ==================
    out.ws <- pls.weights(X, IDM, blocks, modes, scheme)
    if (is.null(out.ws)) stop("The pls algorithm is non convergent") 
    out.weights <- round(out.ws[[1]], 4)
    cor.XY <- cor(X, X%*%out.ws[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
         w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    if (scaled) {
        Y.lvs <- X %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs)
        Z.lvs <- Y.lvs
    } else
    {   # expressing LVs in scale of MVs
        Y.lvs <- DM %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs)
        mv.ran <- range(DM)
        lv.ran <- range(Y.lvs)
        Z.lvs <- (mv.ran[2]-mv.ran[1])/(lv.ran[2]-lv.ran[1]) * Y.lvs
        if (range(Z.lvs)[2] > mv.ran[2]) 
            Z.lvs <- Z.lvs - (range(Z.lvs)[2]-mv.ran[2])
        if (range(Z.lvs)[2] < mv.ran[2])        
            Z.lvs <- Z.lvs + (mv.ran[1] - range(Z.lvs)[1]) 
        Z.lvs <- round(Z.lvs, 4)
    }
    dimnames(Y.lvs) <- list(rownames(X), lvs.names)
    dimnames(Z.lvs) <- list(rownames(X), lvs.names)
    # ============ Stage 2: Path coefficients and total effects ==========
    pathmod <- pls.paths(IDM, Y.lvs, plsr)
    innmod <- pathmod[[1]]
    Path <- pathmod[[2]]
    R2 <- pathmod[[3]]
    Path.efs <- pls.efects(Path)
    # ========== Stage 3: Measurement loadings and communalities =========
    loadcomu <- pls.loads(X, Y.lvs, blocks)    
    loads <- loadcomu[[1]]
    comu <- loadcomu[[2]]
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1  
    redun <- rep(0,mvs)
    for (j in 1:lvs)
        if (endo[j]==1)
            redun[blocklist==j] <- comu[blocklist==j] * R2[j]
    # ========================= Measurement model ========================
    outcor <- as.list(1:lvs)
    outmod <- as.list(1:lvs)
    for (j in 1:lvs)
    {
        aux <- which(blocklist==j)
        outmod[[j]] <- round(cbind(weights=out.weights[aux], std.loads=loads[aux], 
                             communal=comu[aux], redundan=redun[aux]), 4)
        outcor[[j]] <- round(cor(DM[,aux], Y.lvs), 4)
    }
    names(outmod) <- lvs.names
    names(outcor) <- lvs.names  
    # ======================== Unidimensionality =========================
    uni.gof <- pls.unidim(DM, blocks, modes)
    unidim <- uni.gof[[1]]# unidimensionality
    gof.om <- uni.gof[[2]]# GoF relative outer model
    # ======================== Summary Inner model =======================
    exo.endo <- rowSums(IDM)
    exo.endo[rowSums(IDM)==0] <- "Exogen"
    exo.endo[rowSums(IDM)!=0] <- "Endogen"
    av.comu <- rep(0,lvs)   # avergae communality
    av.redu <- rep(0,lvs)   # average redundancy
    ave <- rep(0, lvs)      # average variance extracted
    for (k in 1:lvs)
    {
        av.comu[k] <- mean(comu[which(blocklist==k)])
        av.redu[k] <- mean(redun[which(blocklist==k)])
        if (modes[k]=="A")
        {
            ave.num <- sum(comu[which(blocklist==k)])
            ave.denom <- sum(comu[which(blocklist==k)]) + sum(1-(comu[which(blocklist==k)]))
            ave[k] <- round(ave.num / ave.denom, 3)
        }
    }
    names(ave) <- lvs.names
    innsum = data.frame(LV.Type=exo.endo, Measure=abbreviate(Mode,5), MVs=blocks, 
           R.square=round(R2,4), Av.Commu=round(av.comu,4), Av.Redun=round(av.redu,4), AVE=ave)
    rownames(innsum) <- lvs.names
    # ============================ GoF Indexes ===========================
    gof.im <- rep(1, lvs)# GoF relative inner model
    for (j in 1:lvs)
    {
        if (endo[j]==1)
        {
             EB <- DM[,blocklist==j]  # endog block
             aux <- which(IDM[j,]==1)
             SB <- DM[,which(blocklist %in% aux)]  # super block             
             gof.im[j] <- round(cancor(EB, SB)$cor[1]^2, 4)
        }
    }
    gof.abs <- pls.gof(comu, R2, blocks, IDM)
    gof.out <- round(mean(av.comu/gof.om), 4)    
    gof.inn <- round(sum(av.redu/gof.im)/sum(endo), 4)
    gof.rel <- round(sqrt(gof.out*gof.inn), 4)
    gof <- data.frame(GoF=c("Absolute", "Relative", "Outer.mod", "Inner.mod"),
                      value=c(round(gof.abs,4), gof.rel, gof.out, gof.inn))
    # Results
    skem <- switch(scheme, "centroid"="centroid", "factor"="factor", "path"="path")
    model <- list(IDM=IDM, blocks=blocks, scheme=skem, modes=modes, scaled=scaled, 
                  boot.val=boot.val, plsr=plsr, obs=nrow(X), br=br)
    res <- list(unidim=unidim, outer.mod=outmod, inner.mod=innmod, latents=Y.lvs, scores=Z.lvs,
               out.weights=out.weights, loadings=loads, path.coefs=Path, r.sqr=R2, 
               outer.cor=outcor, inner.sum=innsum, effects=Path.efs, gof=gof, data=DM, model=model)
    # ========================= Bootstrap Validation =========================
    if (boot.val) 
    {
        if (nrow(X) <= 10) {
            warning("Bootstrapping stopped: very few cases.") 
        } else 
        { 
            n.efs <- nrow(Path.efs)
            res.boot <- pls.boot(DM, IDM, blocks, modes, scheme, scaled, br, plsr)
        }
        res <- list(unidim=unidim, outer.mod=outmod, inner.mod=innmod, latents=Y.lvs, scores=Z.lvs,
            out.weights=out.weights, loadings=loads, path.coefs=Path, r.sqr=R2, outer.cor=outcor, 
            inner.sum=innsum, effects=Path.efs, gof=gof, boot=res.boot, data=DM, model=model)
    } # end 'if' bootstrapping
    # --------------------------------------------------------------------
    class(res) <- "plspm"
    return(res)
}

