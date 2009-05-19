`plspm` <-
function(x, inner.mat, sets, modes=NULL, scheme="centroid", 
                          scaled=TRUE, boot.val=FALSE, plsr=FALSE)
{
    # =========================== ARGUMENTS ==============================
    # x: a numeric matrix or data.frame containing the manifest variables
    # inner.mat: a square boolean matrix indicating the path relations
    # sets: a list of vectors of indices indicating the manifest variables for each block
    # modes: a character vector indicating the type of measurement for each latent variable 
    #"A" for reflective measurement, "B" for formative measurement
    # scheme: the inner weighting scheme to be used: "cenrtoid", "factor", or "path"
    # scaled: a logical value indicating whether scale data is donde
    # boot.val:a logical value indicating whether bootstrap validation is done  
    # plsr: a logical value for calculating path coeffs by pls-regression
    
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
        warning("Invalid 'boot.val' argument. No bootstrap validation is done.")
        boot.val <- FALSE
    }   
    if (!is.logical(plsr)) plsr<-FALSE
    # ========================== INPUTS SETTING ==========================
    IDM <- inner.mat
    dimnames(IDM) <- list(lvs.names, lvs.names)
    lvs <- nrow(IDM)
    lvs.names <- rownames(IDM)
    mvs <- sum(unlist(lapply(sets, length)))
    blocks <- unlist(lapply(sets, length))
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
    out.ws <- pls.weights(X, IDM, sets, modes, scheme, blocklist)
    if (is.null(out.ws)) stop("The pls algorithm is non convergent") 
    out.weights <- round(out.ws[[1]], 4)
    cor.XY <- cor(X, X%*%out.ws[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
         w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    if (scaled) {
        Y.lvs <- round(X %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs), 4)
        Z.lvs <- Y.lvs
    } else
    {   # expressing LVs in scale of MVs
        Y.lvs <- round(DM %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs), 4)
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
    loadcomu <- pls.loads(X, Y.lvs, blocklist, modes)    
    loads <- round(loadcomu[[1]], 4)
    comu <- round(loadcomu[[2]], 4)
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1  
    redun <- rep(0,mvs)
    for (k in 1:lvs)
        if (endo[k]==1)
            redun[which(blocklist==k)] <- comu[which(blocklist==k)] * R2[k]
    # ========================= Measurement model ========================
    outcor <- as.list(1:lvs)
    outmod <- as.list(1:lvs)
    for (k in 1:lvs)
    {
        aux <- which(blocklist==k)
        outmod[[k]] <- round(cbind(weights=out.weights[aux], std.loads=loads[aux], 
                             communal=comu[aux], redundan=redun[aux]), 4)
        outcor[[k]] <- round(cor(DM[,aux], Y.lvs), 4)
    }
    names(outmod) <- lvs.names
    names(outcor) <- lvs.names  
    # ======================== Unidimensionality =========================
    gof.om <- rep(1, lvs)# GoF relative outer model
    Alpha <- rep(1, lvs)# Cronbach's Alpha for each block
    Rho <- rep(1, lvs)# D.G. Rho for each block
    eig.1st <- rep(1,lvs)# first eigenvalue
    eig.2nd <- rep(0,lvs)# second eigenvalue
    for (aux in 1:lvs) 
    {      
        if (blocks[aux] != 1) 
        { 
            X.uni <- scale(DM[,which(blocklist==aux)])            
            if (nrow(X.uni)<ncol(X.uni)) {
                acp <- princomp(t(X.uni)) 
                gof.om[aux] <- round(mean(cor(t(X.uni),acp$scores[,1])^2),4)
                if (modes[aux]=="A") {
                    p = ncol(X.uni)
                    a.denom <- var(rowSums(X.uni))
                    a.numer <- 2*sum(cor(X.uni)[lower.tri(cor(X.uni))])
                    alpha <- round((a.numer / a.denom) * (p/(p-1)), 4)
                    Alpha[aux] <- alpha
                    numer.rho <- colSums(cor(t(X.uni), acp$scores[,1]))^2
                    denom.rho <- numer.rho + (p - colSums(cor(t(X.uni), acp$scores[,1])^2) )
                    Rho[aux] <- round(numer.rho / denom.rho, 4)
                }
            }
            if (nrow(X.uni)>=ncol(X.uni)) {
                acp <- princomp(X.uni)
                gof.om[aux]<-round(mean(cor(X.uni,acp$scores[,1])^2),4)
                if (modes[aux]=="A") {
                    p = ncol(X.uni)
                    a.denom <- var(rowSums(X.uni))
                    a.numer <- 2*sum(cor(X.uni)[lower.tri(cor(X.uni))])
                    alpha <- round((a.numer / a.denom) * (p/(p-1)), 4)
                    Alpha[aux] <- alpha
                    numer.rho <- colSums(cor(X.uni, acp$scores[,1]))^2
                    denom.rho <- numer.rho + (p - colSums(cor(X.uni, acp$scores[,1])^2) )
                    Rho[aux] <- round(numer.rho / denom.rho, 4)
                }
            }
            eig.1st[aux] <- round(acp$sdev[1]^2, 2)
            eig.2nd[aux] <- round(acp$sdev[2]^2, 2)
            if (modes[aux]=="B") {
                Alpha[aux] <- 0
                Rho[aux] <- 0
            }
        }
    }
    unidim <- data.frame(Type.measure=Mode, MVs=blocks, eig.1st,
                         eig.2nd, C.alpha=Alpha, DG.rho=Rho)
    rownames(unidim) <- lvs.names
    # ======================== Summary Inner model =======================
    exo.endo <- rowSums(IDM)
    exo.endo[rowSums(IDM)==0] <- "Exogen"
    exo.endo[rowSums(IDM)!=0] <- "Endogen"
    av.comu <- rep(0,lvs)   # avergae communality
    av.redu <- rep(0,lvs)   # average redundancy
    ave <- rep(0, lvs)      # average variance extracted
    for (k in 1:lvs)
    {
        av.comu[k] <- round(mean(comu[which(blocklist==k)]), 3)
        av.redu[k] <- round(mean(redun[which(blocklist==k)]), 3)
        if (modes[k]=="A")
        {
            ave.num <- sum(comu[which(blocklist==k)])
            ave.denom <- sum(comu[which(blocklist==k)]) + sum(1-(comu[which(blocklist==k)]))
            ave[k] <- round(ave.num / ave.denom, 3)
        }
    }
    names(ave) <- lvs.names
    innsum = data.frame(LV.Type=exo.endo, Measure=abbreviate(Mode,5), MVs=blocks, 
                        R.square=R2, Av.Commu=av.comu, Av.Redun=av.redu, AVE=ave)
    rownames(innsum) <- lvs.names
    # ============================ GoF Indexes ===========================
    gof.im <- rep(1, lvs)# GoF relative inner model
    for (k in 1:lvs)
    {
        if (endo[k]==1)
        {
             EB <- DM[,which(blocklist==k)]  # endog block
             aux <- which(IDM[k,]==1)
             SB <- DM[,which(blocklist %in% aux)]  # super block             
             gof.im[k] <- round(cancor(EB, SB)$cor[1]^2, 4)
        }
    }
    gof.abs <- round(sqrt(mean(av.comu)*(sum(R2)/sum(endo))), 4)
    gof.out <- round(mean(av.comu/gof.om), 4)    
    gof.inn <- round(sum(av.redu/gof.im)/sum(endo), 4)
    gof.rel <- round(sqrt(gof.out*gof.inn), 4)
    gof <- data.frame(GoF=c("Absolute", "Relative", "Outer.mod", "Inner.mod"),
                      value=c(gof.abs, gof.rel, gof.out, gof.inn))
    # Results
    model <- list(IDM, blocks, scheme, Mode, scaled, obs=nrow(DM), boot.val, sets, plsr)
    res <- list(unidim=unidim, outer.mod=outmod, inner.mod=innmod, latents=Y.lvs, scores=Z.lvs,
               out.weights=out.weights, loadings=loads, path.coefs=Path, r.sqr=R2, 
               outer.cor=outcor, inner.sum=innsum, effects=Path.efs, model=model, gof=gof)
    # ========================= Bootstrap Validation =========================
    if (boot.val) 
    {
        if (nrow(X) <= 10) warning("Bootstrapping stopped: very few cases.")
        if (nrow(X) > 10)
        {
            bootnum <- 200
            path.labs <- NULL
            efs.labs <- NULL
            for (j in 1:lvs)
                for (i in j:lvs)
                     if (IDM[i,j]==1) 
                         path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))
            endo <- rowSums(IDM)
            endo[endo!=0] <- 1
            path.orig <- as.vector(Path[which(IDM==1)])
            r2.orig <- R2[which(endo==1)]
            WEIGS <- matrix(NA, bootnum, mvs)
            LOADS <- matrix(NA, bootnum, mvs)
            PATHS <- matrix(NA, bootnum, sum(IDM))
            TOEFS <- matrix(NA, bootnum, nrow(Path.efs))
            RSQRS <- matrix(NA, bootnum, sum(endo))
            for (i in 1:bootnum)
            {
                boot.obs <- sample(1:nrow(X), nrow(X), replace=TRUE)
                X.boot <- X[boot.obs,]
                w.boot <- pls.weights(X.boot, IDM, sets, modes, scheme, blocklist)
                if (is.null(w.boot)) stop("Bootstrapping failed") 
                WEIGS[i,] <- w.boot[[1]]
                Y.boot <- X.boot %*% w.boot[[2]]
                pathmod <- pls.paths(IDM, Y.boot, plsr)
                P.boot <- pathmod[[2]]
                Toef.boot <- pls.efects(P.boot)
                PATHS[i,] <- as.vector(P.boot[which(IDM==1)])
                TOEFS[i,] <- Toef.boot[,4]
                RSQRS[i,] <- pathmod[[3]][which(endo==1)]
                l.boot <- pls.loads(X.boot, Y.boot, blocklist, modes)    
                LOADS[i,] <- l.boot[[1]]
            }
            # Outer weights
            colnames(WEIGS) <- mvs.names
            t.wb <- matrix(NA, mvs, 2)
            for (j in 1:ncol(WEIGS))
                t.wb[j,] <- round(unlist(t.test(WEIGS[,j],mu=out.weights[j])[c(1,3)]), 4)
            WB <- data.frame(Original=out.weights, Mean.Boot=round(apply(WEIGS,2,mean), 4), 
                        Std.Err=round(apply(WEIGS,2,sd),4), t.statis=t.wb[,1], p.value=t.wb[,2],
                        perc.025=round(apply(WEIGS,2,mean) - qt(.975,nrow(X))*apply(WEIGS,2,sd), 4), 
                        perc.975=round(apply(WEIGS,2,mean) + qt(.975,nrow(X))*apply(WEIGS,2,sd), 4))
            # Loadings
            colnames(LOADS) <- mvs.names
            t.lb <- matrix(NA, mvs, 2)
            for (j in 1:ncol(LOADS)) {
                if (loads[j]!=1)
                    t.lb[j,] <- round(unlist(t.test(LOADS[,j],mu=loads[j])[c(1,3)]), 4)
                else
                    t.lb[j,] <- c(0,1)
            }
            LB <- data.frame(Original=loads, Mean.Boot=round(apply(LOADS,2,mean), 4), 
                        Std.Err=round(apply(LOADS,2,sd),4), t.statis=t.lb[,1], p.value=t.lb[,2],
                        perc.025=round(apply(LOADS,2,mean) - qt(.975,nrow(X))*apply(LOADS,2,sd), 4), 
                        perc.975=round(apply(LOADS,2,mean) + qt(.975,nrow(X))*apply(LOADS,2,sd), 4))
    
            # Path Coefficients
            colnames(PATHS) <- path.labs
            t.pb <- matrix(NA, sum(IDM), 2)
            for (j in 1:ncol(PATHS))
                t.pb[j,] <- round(unlist(t.test(PATHS[,j],mu=path.orig[j])[c(1,3)]), 4)
            PB <- data.frame(Original=path.orig, Mean.Boot=round(apply(PATHS,2,mean), 4), 
                        Std.Err=round(apply(PATHS,2,sd),4), t.statis=t.pb[,1], p.value=t.pb[,2],
                        perc.025=round(apply(PATHS,2,mean) - qt(.975,nrow(X))*apply(PATHS,2,sd), 4), 
                        perc.975=round(apply(PATHS,2,mean) + qt(.975,nrow(X))*apply(PATHS,2,sd), 4))
            # Total Effects
            colnames(TOEFS) <- Path.efs[,1]
            t.teb <- matrix(NA, nrow(Path.efs), 2)
            for (j in 1:ncol(TOEFS))
                t.teb[j,] <- round(unlist(t.test(TOEFS[,j],mu=Path.efs[j,4])[c(1,3)]), 4)
            TE <- data.frame(Original=Path.efs[,4], Mean.Boot=round(apply(TOEFS,2,mean), 4), 
                        Std.Err=round(apply(TOEFS,2,sd), 4), t.statis=t.teb[,1], p.value=t.teb[,2],
                        perc.025=round(apply(TOEFS,2,mean) - qt(.975,nrow(X))*apply(TOEFS,2,sd), 4), 
                        perc.975=round(apply(TOEFS,2,mean) + qt(.975,nrow(X))*apply(TOEFS,2,sd), 4))
            # R Squared
            colnames(RSQRS) <- lvs.names[endo==1]
            t.rb <- matrix(NA, sum(endo), 2)
            for (j in 1:ncol(RSQRS))
                t.rb[j,] <- round(unlist(t.test(RSQRS[,j],mu=r2.orig[j])[c(1,3)]), 4)
            RB <- data.frame(Original=r2.orig, Mean.Boot=round(apply(RSQRS, 2, mean), 3), 
                        Std.Err=round(apply(RSQRS,2,sd),3), t.statis=t.rb[,1], p.value=t.rb[,2],
                        perc.025=round(apply(RSQRS,2,mean) - qt(.975,nrow(X))*apply(RSQRS,2,sd), 3), 
                        perc.975=round(apply(RSQRS,2,mean) + qt(.975,nrow(X))*apply(RSQRS,2,sd), 3))
            # Bootstrap Results
            res.boot <- list(weights=WB, loadings=LB, paths=PB, rsq=RB, total.efs=TE)
            res <- list(unidim=unidim, outer.mod=outmod, inner.mod=innmod, latents=Y.lvs, scores=Z.lvs,
                out.weights=out.weights, loadings=loads, path.coefs=Path, r.sqr=R2, outer.cor=outcor, 
                inner.sum=innsum, effects=Path.efs, gof=gof, boot=res.boot, model=model)
        }
    } # end 'if' bootstrapping
    # --------------------------------------------------------------------
    class(res) <- "plspm"
    return(res)
}

