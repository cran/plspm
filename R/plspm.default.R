`plspm.default` <-
function(x, inner.mat, sets, modes=NULL, scheme="factor", 
                          scaled=TRUE, boot.val=FALSE, plsr=FALSE)
{
    # =========================== ARGUMENTS ==============================
    # x: a numeric matrix or data.frame containing the manifest variables
    # inner.mat: a square boolean matrix indicating the path relations
    # sets: a list of vectors of indices indicating the manifest variables for each block
    # modes: a character vector indicating the type of measurement for each latent variable 
    #"A" for reflective measurement, "B" for formative measurement
    # scheme: the inner weighting scheme to be used: "factor", "centroid", or "path"
    # scaled: a logical value indicating whether scale data is donde
    # boot.val:a logical value indicating whther bootstrap validation is done  
    
    # ======================== Internal functions ========================
    plsr1 <- function(x, y, nc=2, scaled=TRUE)
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
    #------------------------------
    pls.weights <- function(X, IDM, sets, scheme, blocklist)
    {
        lvs <- nrow(IDM)
        mvs <- ncol(X)
        lvs.names <- rownames(IDM)
        mvs.names <- colnames(X)
        # initialize arbitrary outer weights (value of '1' for all weights)
        out.ws <- sets
        for (k in 1:length(sets))
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
            # computing outer weights 'w'
            for (k in 1:lvs)
            {
                X.blok = X[,which(blocklist==k)] 
                if (modes[k]=="A")# reflective way
                    ODM[which(blocklist==k),k] = solve(t(Z[,k])%*%Z[,k])%*%Z[,k] %*% X.blok
                if (modes[k]=="B")# formative way
                    ODM[which(blocklist==k),k] = solve.qr(qr(X.blok),Z[,k])
            }
            W <- ODM %*% diag(1/sd(X %*% ODM),lvs,lvs)
            w.new = rowSums(W)                
            w.dif <- sum((w.old - w.new)^2)  # difference of out.weights 
            w.old <- w.new
            if (sum(w.dif^2)<1e-07 || itermax==200) break
            itermax <- itermax + 1
        } # end repeat       
        res.ws <- list(w.new, W)
        if (itermax==200) res.ws=NULL
        return(res.ws)
    } # end function pls.weights
    #------------------------------
    pls.paths <- function(IDM, Y.lvs, plsr)
    {
        lvs.names <- colnames(IDM)
        endo = rowSums(IDM)
        endo[endo!=0] <- 1  # vector indicating endogenous LVs
        innmod <- as.list(1:sum(endo))
        Path <- IDM
        residuals <- as.list(1:sum(endo))
        R2 <- rep(0,nrow(IDM))
        aux1 <- 1
        for (aux in 1:nrow(IDM)) 
        {
            if (endo[aux] == 1) # endogenous LV
            {
                c <- which(IDM[aux,1:aux]==1)   
                if (length(c)>1 & plsr) {               
                    path.lm <- plsr1(Y.lvs[,c], Y.lvs[,aux])
                    Path[aux,c] <- path.lm$coeffs
                    residuals[[aux1]] <- path.lm$resid
                    R2[aux] <- path.lm$R2[1]
                    inn.val <- round(c(path.lm$R2[1], path.lm$cte, path.lm$coeffs), 3)
                    inn.lab <- c("R2", "Intercept", paste(rep("path_",length(c)),names(c),sep=""))
                    innmod[[aux1]] <- data.frame(concept=inn.lab, value=inn.val)
                    aux1 <- aux1 + 1   
                }
                if (length(c)==1 | !plsr) {
                    path.lm <- summary(lm(Y.lvs[,aux] ~ Y.lvs[,c]))
                    Path[aux,c] <- round(path.lm$coef[-1,1], 4)
                    residuals[[aux1]] <- path.lm$residuals  
                    R2[aux] <- round(path.lm$r.squared, 3)
                    inn.val <- round(c(path.lm$r.squared, path.lm$coef[,1]), 3)
                    inn.lab <- c("R2", "Intercept", paste(rep("path_",length(c)),names(c),sep=""))
                    innmod[[aux1]] <- data.frame(concept=inn.lab, value=inn.val)
                    aux1 <- aux1 + 1   
                }
            }
        }
        names(innmod) <- lvs.names[endo!=0]  
        names(R2) <- lvs.names
        res.paths <- list(innmod, Path, R2)
        return(res.paths)
    } # end function 'pls.paths'
    #------------------------------
    pls.loads <- function(X, Y.lvs, blocklist, modes)
    {
        lvs <- length(modes)
        mvs <- ncol(X)
        loads <- rep(NA, mvs)
        comu <- rep(NA, mvs)
        for (k in 1:lvs)
        {
            X.blok <- X[,which(blocklist==k)] 
            loads[which(blocklist==k)] <- cor(X.blok, Y.lvs[,k])
            comu[which(blocklist==k)] <- cor(X.blok, Y.lvs[,k])^2
        }
        names(loads) <- colnames(X)  
        names(comu) <- colnames(X)
        res.loads <- list(loads, comu)
        return(res.loads)
    } # end function 'pls.loads'
    #------------------------------
    pls.efects <- function(Path)
    {
        lvs <- nrow(IDM)
        lvs.names <- rownames(IDM)
        path.efects <- as.list(1:(lvs-1))
        path.efects[[1]] <- Path
        for (k in 2:(lvs-1))
            path.efects[[k]] <- round(path.efects[[k-1]] %*% Path, 4)
        ind.paths <- matrix(0, lvs, lvs)
        for (k in 2:length(path.efects))
            ind.paths <- ind.paths + path.efects[[k]]
        total.paths <- Path + ind.paths
        efs.labs <- NULL
        dir.efs <- NULL
        ind.efs <- NULL
        tot.efs <- NULL
        for (j in 1:lvs)
            for (i in j:lvs)
                 if (total.paths[i,j]!=0) 
                 {
                     efs.labs <- c(efs.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))
                     dir.efs <- c(dir.efs, Path[i,j])# direct effects
                     ind.efs <- c(ind.efs, ind.paths[i,j])# indirect effects
                     tot.efs <- c(tot.efs, total.paths[i,j])# total effects
                 }
        Effects <- data.frame(relationships=efs.labs, dir.effects=dir.efs, 
                              ind.effects=ind.efs, tot.effects=tot.efs)
        return(Effects)
    }        
    # -------------------------------------------------------------------- 

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
    if (!is.na(pmatch(scheme, "factor"))) 
        scheme <- "factor"
    SCHEMES <- c("factor", "centroid", "path")
    scheme <- pmatch(scheme, SCHEMES)
    if (is.na(scheme)) {
        warning("Invalid argument 'scheme'. Default 'scheme=factor' is used.")   
        scheme <- "factor"
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
    lvs <- nrow(IDM)
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
        mvs.names[which(blocklist==k)] <- colnames(x[,sets[[k]]])
    }
    dimnames(DM) <- list(rownames(x), mvs.names)
    # apply the selected scaling
    if(scaled) X=scale(DM) else X=scale(DM, scale=FALSE)
    # ==================== Stage 1: Iterative procedure ==================
    out.ws <- pls.weights(X, IDM, sets, scheme, blocklist)
    if (is.null(out.ws)) stop("The algorithm is non convergent") 
    out.weights <- round(out.ws[[1]], 4)
    if (scaled) {
        Y.lvs <- round(X %*% out.ws[[2]], 4)
        Z.lvs <- Y.lvs
    }
    if (!scaled) {   # expressing LVs in scale of MVs
        Y.lvs <- round(DM %*% out.ws[[2]], 4)
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
    model <- list(IDM, blocks, scheme, Mode, scaled, obs=nrow(DM), boot.val)
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
                w.boot <- pls.weights(X.boot, IDM, sets, scheme, blocklist)
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
            for (j in 1:ncol(LOADS))
                t.lb[j,] <- round(unlist(t.test(LOADS[,j],mu=loads[j])[c(1,3)]), 4)
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

