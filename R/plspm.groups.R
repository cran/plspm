`plspm.groups` <-
function(x, pls, g, method="bootstrap", reps=NULL)
{
    # =========================== ARGUMENTS ==============================
    # x: a numeric matrix or data.frame containing the manifest variables
    # pls: an object of class "plspm"
    # g: a factor with 2 levels indicating the groups to be compared
    # method: the method to be used: "bootstrap", or "permutation"
    # reps: number of bootstrap resamples or number of permutations

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
            Z <- scale(Z)
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
            if (sum(w.dif^2)<1e-07 || itermax==300) break
            itermax <- itermax + 1
        } # end repeat       
        res.ws <- list(w.new, W)
        if (itermax==300) res.ws=NULL
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
    # -------------------------------------------------------------------- 

    # ==================== Checking function arguments ===================
    if (!is.matrix(x) && !is.data.frame(x))
        stop("Invalid object 'x'. Must be a numeric matrix or data frame.")
    if (is.null(rownames(x))) 
        rownames(x) <- 1:nrow(x)
    if (is.null(colnames(x))) 
        colnames(x) <- paste("MV", 1:ncol(x), sep="")
    if (class(pls)!="plspm") 
        stop("argument 'pls' must be an object of class 'plspm'")
    if (nrow(x)!=pls$model[[6]])
        stop("arguments 'x' and 'pls' are incompatible")
    if (!is.factor(g)) stop("argument 'g' must be a factor")
    ng <- nlevels(g)
    if (ng > 2) stop("argument 'g' must contain only 2 levels") 
    if (!is.na(pmatch(method, "bootstrap"))) 
        method <- "bootstrap"
    METHODS <- c("bootstrap", "permutation")
    method <- pmatch(method, METHODS)
    if (is.na(method)) {
        warning("Invalid argument 'method'. Default 'method=bootstrp' is used.")   
        method <- "bootstrap"
    }
    if (is.null(reps) | length(reps)>1) reps<-100
    if (!is.numeric(reps) | floor(reps)<=0) reps<-100
    # ========================== INPUTS SETTING ==========================

    IDM <- pls$model[[1]]
    sets <- pls$model[[8]]
    scheme <- switch(pls$model[[3]], "1"="factor", "2"="centroid", "3"="path")
    plsr <- pls$model[[9]]
    lvs <- nrow(IDM)
    lvs.names <- rownames(IDM)
    mvs <- sum(unlist(lapply(sets, length)))
    blocks <- pls$model[[2]]
    blocklist <- sets
    for (k in 1:length(sets))
         blocklist[[k]] <- rep(k,blocks[k])
    blocklist <- unlist(blocklist)
    Mode <- pls$model[[4]]
    modes <- Mode
    modes[Mode=="Reflective"] <- "A"
    modes[Mode=="Formative"] <- "B"
    scaled <- pls$model[[5]]
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

    # ====================== Global model estimation =====================
    out.ws <- pls.weights(X, IDM, sets, scheme, blocklist)
    if (is.null(out.ws)) stop("The pls algorithm is non convergent") 
    cor.XY <- cor(X, X%*%out.ws[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
         w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    if (scaled) 
        Y.lvs <- round(X %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs), 4)
    else   
        Y.lvs <- round(DM %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs), 4)
    dimnames(Y.lvs) <- list(rownames(X), lvs.names)
    # Path coefficients
    pathmod <- pls.paths(IDM, Y.lvs, plsr)
    innmod <- pathmod[[1]]
    Path.global <- pathmod[[2]]
    R2.global <- pathmod[[3]]
    endo = rowSums(IDM)
    endo[endo!=0] <- 1  # vector indicating endogenous LVs
    path.labs <- NULL
    efs.labs <- NULL
    for (j in 1:lvs)
        for (i in j:lvs)
             if (IDM[i,j]==1) 
                 path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))
    path.global <- as.vector(Path.global[which(IDM==1)])
    names(path.global) <- path.labs
 
    # ====================== Group1 model estimation =====================
    g1.lab <- levels(g)[1]
    group1 <- which(g==levels(g)[1])
    wgs.g1 <- pls.weights(X[group1,], IDM, sets, scheme, blocklist)
    if (is.null(wgs.g1)) stop("The algorithm is non convergent") 
    cor.XY <- cor(X[group1,], X[group1,]%*%wgs.g1[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
         w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    if (scaled) 
        Y1.lvs <- round(X[group1,] %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs), 4)
    else   
        Y1.lvs <- round(DM[group1,] %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs), 4)
    dimnames(Y1.lvs) <- list(rownames(X[group1,]), lvs.names)
    # Path coefficients 
    pathmod.g1 <- pls.paths(IDM, Y1.lvs, plsr)
    innmod.g1 <- pathmod.g1[[1]]
    Path.g1 <- pathmod.g1[[2]]
    R2.g1 <- pathmod.g1[[3]]    
    path.g1 <- as.vector(Path.g1[which(IDM==1)])
    names(path.g1) <- path.labs

    # ====================== Group2 model estimation =====================
    g2.lab <- levels(g)[2]
    group2 <- which(g==levels(g)[2])
    wgs.g2 <- pls.weights(X[group2,], IDM, sets, scheme, blocklist)
    if (is.null(wgs.g2)) stop("The algorithm is non convergent") 
    cor.XY <- cor(X[group2,], X[group2,]%*%wgs.g2[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
         w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    if (scaled) 
        Y2.lvs <- round(X[group2,] %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs), 4)
    else   
        Y2.lvs <- round(DM[group2,] %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs), 4)
    dimnames(Y2.lvs) <- list(rownames(X[group2,]), lvs.names)
    # Path coefficients 
    pathmod.g2 <- pls.paths(IDM, Y2.lvs, plsr)
    innmod.g2 <- pathmod.g2[[1]]
    Path.g2 <- pathmod.g2[[2]]
    R2.g2 <- pathmod.g2[[3]]    
    path.g2 <- as.vector(Path.g2[which(IDM==1)])
    names(path.g2) <- path.labs

    # ====================== Group Comparison =====================
    dif.orig <- abs(path.g1-path.g2)
    nb <- round(reps)
    ng1 <- length(group1)
    ng2 <- length(group2)

    if (method==1)
    {
        BG1 <- matrix(0, nb, sum(IDM))
        BG2 <- BG1
        for (i in 1:nb)
        {
            samg1 <- sample(group1, ng1, replace=TRUE) 
            samg2 <- sample(group2, ng2, replace=TRUE)
            wgs.g1 <- pls.weights(X[samg1,], IDM, sets, scheme, blocklist)
            wgs.g2 <- pls.weights(X[samg2,], IDM, sets, scheme, blocklist)
            if (is.null(wgs.g1)) stop("Non convergence in bootstrap samples") 
            if (is.null(wgs.g2)) stop("Non convergence in bootstrap samples") 
            cor.XY <- cor(X[samg1,], X[samg1,]%*%wgs.g1[[2]])
            w.sig <- rep(NA,lvs)
            for (k in 1:lvs) 
                 w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
            if (scaled) 
                Y1.lvs <- round(X[samg1,] %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs), 4)
            else   
                Y1.lvs <- round(DM[samg1,] %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs), 4)        
            cor.XY <- cor(X[samg2,], X[samg2,]%*%wgs.g2[[2]])
            w.sig <- rep(NA,lvs)
            for (k in 1:lvs) 
                 w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
            if (scaled) 
                Y2.lvs <- round(X[samg2,] %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs), 4)
            else   
                Y2.lvs <- round(DM[samg2,] %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs), 4)
            pathmod.g1 <- pls.paths(IDM, Y1.lvs, plsr)
            paths.g1 <- pathmod.g1[[2]]    
            pathmod.g2 <- pls.paths(IDM, Y2.lvs, plsr)
            paths.g2 <- pathmod.g2[[2]]
            BG1[i,] <- as.vector(paths.g1[which(IDM==1)])
            BG2[i,] <- as.vector(paths.g2[which(IDM==1)])
        }    
        path.difs <- abs(apply(BG1,2,mean) - apply(BG2,2,mean))
        SE1 <- apply(BG1, 2, var)
        SE2 <- apply(BG2, 2, var)
        names(path.global) <- path.labs
        t.stat <- rep(NA, sum(IDM))
        k1 <- ((ng1-1)^2)/(ng1+ng2-2)
        k2 <- ((ng2-1)^2)/(ng1+ng2-2)
        k3 <- sqrt(1/ng1 + 1/ng2)
        for (i in 1:sum(IDM))         
             t.stat[i] <- path.difs[i] / (sqrt(k1*SE1[i]+k2*SE2[i]) * k3)        
        p.val <- pt(t.stat, ng1+ng2-2, lower.tail=FALSE)
        res <- round(cbind(path.global, path.g1, path.g2, dif.orig, 
                      t.stat, df=rep(ng1+ng2-2,sum(IDM)), p.val), 4)
        colnames(res) <- c("global", paste(rep("group",2),levels(g),sep="."), 
                             "diff.abs", "t.stat", "deg.fr", "p.value")
    } else
    {
        dif.perm <- matrix(0, nb, sum(IDM))
        for (i in 1:nb)# multigroup permutation
        {
            permu <- sample(1:(ng1+ng2), ng1+ng2)
            samg1 <- permu[1:ng1]
            samg2 <- permu[(ng1+1):(ng1+ng2)]
            wgs.g1 <- pls.weights(X[samg1,], IDM, sets, scheme, blocklist)
            wgs.g2 <- pls.weights(X[samg2,], IDM, sets, scheme, blocklist)
            if (is.null(wgs.g1)) stop("Non convergence in bootstrap samples") 
            if (is.null(wgs.g2)) stop("Non convergence in bootstrap samples") 
            cor.XY <- cor(X[samg1,], X[samg1,]%*%wgs.g1[[2]])
            w.sig <- rep(NA,lvs)
            for (k in 1:lvs) 
                 w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
            if (scaled) 
                Y1.lvs <- round(X[samg1,] %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs), 4)
            else   
                Y1.lvs <- round(DM[samg1,] %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs), 4)        
            cor.XY <- cor(X[samg2,], X[samg2,]%*%wgs.g2[[2]])
            w.sig <- rep(NA,lvs)
            for (k in 1:lvs) 
                 w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
            if (scaled) 
                Y2.lvs <- round(X[samg2,] %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs), 4)
            else   
                Y2.lvs <- round(DM[samg2,] %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs), 4)
            pathmod.g1 <- pls.paths(IDM, Y1.lvs, plsr)
            paths.g1 <- pathmod.g1[[2]]    
            pathmod.g2 <- pls.paths(IDM, Y2.lvs, plsr)
            paths.g2 <- pathmod.g2[[2]]
            pp1 <- as.vector(paths.g1[which(IDM==1)])
            pp2 <- as.vector(paths.g2[which(IDM==1)])
            dif.perm[i,] <- abs(pp1 - pp2)
        }   
        s.perm <- dif.orig 
        for (i in 1:sum(IDM))         
            s.perm[i] <- length(which(dif.perm[,i]>dif.orig[i]))
        p.val <- (1/(nb+1))*s.perm 
        res <- round(cbind(path.global, path.g1, path.g2, dif.orig, p.val), 4)
        colnames(res) <- c("global", paste(rep("group",2),levels(g),sep="."), 
                             "diff.abs", "p.value")  
    }
    met <- switch(method, "1"="bootstrap", "2"="permutation")
    settings <- c(scaled=scaled, scheme=scheme, method=met)
    res <- list(test=res, global=innmod, group1=innmod.g1, group2=innmod.g2, 
                settings=settings, reps=reps)
    class(res) <- "plspm.groups"
    return(res)
}

