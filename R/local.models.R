local.models <-
function(pls, y, scheme=NULL, scaled=NULL, boot.val=FALSE, br=NULL)
{
    # ======================== local.models function ======================
    # Function to calculate PLS-PM for global and local models
    # using a "rebus" object or a vector (with memberships or categories)
    # =========================== arguments ===============================
    # pls: object of class "plspm"
    # y: can be an object of class "rebus", a numeric vector, or a factor
    # scheme: a character string indicating the inner weighting scheme 
    #         to be used: "factor", "centroid", or "path"
    # scaled: a logical value indicating whether scale data is performed
    # boot.val:a logical value indicating whether bootstrap validation is done 
    # br: an integer indicating the number of bootstraps resamples, used 
    #     only when boot.val=TRUE, (100 <= br <= 1000)

    # ==================== Checking function arguments ====================
    if (class(pls)!="plspm") 
        stop("argument 'pls' must be an object of class 'plspm'")
    if (!is.element(class(y), c("rebus","integer","factor")))   
        stop("argument 'y' must be of class 'rebus', 'integer' or 'factor'")
    if (class(y)=="rebus") {
        if (length(y$segments)!=nrow(pls$data))
            stop("arguments 'pls' and 'y' are incompatible")
    } else {
        if (length(y)!=nrow(pls$data))
            stop("arguments 'pls' and 'y' are incompatible")
    }
    if (is.null(scheme))
        scheme <- pls$model[[3]]
    if (is.null(pmatch(scheme, "centroid"))) 
        scheme <- "centroid"
    SCHEMES <- c("centroid", "factor")
    scheme <- pmatch(scheme, SCHEMES)
    if (is.na(scheme)) {
        warning("Invalid argument 'scheme'. Default 'scheme=centroid' is used.")   
        scheme <- "centroid"
    }
    if (is.null(scaled) || !is.logical(scaled))
        scaled <- pls$model[[5]]# type of scaling
    if (!is.logical(boot.val)) {
        warning("Invalid argument 'boot.val'. No bootstrap validation is done.")
        boot.val <- FALSE
    }   
    if (boot.val) {
        if (!is.null(br)) {        
            if (mode(br)!="numeric" || length(br)!=1 || (br%%1)!=0 ||
                br<100 || br>1000) {
                warning("Invalid argument 'br'. Default 'br=100' is used.")   
                br <- 100
            } 
        } else
            br <- 100
    }

    # ========================== INPUTS SETTING ==========================
    IDM <- pls$model[[1]]# Inner Design Matrix
    blocks <- pls$model[[2]]# cardinality of blocks
    modes <- pls$model[[4]]# measurement modes    
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
    end.ind <- cumsum(blocks)
    ini.ind <- cumsum(blocks) - blocks + 1
    new.sets <- as.list(1:lvs)
    for (j in 1:lvs)
        new.sets[[j]] <- ini.ind[j]:end.ind[j]
    if (class(y)=="rebus") {
        segments <- as.factor(y$segments)
    } else {
        segments <- as.factor(y)
    }
    n.clus <- length(table(segments))

    # ============ final models computation (global and local models) ============
    skem <- switch(scheme, "centroid"="centroid", "factor"="factor")
    final.mod <- as.list(1:(n.clus+1))# final plspm models
    for (k in 1:(n.clus+1))
    {
        if (k==1) {
            # global model
            X <- DM
            final.mod[[1]] <- plspm(X, IDM, new.sets, modes, skem, scaled, boot.val, br)
        } else
        {
            units.k <- which(segments==levels(segments)[k-1])
            # local models
            X.k <- DM[units.k,]
            final.mod[[k]] <- plspm(X.k, IDM, new.sets, modes, skem, scaled, boot.val, br)
        }
    }
    names(final.mod) <- c("glob.model",paste(rep("loc.model",n.clus), 1:n.clus, sep="."))
    class(final.mod) <- "local.models"
    return(final.mod)
}

