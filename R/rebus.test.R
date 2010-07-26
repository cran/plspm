rebus.test <-
function(pls, reb)
{
    # ======================= rebus.test function ========================
    # Function to perform tests for multi-group comparison from the 
    # classes obtained by REBUS
    # =========================== arguments ==============================
    # pls: object of class "plspm"
    # reb: object of class "rebus" 

    # ==================== Checking function arguments ===================
    if (class(pls)!="plspm") 
        stop("argument 'pls' must be an object of class 'plspm'")
    if (any(pls$model[[4]]!="A"))# checking reflective modes
        stop("REBUS only works for reflective modes")
    if (!pls$model[[5]])# checking scaled data
        stop("REBUS only works with scaled='TRUE'")
    if (class(reb)!="rebus") 
        stop("argument 'reb' must be an object of class 'rebus'")
    if (length(reb$segments)!=nrow(pls$data))
        stop("arguments 'pls' and 'reb' are incompatible")
    if (length(table(reb$segments))>6)
        stop("the number of classes in 'rebus.test' is limited to 6")

    # ========================== INPUTS SETTING ==========================
    IDM <- pls$model$IDM# Inner Design Matrix
    blocks <- pls$model$blocks# cardinality of blocks
    scheme <- pls$model$scheme# inner weighting scheme
    modes <- pls$model$modes# measurement modes
    scaled <- pls$model$scaled# type of scaling
    plsr <- pls$model$plsr# pls-regression
    tol <- pls$model$tol# tolerance criterion
    iter <- pls$model$iter# max num iterations
    DM <- pls$data
    lvs <- nrow(IDM)
    lvs.names <- rownames(IDM)
    mvs <- sum(blocks)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    # data scaling (standardized data)
    sd.X <- sqrt((nrow(DM)-1)/nrow(DM)) * apply(DM, 2, sd)
    X <- scale(DM, scale=sd.X)
    n.clus <- length(table(reb$segments))
    # multi-group comparison
    ic <- NULL
    ec <- NULL
    for (i in 1:(n.clus-1))
    {
        ic <- c(ic, rep(i,(n.clus-i)))
        ec <- c(ec, seq((i+1),n.clus))
    }
    gp.index <- cbind(ic,ec)
    gp.test <- as.list(1:nrow(gp.index))
    for (i in 1:nrow(gp.index))
    { 
        a <- which(reb$segments%in%gp.index[i,])
        g <- as.factor(reb$segments[a])
        gp.test[[i]] <- .pls.locals.test(DM[a,], pls, g)
    }
    names(gp.test) <- paste(rep("test",nrow(gp.index)),gp.index[,1],gp.index[,2],sep="_")
    class(gp.test) <- "rebus.test"
    return(gp.test)
}

