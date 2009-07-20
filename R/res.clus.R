`res.clus` <-
function(pls)
{
    # ========================= res.clus function ========================
    # Function to calculate communality and structural residuals that
    # will be used for the Response-Based Units Segmentation (REBUS)
    # in the "rebus" function
    # =========================== ARGUMENTS ==============================
    # pls: object of class "plspm"
    ## pls$model <- list(IDM, blocks, scheme, modes, 
    ##                   scaled, boot.val, plsr, obs, br)

    # ==================== Checking function arguments ===================
    if (class(pls)!="plspm") 
        stop("An object of class 'plspm' was expected")
    if (any(pls$model[[4]]!="A"))# checking reflective modes
        stop("REBUS only works for reflective modes")
    if (!pls$model[[5]])# checking scaled data
        stop("REBUS only works with scaled='TRUE'")

    # ========================== INPUTS SETTING ==========================
    IDM <- pls$model[[1]]# Inner Design Matrix
    blocks <- pls$model[[2]]# cardinality of blocks
    scheme <- pls$model[[3]]# inner weighting scheme
    modes <- pls$model[[4]]# measurement modes
    scaled <- pls$model[[5]]# type of scaling
    plsr <- pls$model[[7]]# pls-regression
    DM <- pls$data
    lvs <- nrow(IDM)
    lvs.names <- rownames(IDM)
    mvs <- sum(blocks)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    # data scaling (standardized data)
    X <- scale(DM)
    
    # ====================== computation of residuals =====================
    Y.lvs <- pls$latents# recovering LV scores from pls
    loads <- pls$loadings# recovering loadings from pls
    PaCo <- pls$path.coefs# recovering path coeffs from pls
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1  # indicator of endogenous LVs
    out.res <- DM# matrix for storing outer resids
    inn.res <- Y.lvs[,endo==1]# matrix for storing inner resids
    # computation of outer residuals
    for (j in 1:lvs)
    {
        q <- which(blocklist==j) 
        X.hat <- Y.lvs[,j] %*% t(loads[q])
        out.res[,q] <- X[,q] - X.hat# outer residuals
    }
    # computationi of inner residuals
    if (sum(endo)!=1)# more than 1 endogenous LV
        Y.hat <- Y.lvs %*% t(PaCo[endo==1,])        
    if (sum(endo)==1)# only 1 endogenous LV
        Y.hat <- Y.lvs %*% PaCo[endo==1,]        
    inn.res <- Y.lvs[,endo==1] - Y.hat# inner residuals
    
    # ====================== cluster analysis =====================
    # hierarchical cluster analysis with Ward method using function "hcluster"
    res <- cbind(out.res, inn.res)    
    res.clus <- hcluster(res, method="euclidean", diag=FALSE, upper=FALSE,
     link="ward", members=NULL, nbproc=2, doubleprecision=TRUE)
    # plot of the dendrogram
    plot(res.clus, main=c("REBUS", "Cluster Dendrogram of Outer and Inner Residuals"),
         hang=-1, cex.main=.9, cex.axis=.5, xlab="Hierarchical Clustering", 
         sub="Ward method", labels=FALSE)
    return(res.clus)
}

