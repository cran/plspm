`pls.loads` <-
function(X, Y.lvs, blocklist, modes)
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

