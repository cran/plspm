`pls.loads` <-
function(X, Y.lvs, blocks)
{
    lvs <- length(blocks)
    mvs <- ncol(X)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    loads <- rep(NA, mvs)
    comu <- rep(NA, mvs)
    for (j in 1:lvs)
        loads[blocklist==j] <- cor(X[,blocklist==j], Y.lvs[,j])
    comu <- loads^2
    names(loads) <- colnames(X)  
    names(comu) <- colnames(X)
    res.loads <- list(loads, comu)
    return(res.loads)
} # end function 'pls.loads'

