`pls.gof` <-
function(comu, R2, blocks, IDM)
{
    lvs <- nrow(IDM)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1  
    n.end <- sum(endo)
    # average of communalities
    comu.aveg <- rep(NA,lvs) 
    R2.aux <- rep(NA,n.end)
    aux <- 0
    for (j in 1:lvs)
        comu.aveg[j] <- mean(comu[blocklist==j]) 
    R2.aux <- R2[endo==1]
    gof <- sqrt(mean(comu.aveg)*mean(R2.aux))
    return(gof)
}

