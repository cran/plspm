`pls.unidim` <-
function(DM, blocks, modes)
{
    lvs <- length(blocks) 
    lvs.names <- names(blocks)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    Mode <- modes
    Mode[modes=="A"] <- "Reflective"
    Mode[modes=="B"] <- "Formative"   
    # Unidimensionality and GoF outer model
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
    res.uni <- list(unidim=unidim, gof.om=gof.om)
    return(res.uni)
}

