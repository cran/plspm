`pls.efects` <-
function(Path)
{
    lvs <- nrow(Path)
    lvs.names <- rownames(Path)
    path.efects <- as.list(1:(lvs-1))
    path.efects[[1]] <- Path
    if (lvs == 2)
    {
        ind.paths <- matrix(c(0,0,0,0),2,2)
        total.paths <- Path
    }
    if (lvs > 2)
    {
        for (k in 2:(lvs-1))
            path.efects[[k]] <- round(path.efects[[k-1]] %*% Path, 4)
        ind.paths <- matrix(0, lvs, lvs)
        for (k in 2:length(path.efects))
            ind.paths <- ind.paths + path.efects[[k]]
        total.paths <- Path + ind.paths
    }
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

