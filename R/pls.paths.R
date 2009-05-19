`pls.paths` <-
function(IDM, Y.lvs, plsr)
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

