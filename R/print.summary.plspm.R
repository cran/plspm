`print.summary.plspm` <-
function(x, ...)
{
    ## x$xxx <- list(IDM, blocks, scheme, modes, scaled, boot.val, plsr, obs, br) 
    Scheme <- x$xxx[[3]]
    boot.sam <- if(is.null(x$xxx[[9]])) "NULL" else x$xxx[[9]]
    if (x$xxx[[5]]) Scale="Standardized Data" else Scale="Raw Data"
    cat("\n")
    cat("PARTIAL LEAST SQUARES PATH MODELING (PLS-PM)", "\n\n")
    cat("----------------------------------------------", "\n")    
    cat("MODEL SPECIFICATION", "\n")
    cat("1   Number of Cases", "\t", x$xxx[[8]], "\n")
    cat("2   Latent Variables", "\t", nrow(x$xxx[[1]]), "\n")
    cat("3   Manifest Variables", "\t", sum(x$xxx[[2]]), "\n")
    cat("4   Scale of Data", "\t", Scale, "\n")
    cat("5   Weighting Scheme", "\t", Scheme, "\n")
    cat("6   Bootstrapping", "\t", x$xxx[[6]], "\n")
    cat("7   Bootstrap samples", "\t", boot.sam, "\n")
    cat("\n")
    cat("---------------------------------------------------", "\n")    
    cat("BLOCKS DEFINITION", "\n")
    print(x$inputs, print.gap=3)
    cat("\n")
    cat("---------------------------------------------------", "\n")    
    cat("BLOCKS UNIDIMENSIONALITY","\n")
    print(x$unidim, print.gap=2)
    cat("\n")
    cat("----------------------------------------------------------", "\n\n")    
    cat("OUTER MODEL","\n")
    OM <- NULL
    om.labs <- NULL
    for (k in 1:length(x$outer.mod)) {        
        OM <- rbind(OM, rep(NA,ncol(x$outer.mod[[k]])), x$outer.mod[[k]])
        om.labs <- c(om.labs, names(x$outer.mod)[k], 
                paste(rep(" ",nrow(x$outer.mod[[k]])), rownames(x$outer.mod[[k]])))
    }
    rownames(OM) <- om.labs
    print(OM, na.print="", print.gap=2)
    cat("\n")
    cat("----------------------------------------------------------", "\n")    
    cat("CORRELATIONS BETWEEN MVs AND LVs","\n")
    Cros <- NULL
    for (k in 1:length(x$outer.mod)) 
        Cros <- rbind(Cros, rep(NA,ncol(x$outer.cor[[k]])), x$outer.cor[[k]])
    rownames(Cros) <- om.labs
    print(Cros, na.print="", print.gap=2)
    cat("\n")
    cat("----------------------------------------------------------", "\n")    
    cat("INNER MODEL","\n")
    print(x$inner.mod, print.gap=3)
    cat("----------------------------------------------------------", "\n")    
    cat("CORRELATIONS BETWEEN LVs","\n")
    print(x$latent.cor, print.gap=2)
    cat("\n")
    cat("----------------------------------------------------------", "\n")    
    cat("SUMMARY INNER MODEL","\n")
    print(x$inner.sum, print.gap=2)
    cat("\n")
    cat("----------------------------------------------------------", "\n") 
    cat("GOODNESS-OF-FIT","\n")
    print(x$gof, print.gap=2)
    cat("\n")
    cat("----------------------------------------------------------", "\n")        
    cat("TOTAL EFFECTS","\n")
    print(x$effects, print.gap=2, digits=4)
    if (length(x)==11)
    {
        cat("\n")
        cat("---------------------------------------------------------", "\n")    
        cat("BOOTSTRAP VALIDATION", "\n")
        for (i in 1:length(x$boot))
        {
             cat(names(x$boot)[i], "\n")
             print(x$boot[[i]], print.gap=2)
             cat("\n")
        }
    }      
    invisible(x)
}

