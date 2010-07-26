print.summary.plspm <-
function(x, ...)
{
    ## x$xxx <- list(IDM, blocks, scheme, modes, scaled, boot.val, plsr, obs, br, tol, iter, n.iter) 
    ## x.fit$xxx <- list(IDM, blocks, scheme, modes, scaled, obs, tol, iter, n.iter) 
    Scheme <- x$xxx[[3]]
    if (x$xxx[[5]]) Scale="Standardized Data" else Scale="Raw Data"
    cat("\n")
    cat("PARTIAL LEAST SQUARES PATH MODELING (PLS-PM)", "\n\n")
    cat("----------------------------------------------------------", "\n")    
    cat("MODEL SPECIFICATION", "\n")
    cat("1   Number of Cases", "\t", x$xxx$obs, "\n")
    cat("2   Latent Variables", "\t", nrow(x$xxx[[1]]), "\n")
    cat("3   Manifest Variables", "\t", sum(x$xxx[[2]]), "\n")
    cat("4   Scale of Data", "\t", Scale, "\n")
    cat("5   Weighting Scheme", "\t", Scheme, "\n")
    cat("6   Tolerance Crit", "\t", x$xxx$tol, "\n")
    cat("7   Max Num Iters", "\t", x$xxx$iter, "\n")
    cat("8   Convergence Iters", "\t", x$xxx$n.iter, "\n")
    if (length(x$xxx)>9)
    {
        boot.sam <- if(is.null(x$xxx[[9]])) "NULL" else x$xxx[[9]]
        cat("9   Paths by PLS-R", "\t", x$xxx[[7]], "\n")
        cat("10  Bootstrapping", "\t", x$xxx[[6]], "\n")
        cat("11  Bootstrap samples", "\t", boot.sam, "\n")
    }
    cat("\n")
    cat("----------------------------------------------------------", "\n")    
    cat("BLOCKS DEFINITION", "\n")
    print(x$inputs, print.gap=3)
    cat("\n")
    cat("----------------------------------------------------------", "\n") 
    if (length(x$xxx)>9) 
    {   
        cat("BLOCKS UNIDIMENSIONALITY","\n")
        print(x$unidim, print.gap=2, digits=3)
        cat("\n")
        cat("----------------------------------------------------------", "\n")    
    }
    cat("OUTER MODEL","\n")
    OM <- NULL
    om.labs <- NULL
    for (k in 1:length(x$outer.mod)) {        
        OM <- rbind(OM, rep(NA,ncol(x$outer.mod[[k]])), x$outer.mod[[k]])
        om.labs <- c(om.labs, names(x$outer.mod)[k], 
                paste(rep(" ",nrow(x$outer.mod[[k]])), rownames(x$outer.mod[[k]])))
    }
    rownames(OM) <- om.labs
    print(OM, na.print="", print.gap=2, digits=3)
    cat("\n")
    cat("----------------------------------------------------------", "\n")    
    if (length(x$xxx)>9)
    {
        cat("CORRELATIONS BETWEEN MVs AND LVs","\n")
        Cros <- NULL
        for (k in 1:length(x$outer.mod)) 
            Cros <- rbind(Cros, rep(NA,ncol(x$outer.cor[[k]])), x$outer.cor[[k]])
        rownames(Cros) <- om.labs
        print(Cros, na.print="", print.gap=2, digits=3)
        cat("\n")
        cat("----------------------------------------------------------", "\n")    
    }
    cat("INNER MODEL","\n")
    print(x$inner.mod, print.gap=3, digits=3)
    if (length(x$xxx)>9)
    {
        cat("----------------------------------------------------------", "\n")    
        cat("CORRELATIONS BETWEEN LVs","\n")
        print(x$latent.cor, print.gap=2, digits=3)
        cat("\n")
        cat("----------------------------------------------------------", "\n")    
        cat("SUMMARY INNER MODEL","\n")
        print(x$inner.sum, print.gap=2, digits=3)
        cat("\n")
        cat("----------------------------------------------------------", "\n") 
        cat("GOODNESS-OF-FIT","\n")
        print(x$gof, print.gap=2, digits=4)
        cat("\n")
        cat("----------------------------------------------------------", "\n")        
        cat("TOTAL EFFECTS","\n")
        print(x$effects, print.gap=2, digits=3)
        if (length(x)==11)
        {
            cat("\n")
            cat("---------------------------------------------------------", "\n")    
            cat("BOOTSTRAP VALIDATION", "\n")
            for (i in 1:length(x$boot))
            {
                 cat(names(x$boot)[i], "\n")
                 print(x$boot[[i]], print.gap=2, digits=3)
                 cat("\n")
            }
        }      
    }
    invisible(x)
}

