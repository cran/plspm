print.local.models <-
function(x, ...)
{
    ## x$model = list(IDM, blocks, scheme, modes, scaled, boot.val, plsr, obs, br) 
    model <- x$glob.model$model
    boot.sam <- if(is.null(model[[9]])) "NULL" else model[[9]]
    if (model[[5]]) Scale="Standardized Data" else Scale="Raw Data"
    n.clus <- length(x) - 1
    Name <- c("$glob.model", paste(rep("$loc.model",n.clus), 1:n.clus, sep="."))
    Description <- c("global model", paste(rep("local model class",n.clus), 1:n.clus, sep=" "))
    res1 <- cbind(Name, Description)    
    cat("This function calculates PLS-PM for global and local models", "\n")
    cat("-----------------------------------------------------------", "\n")    
    cat("  Number of classes", "\t", length(x)-1,"\n")
    cat("  Scale of Data", "\t", Scale, "\n")
    cat("  Weighting Scheme", "\t", model[[3]], "\n")
    cat("  Bootstrapping", "\t", model[[6]], "\n")
    cat("  Bootstrap samples", "\t", boot.sam, "\n\n")
    cat("\n")
    cat("PLS-PM models in the following objects", "\n")
    rownames(res1) <- 1:nrow(res1)
    print(res1, print.gap=3, justify="right")
}

