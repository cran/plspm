`print.plspm` <-
function(x, ...)
{
    Name <- c("$unidim", "$outer.mod", "$inner.mod", "$latents", "$scores", "$out.weights",  
             "$loadings", "$path.coefs", "$r.sqr", "$outer.cor", "$inner.sum", "$gof", "$effects")
    Description <- c("unidimensionality", "outer model", "inner model", "scaled LVs (variance=1)", 
       "re-scaled LVs for scaled=FALSE", "outer weights", "measurement loadings", 
       "path coefficients", "R-squared", "outer correlations", "summary inner model", 
       "goodnes-of-fit", "total effects")
    if (length(x)==15) {
       Name <- c(Name, "$boot", "$boot$weights", "$boot$loadings", "$boot$paths",
                 "$boot$efects", "$boot$rsq")
       Description <- c(Description, "bootstrap results", "bootstrap weights", 
              "bootstrap loadings", "bootstrap path.coefs", "bootstrap total effects", "bootstrap R2")
    }
    res1 <- cbind(Name, Description)
    cat("\n")
    cat("PARTIAL LEAST SQUARES PATH MODELING (PLS-PM)", "\n")
    cat("----------------------------------------------", "\n")    
    cat("Results available in the following objects:", "\n\n")
    rownames(res1) <- 1:nrow(res1)
    print(res1, print.gap=3, justify="right")
    cat("----------------------------------------------", "\n")
    cat("You can also use the function 'summary'", "\n\n")    
    invisible(x)
}

