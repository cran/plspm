.rec.hclust <-
function(index, lwd=1, lty=1, col="black")
{
    # index: index of the current tree to draw
    members <- get('members', envir= ._a2r_envir) 
    bottom  <- get('bottom',  envir= ._a2r_envir) 
    if (index<0){ # it is a leaf
        if(is.null(members)){
           ._a2r_counter <<- ._a2r_counter + 1
           return(list(x=._a2r_counter, n=1))
        }
        else{
            cc <- ._a2r_counter
            mm <- members[-index]
            polygon(x=c(cc, cc+mm/2, cc+mm), y=c(bottom, 0, bottom),
                    col=col, border = col, lwd=lwd)
            ._a2r_counter <<- ._a2r_counter + mm
            return(list(x=cc+mm/2, n=mm))
        }
    }    
    h.m   <- ._a2r_hclu$height[index]
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~ do left
    index.l  <- ._a2r_hclu$merge[index,1]    
    h.l <- if(index.l<0) 0 else ._a2r_hclu$height[index.l]
    if (h.l<._a2r_height_cut & h.m > ._a2r_height_cut){
        ._a2r_group <<- ._a2r_group + 1
        col.l <- get("col.down",envir=._a2r_envir)[._a2r_group]
        lwd.l <- get("lwd.down",envir=._a2r_envir)
        lty.l <- get("lty.down",envir=._a2r_envir)
    }
    else{
        col.l <- col
        lwd.l <- lwd
        lty.l <- lty
    }
    out.l   <- .rec.hclust(index.l, col=col.l, lty=lty.l, lwd=lwd.l)
    x.l     <- out.l$x
    n.l     <- out.l$n
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~ do right
    index.r  <- ._a2r_hclu$merge[index,2]
    h.r <- if(index.r<0) 0 else ._a2r_hclu$height[index.r]
    if (h.r<._a2r_height_cut & h.m > ._a2r_height_cut){
        ._a2r_group <<- ._a2r_group + 1
        col.r <- get("col.down",envir=._a2r_envir)[._a2r_group]
        lwd.r <- get("lwd.down",envir=._a2r_envir)
        lty.r <- get("lty.down",envir=._a2r_envir)
    }
    else{
        col.r <- col
        lwd.r <- lwd
        lty.r <- lty
    }
    out.r   <- .rec.hclust(index.r, col=col.r, lty=lty.r, lwd=lwd.r)
    x.r     <- out.r$x
    n.r     <- out.r$n
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~ draw what you have to draw    
    type <- get("type",envir=._a2r_envir)
    x.m  <- (x.r + x.l) / 2  
    n    <- n.r + n.l
    x.b  <- (n.r * x.r + n.l * x.l) / n      
    knot.pos <- get("knot.pos",envir=._a2r_envir)     
    x <- switch(knot.pos, mean=x.m, left=x.l, right= x.r,
            random = x.l + runif(1)*(x.r-x.l), bary=x.b)
            
    if (type=="rectangle"){
        segments(x0  = c(x.l, x.l, x.r),
                 x1  = c(x.l, x.r, x.r),
                 y0  = c(h.l, h.m, h.r),
                 y1  = c(h.m, h.m, h.m),
                 col = col,
                 lty = lty,
                 lwd = lwd)
    }
    if (type =="triangle"){
        segments(x0  = c(x.l, x.r),
                 x1  = c(x  , x),
                 y0  = c(h.l, h.r),
                 y1  = c(h.m, h.m),
                 col = col,
                 lty = lty,
                 lwd = lwd)
    }
                        
    list(x=x, n=n)
}

