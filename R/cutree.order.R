cutree.order <-
function(hclu, k=NULL, h=NULL)
{  
    coupe <- cutree(hclu, k=k, h=h)
    coupe.or <- coupe[hclu$order]
    coupe.out<- rep(NA,length(coupe))
    j <- 1 
    k <- coupe.or[1]
    for (i in 1:length(coupe))
    {
        if (coupe.or[i]==k) next
        else {
            coupe.out[which(coupe==k)] <- j
            j <- j + 1
            k <- coupe.or[i]
        }
    }
    coupe.out[is.na(coupe.out)] <- j
    names(coupe.out) <- names(coupe)
    coupe.out
}

