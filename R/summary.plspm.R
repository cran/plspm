`summary.plspm` <-
function(object, ...)
{
    ## $model <- list(IDM, blocks, scheme, Mode, scaled, obs=nrow(DM), boot.val)
    y <- object
    IDM <- y$model[[1]]
    blocks <- y$model[[2]]
    exo.endo <- rowSums(IDM)
    exo.endo[rowSums(IDM)==0] <- "Exogenous"
    exo.endo[rowSums(IDM)!=0] <- "Endogenous"
    blocklist <- as.list(1:sum(blocks))
    for (k in 1:length(blocks))
         blocklist[[k]] = rep(k,blocks[k])
    blocklist <- unlist(blocklist)
    inputs <- data.frame(Block=rownames(IDM), Type=exo.endo, 
                  NMVs=y$model[[2]], Mode=y$model[[4]])
    rownames(inputs) <- 1:length(exo.endo)
    lat.cor <- round(cor(y$latents), 4)
    res <- list(inputs=inputs, unidim=y$unidim, outer.mod=y$outer.mod, 
              outer.cor =y$outer.cor, inner.mod=y$inner.mod, latent.cor=lat.cor,
              inner.sum=y$inner.sum, gof=y$gof, effects=y$effects, xxx=y$model)
    if (y$model[[7]])
        res <- list(inputs=inputs, unidim=y$unidim, outer.mod=y$outer.mod, 
                outer.cor =y$outer.cor, inner.mod=y$inner.mod, latent.cor=lat.cor,
                inner.sum=y$inner.sum, gof=y$gof, effects=y$effects, boot=y$boot, xxx=y$model)
    class(res) <- "summary.plspm"
    res
}

