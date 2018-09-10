MCMCdefs.RW.block <- list('RW_block' = quote({
    ## *********************************************************************
    ##  block sampler for intercept and other covariates
    ## *********************************************************************
    customSpec <- configureMCMC(Rmodel)
    customSpec$removeSamplers(c("beta0","beta"))
    customSpec$addSampler(c("beta0","beta"), type = "RW_block")
    customSpec$removeSamplers(c("mean.p", "alpha"))
    customSpec$addSampler(c("mean.p", "alpha"), type = "RW_block")
    customSpec
}))



MCMCdefs.AFSS.block <- list('AFSS_block' = quote({
    ## *********************************************************************
    ##  block sampler for species random effects for each species, AFSS
    ## *********************************************************************
    customSpec <- configureMCMC(Rmodel)
    customSpec <- configureMCMC(Rmodel)
    customSpec$removeSamplers(c("beta0","beta"))
    customSpec$addSampler(c("beta0","beta"), type = "AF_slice")
    customSpec$removeSamplers(c("mean.p", "alpha"))
    customSpec$addSampler(c("mean.p", "alpha"), type = "AF_slice")
    customSpec
}))



MCMCdefs.slice <- list('jags_like_nimble' = quote({
    ## ***************************************************************
    ## jags sampler specifications via nimble. Slice samplers on
    ## continous nodes.
    ## ***************************************************************
    customSpec <- configureMCMC(Rmodel)
    ## find node names of each species for random effects
    base.names <- c("alpha", "eps", "beta")
    all.nodes <- Rmodel$getNodeNames(stochOnly=TRUE)
    cont.nodes <- all.nodes[unlist(sapply(base.names,
                                                      function(x) grep(x, all.nodes)))]
    for(node in cont.nodes){
        customSpec$removeSamplers(node, print=FALSE)
        customSpec$addSampler(target = node,
                              type = "slice")
    }
    customSpec
}))

