MCMCdefs.RW.block <- list('RW_block' = quote({
    ## ## ************************************************************
    ## block together corresponging phi and gamma for each year with
    ## random walk sampler
    ## ***************************************************************
    customSpec <- configureMCMC(Rmodel)
    parms.phi <- Rmodel$expandNodeNames("phi")
    parms.gam <- Rmodel$expandNodeNames("gamma")
    parms.p <- Rmodel$expandNodeNames("p")
    phi.gam.p <- cbind(parms.phi, parms.gam, parms.p[-length(parms.p)])
    ## find node names for random effects
    for(i in 1:nrow(phi.gam.p)){
        customSpec$removeSamplers(phi.gam.p[i,], print=FALSE)
        customSpec$addSampler(target = phi.gam.p[i,],
                              type = "RW_block")
    }
    customSpec
}))



MCMCdefs.AFSS.block <- list('AFSS_block' = quote({
    ## ***************************************************************
    ## block together corresponging phi and gamma for each year with
    ## automated factor slice sampler
    ## ***************************************************************
    customSpec <- configureMCMC(Rmodel)
    parms.phi <- Rmodel$expandNodeNames("phi")
    parms.gam <- Rmodel$expandNodeNames("gamma")
    parms.p <- Rmodel$expandNodeNames("p")
    phi.gam.p <- cbind(parms.phi, parms.gam, parms.p[-length(parms.p)])
    ## find node names for random effects
    for(i in 1:nrow(phi.gam.p)){
        customSpec$removeSamplers(phi.gam.p[i,], print=FALSE)
        customSpec$addSampler(target = phi.gam.p[i,],
                              type = "AF_slice")
    }
    customSpec
}))



MCMCdefs.slice <- list('jags_like_nimble' = quote({
    ## ***************************************************************
    ## jags sampler specifications via nimble. Slice samplers on
    ## continous nodes.
    ## ***************************************************************
    customSpec <- configureMCMC(Rmodel)
    parms.phi <- Rmodel$expandNodeNames("phi")
    parms.gam <- Rmodel$expandNodeNames("gamma")
    parms.p <- Rmodel$expandNodeNames("p")
    phi.gam.p <- c(parms.phi, parms.gam, parms.p)
    ## find node names for random effects
    for(i in 1:length(phi.gam.p)){
        customSpec$removeSamplers(phi.gam.p[i], print=FALSE)
        customSpec$addSampler(target = phi.gam.p[i],
                              type = "slice")
    }
    customSpec
}))

