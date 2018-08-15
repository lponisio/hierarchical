MCMCdefs.RW.block <- list('RW_block' = quote({
    ## block together corresponging phi and gamma for each year with
    ## random walk sampler
    customSpec <- configureMCMC(Rmodel)
    parms.phi <- Rmodel$getNodeNames(
                            includeData = FALSE)[grepl("^phi",
                            Rmodel$getNodeNames(includeData = FALSE))]
    parms.gam <- Rmodel$getNodeNames(
                            includeData = FALSE)[grepl("^gamma",
                            Rmodel$getNodeNames(includeData = FALSE))]
    phi.gam <- cbind(parms.phi, parms.gam)[-1,]
    ## find node names for random effects
    for(i in 1:nrow(phi.gam)){
        customSpec$removeSamplers(phi.gam[i,], print=FALSE)
        customSpec$addSampler(target = phi.gam[i,],
                              type = "RW_block")
    }
    customSpec
}))



MCMCdefs.AFSS.block <- list('AFSS_block' = quote({
    ## block together corresponging phi and gamma for each year with
    ## automated factor slice sampler
    customSpec <- configureMCMC(Rmodel)
    parms.phi <- Rmodel$getNodeNames(
                            includeData = FALSE)[grepl("^phi",
                            Rmodel$getNodeNames(includeData = FALSE))]
    parms.gam <- Rmodel$getNodeNames(
                            includeData = FALSE)[grepl("^gamma",
                            Rmodel$getNodeNames(includeData = FALSE))]
    phi.gam <- cbind(parms.phi, parms.gam)[-1,]
    ## find node names for random effects
    for(i in 1:nrow(phi.gam)){
        customSpec$removeSamplers(phi.gam[i,], print=FALSE)
        customSpec$addSampler(target = phi.gam[i,],
                              type = "AF_slice")
    }
    customSpec
}))



MCMCdefs.slice <- list('jags_like_nimble' = quote({
    ## jags sampler specifications via nimble. Slice samplers on
    ## continous nodes.
    customSpec <- configureMCMC(Rmodel)
    parms.phi <- Rmodel$getNodeNames(
                            includeData = FALSE)[grepl("phi",
                            Rmodel$getNodeNames(includeData = FALSE))]
    parms.gam <- Rmodel$getNodeNames(
                            includeData = FALSE)[grepl("gamma",
                            Rmodel$getNodeNames(includeData = FALSE))]
    ## parms.mu <- Rmodel$getNodeNames(
    ##                         includeData = FALSE)[grepl("^mu",
    ##                          Rmodel$getNodeNames(includeData = FALSE))]
    ## parms.sigma <- Rmodel$getNodeNames(
    ##                         includeData = FALSE)[grepl("^sigma",
    ##                         Rmodel$getNodeNames(includeData = FALSE))]
    phi.gam <- c(parms.phi, parms.gam)
    phi.gam <- phi.gam[!grepl("tau", phi.gam)]
    phi.gam <- phi.gam[!grepl("lifted", phi.gam)]
    ## find node names for random effects
    for(i in 1:length(phi.gam)){
        customSpec$removeSamplers(phi.gam[i], print=FALSE)
        customSpec$addSampler(target = phi.gam[i],
                              type = "slice")
    }
    customSpec
}))


MCMCdefs.subsamp <- list('nimbleSubsamp' = quote({
    customSpec <- configureMCMC(Rmodel)
    customSpec$removeSamplers('z')
    customSpec$addSampler('z', type = 'sampler_latentSub',
                          control = list(leaveOutProportion = 0.6,
                                         control = list()))
    customSpec
}))
