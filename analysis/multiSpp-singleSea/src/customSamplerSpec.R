

MCMCdefs.RW.block <- list('block_RW' = quote({
    ## *********************************************************************
    ##  block sampler for species random effects for each species, RW
    ## *********************************************************************
    customSpec <- configureMCMC(Rmodel)
    sigma.nodes <-  Rmodel$getNodeNames(stochOnly=TRUE)[grepl("sigma",
                                                              Rmodel$getNodeNames(stochOnly=TRUE))]
    if(length(sigma.nodes) > 0){
        customSpec$removeSamplers(sigma.nodes, print=FALSE)
        for(node in sigma.nodes){
            customSpec$addSampler(target = node,
                                  type = "slice")
        }
    }
    ## find node names of each species for random effects
    base.names <- c("a1", "a2", "a3", "a4", "b1", "b2",
                    "u.cato", "u.fcw", "v.cato", "v.fcw")
    exp.names.list <- list()
    for(bn in base.names){
        exp.names.list[[bn]] <- Rmodel$expandNodeNames(bn)
    }
    if(all(sapply(exp.names.list, length) == 1)){
        ## base.names <- c("a1", "a2", "a3", "a4", "b1", "b2",
        ##                 "cato.occ.mean", "fcw.occ.mean", "cato.det.mean", "fcw.det.mean")
        customSpec$removeSamplers(Rmodel$expandNodeNames(base.names), print=FALSE)
        customSpec$addSampler(target = Rmodel$expandNodeNames(base.names),
                              type = "RW_block")
    } else{
        for(i in 1:length(exp.names.list[[1]])){
            blocknames <- unlist(lapply(exp.names.list, function(x) x[i]))
            customSpec$removeSamplers(blocknames, print=FALSE)
            customSpec$addSampler(target = blocknames,
                                  type = "RW_block")
        }
    }
    customSpec
}))



MCMCdefs.AFSS.block <- list('block_AFSS' = quote({
    ## *********************************************************************
    ##  block sampler for species random effects for each species, AFSS
    ## *********************************************************************
    customSpec <- configureMCMC(Rmodel)
    sigma.nodes <-  Rmodel$getNodeNames(stochOnly=TRUE)[grepl("sigma",
                                                              Rmodel$getNodeNames(stochOnly=TRUE))]
    if(length(sigma.nodes) > 0){
        customSpec$removeSamplers(sigma.nodes, print=FALSE)
        for(node in sigma.nodes){
            customSpec$addSampler(target = node,
                                  type = "slice")
        }
    }
    ## find node names of each species for random effects
    base.names <- c("a1", "a2", "a3", "a4", "b1", "b2",
                    "u.cato", "u.fcw", "v.cato", "v.fcw")
    exp.names.list <- list()
    for(bn in base.names){
        exp.names.list[[bn]] <- Rmodel$expandNodeNames(bn)
    }
    if(all(sapply(exp.names.list, length) == 1)){
        ## base.names <- c("a1", "a2", "a3", "a4", "b1", "b2",
        ##                 "cato.occ.mean", "fcw.occ.mean", "cato.det.mean", "fcw.det.mean")
        customSpec$removeSamplers(Rmodel$expandNodeNames(base.names), print=FALSE)
        customSpec$addSampler(target = Rmodel$expandNodeNames(base.names),
                              type = "AF_slice")
    } else{
        for(i in 1:length(exp.names.list[[1]])){
            blocknames <- unlist(lapply(exp.names.list, function(x) x[i]))
            customSpec$removeSamplers(blocknames, print=FALSE)
            customSpec$addSampler(target = blocknames,
                                  type = "AF_slice")
        }
    }
    customSpec
}))



MCMCdefs.slice <- list('jags_like_nimble' = quote({
    ## ***************************************************************
    ## jags sampler specifications via nimble. Slice samplers on
    ## continous nodes.
    ## ***************************************************************
    customSpec <- configureMCMC(Rmodel)
    sigma.nodes <-  Rmodel$getNodeNames(stochOnly=TRUE)[grepl("sigma",
                                                              Rmodel$getNodeNames(stochOnly=TRUE))]
    if(length(sigma.nodes) > 0){
        customSpec$removeSamplers(sigma.nodes, print=FALSE)
        for(node in sigma.nodes){
            customSpec$addSampler(target = node,
                                  type = "slice")
        }
    }

    mu.nodes <-  Rmodel$getNodeNames(stochOnly=TRUE)[grepl("mu",
                                                           Rmodel$getNodeNames(stochOnly=TRUE))]
    if(length(mu.nodes) > 0){
        customSpec$removeSamplers(mu.nodes, print=FALSE)
        for(node in mu.nodes){
            customSpec$addSampler(target = node,
                                  type = "slice")
        }
    }
    ## find node names of each species for random effects
    base.names <- c("a1", "a2", "a3", "a4", "b1", "b2")
    cont.nodes <- Rmodel$expandNodeNames(base.names)
    for(node in cont.nodes){
        customSpec$removeSamplers(node, print=FALSE)
        customSpec$addSampler(target = node,
                              type = "slice")
    }
    customSpec
}))

## ## remove the samples, add block samplers
## MCMCdefs.RW.block2 <- list('RW_block2' = quote({
##     ## *********************************************************************
##     ## ##  block sampler for species random effects for each
##     ## "type", RW
##     ## *********************************************************************

##     customSpec <- configureMCMC(Rmodel)
##     ## find node names for random effects
##     sp.parms.a <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^a",
##                   Rmodel$getNodeNames(includeData = FALSE))]
##     sp.parms.b <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^b",
##                   Rmodel$getNodeNames(includeData = FALSE))]
##     sp.parms.u <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^u",
##                   Rmodel$getNodeNames(includeData = FALSE))]
##     sp.parms.v <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^v",
##                   Rmodel$getNodeNames(includeData = FALSE))]
##     customSpec$removeSamplers(c(sp.parms.a, sp.parms.b, sp.parms.u,
##                                 sp.parms.v),
##                               print=FALSE)

##     customSpec$addSampler(target = sp.parms.a, type = "RW_block")
##                           customSpec$addSampler(target =
##                           sp.parms.b, type = "RW_block")
##                           customSpec$addSampler(target =
##                           sp.parms.u, type = "RW_block")
##                           customSpec$addSampler(target =
##                           sp.parms.v, type = "RW_block") customSpec
##                           }))



## MCMCdefs.AFSS.block2 <- list('AFSS_block2' = quote({
##     ## *********************************************************************
##     ## ##  block sampler for species random effects for each
##     ## "type", AFSS
##     ## *********************************************************************

##     customSpec <- configureMCMC(Rmodel)
##     ## find node names for random effects
##     sp.parms.a <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^a",
##                   Rmodel$getNodeNames(includeData = FALSE))]
##     sp.parms.b <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^b",
##                   Rmodel$getNodeNames(includeData = FALSE))]
##     sp.parms.u <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^u",
##                   Rmodel$getNodeNames(includeData = FALSE))]
##     sp.parms.v <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^v",
##                   Rmodel$getNodeNames(includeData = FALSE))]
##     customSpec$removeSamplers(c(sp.parms.a, sp.parms.b, sp.parms.u,
##                                 sp.parms.v),
##                               print=FALSE)

##     customSpec$addSampler(target = sp.parms.a,
##                           type = "AF_slice")
##     customSpec$addSampler(target = sp.parms.b,
##                           type = "AF_slice")
##     customSpec$addSampler(target = sp.parms.u,
##                           type = "AF_slice")
##     customSpec$addSampler(target = sp.parms.v,
##                           type = "AF_slice")
##     customSpec
## }))
