

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
    base.names.phi <- c("phi.0", "phi.hr.area", "phi.nat.area",
                        "phi.fra")
    base.names.gam <- c("gam.0", "gam.hr.area", "gam.nat.area",
                        "gam.fra")

    exp.names.list.phi <- list()
    for(bn in base.names.phi){
        exp.names.list.phi[[bn]] <- Rmodel$expandNodeNames(bn)
    }

    exp.names.list.gam <- list()
    for(bn in base.names.gam){
        exp.names.list.gam[[bn]] <- Rmodel$expandNodeNames(bn)
    }
    if(all(sapply(exp.names.list.phi, length) == 1)){
        customSpec$removeSamplers(Rmodel$expandNodeNames(base.names.phi), print=FALSE)
        customSpec$addSampler(target = Rmodel$expandNodeNames(base.names.phi),
                              type = "RW_block")
        customSpec$removeSamplers(Rmodel$expandNodeNames(base.names.gam), print=FALSE)
        customSpec$addSampler(target = Rmodel$expandNodeNames(base.names.gam),
                              type = "RW_block")
    } else{
        for(i in 1:length(exp.names.list.phi[[1]])){
            blocknames.phi <- unlist(lapply(exp.names.list.phi,
                                            function(x) x[i]))
            blocknames.gam <-  unlist(lapply(exp.names.list.gam,
                                             function(x) x[i]))

            customSpec$removeSamplers(blocknames.phi, print=FALSE)
            customSpec$addSampler(target = blocknames.phi,
                                  type = "RW_block")
            customSpec$removeSamplers(blocknames.gam, print=FALSE)
            customSpec$addSampler(target = blocknames.gam,
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
    base.names.phi <- c("phi.0", "phi.hr.area", "phi.nat.area",
                        "phi.fra")
    base.names.gam <- c("gam.0", "gam.hr.area", "gam.nat.area",
                        "gam.fra")

    exp.names.list.phi <- list()
    for(bn in base.names.phi){
        exp.names.list.phi[[bn]] <- Rmodel$expandNodeNames(bn)
    }

    exp.names.list.gam <- list()
    for(bn in base.names.gam){
        exp.names.list.gam[[bn]] <- Rmodel$expandNodeNames(bn)
    }


    if(all(sapply(exp.names.list.phi, length) == 1)){
        customSpec$removeSamplers(Rmodel$expandNodeNames(base.names.phi), print=FALSE)
        customSpec$addSampler(target = Rmodel$expandNodeNames(base.names.phi),
                              type = "AF_slice")
        customSpec$removeSamplers(Rmodel$expandNodeNames(base.names.gam), print=FALSE)
        customSpec$addSampler(target = Rmodel$expandNodeNames(base.names.gam),
                              type = "AF_slice")
    } else{
        for(i in 1:length(exp.names.list.phi[[1]])){
            blocknames.phi <- unlist(lapply(exp.names.list.phi,
                                            function(x) x[i]))
            blocknames.gam <-  unlist(lapply(exp.names.list.gam,
                                             function(x) x[i]))

            customSpec$removeSamplers(blocknames.phi, print=FALSE)
            customSpec$addSampler(target = blocknames.phi,
                                  type = "AF_slice")
            customSpec$removeSamplers(blocknames.gam, print=FALSE)
            customSpec$addSampler(target = blocknames.gam,
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
    ## find node names of each species for random effects
    base.names <- c("phi.0", "phi.hr.area", "phi.nat.area",
                    "phi.fra",
                    "phi.k", "phi.B",
                    "phi.hr.area.fra",
                    "phi.nat.area.fra",
                    "phi.hr.area.k",
                    "phi.nat.area.k",
                    "phi.hr.area.B",
                    "phi.nat.area.B",
                    "gam.0", "gam.hr.area", "gam.nat.area",
                    "gam.fra",
                    "gam.k", "gam.B",
                    "gam.hr.area.fra",
                    "gam.nat.area.fra",
                    "gam.hr.area.k",
                    "gam.nat.area.k",
                    "gam.hr.area.B",
                    "gam.nat.area.B",
                    "p.0", "p.day.1",
                    "p.day.2")
    cont.nodes <- Rmodel$expandNodeNames(base.names)
    for(node in cont.nodes){
        customSpec$removeSamplers(node, print=FALSE)
        customSpec$addSampler(target = node,
                              type = "slice")
    }
    customSpec
}))

