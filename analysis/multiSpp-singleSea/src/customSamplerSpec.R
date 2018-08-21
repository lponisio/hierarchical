## *********************************************************************
##  block sampler for species random effects for each species
## *********************************************************************

MCMCdefs.RW.block1 <- list('RW_block1' = quote({
  customSpec <- configureMCMC(Rmodel)
  ## find node names of each species for random effects
  base.names <- c("a1", "a2", "a3", "a4", "b1", "b2", "u.cato",
                  "u.fcw", "v.cato", "v.fcw" )
  exp.names.list <- list()
  for(bn in base.names){
    exp.names.list[[bn]] <- Rmodel$expandNodeNames(bn)
  }
  for(i in 1:length(exp.names.list[[1]])){
    blocknames <- unlist(lapply(exp.names.list, function(x) x[i]))
    customSpec$removeSamplers(blocknames, print=FALSE)
    customSpec$addSampler(target = blocknames,
                          type = "RW_block")
  }
  customSpec
}))



MCMCdefs.AFSS.block1 <- list('AFSS_block1' = quote({
  customSpec <- configureMCMC(Rmodel)
  ## find node names of each species for random effects
  base.names <- c("a1", "a2", "a3", "a4", "b1", "b2", "u.cato",
                  "u.fcw", "v.cato", "v.fcw" )
  exp.names.list <- list()
  for(bn in base.names){
    exp.names.list[[bn]] <- Rmodel$expandNodeNames(bn)
  }
  for(i in 1:length(exp.names.list[[1]])){
    blocknames <- unlist(lapply(exp.names.list, function(x) x[i]))
    customSpec$removeSamplers(blocknames, print=FALSE)
    customSpec$addSampler(target = blocknames,
                          type = "AF_slice")
  }
  customSpec
}))


## *********************************************************************
## ##  block sampler for species random effects for each
## "type"
## *********************************************************************

## remove the samples, add block samplers
MCMCdefs.RW.block2 <- list('RW_block2' = quote({
  customSpec <- configureMCMC(Rmodel)
  ## find node names for random effects
  sp.parms.a <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^a",
                                      Rmodel$getNodeNames(includeData = FALSE))]
  sp.parms.b <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^b",
                                      Rmodel$getNodeNames(includeData = FALSE))]
  sp.parms.u <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^u",
                                      Rmodel$getNodeNames(includeData = FALSE))]
  sp.parms.v <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^v",
                                      Rmodel$getNodeNames(includeData = FALSE))]
  customSpec$removeSamplers(c(sp.parms.a, sp.parms.b, sp.parms.u,
                              sp.parms.v),
                            print=FALSE)

  customSpec$addSampler(target = sp.parms.a,
                        type = "RW_block")
  customSpec$addSampler(target = sp.parms.b,
                        type = "RW_block")
  customSpec$addSampler(target = sp.parms.u,
                        type = "RW_block")
  customSpec$addSampler(target = sp.parms.v,
                        type = "RW_block")
  customSpec
}))



MCMCdefs.AFSS.block2 <- list('AFSS_block2' = quote({
  customSpec <- configureMCMC(Rmodel)
  ## find node names for random effects
  sp.parms.a <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^a",
                                      Rmodel$getNodeNames(includeData = FALSE))]
  sp.parms.b <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^b",
                                      Rmodel$getNodeNames(includeData = FALSE))]
  sp.parms.u <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^u",
                                      Rmodel$getNodeNames(includeData = FALSE))]
  sp.parms.v <- Rmodel$getNodeNames(includeData = FALSE)[grepl("^v",
                                      Rmodel$getNodeNames(includeData = FALSE))]
  customSpec$removeSamplers(c(sp.parms.a, sp.parms.b, sp.parms.u,
                              sp.parms.v),
                            print=FALSE)

  customSpec$addSampler(target = sp.parms.a,
                        type = "AF_slice")
  customSpec$addSampler(target = sp.parms.b,
                        type = "AF_slice")
  customSpec$addSampler(target = sp.parms.u,
                        type = "AF_slice")
  customSpec$addSampler(target = sp.parms.v,
                        type = "AF_slice")
  customSpec
}))
