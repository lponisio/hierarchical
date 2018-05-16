## This file contains work by Perry started 5/16/18.
## The goal is to explore sampling strategies for this model.

source("src/dataGen.R")
set.seed(444)
data <- genDynamicOccData()
model.input <- prepModDataOcc(data)

library(nimble)


ss.ms.occ <- nimbleCode({
  ## Specify priors
  psi1 ~ dunif(0, 1)

  for(year in 1:(nyear-1)){
    phi[year] ~ dunif(0, 1)
    gamma[year] ~ dunif(0, 1)
    p[year] ~ dunif(0, 1)
  }
  p[nyear] ~ dunif(0, 1)

  ## Ecological submodel: Define state conditional on parameters
  for (site in 1:nsite){
    z[site,1] ~ dbern(psi1)
    for (year in 2:nyear){
        muZ[site,year]<- z[site,year-1]*phi[year-1] +
            (1-z[site,year-1])*gamma[year-1]
      z[site,year] ~ dbern(muZ[site,year])
    }
  }

  ## Observation model
  for (site in 1:nsite){
    for (rep in 1:nrep){
      for (year in 1:nyear){
        muy[site,rep,year] <- z[site,year]*p[year]
        y[site,rep,year] ~ dbern(muy[site,rep,year])
      }
    }
  }

})

input1 <- c(code=ss.ms.occ,
            model.input)


## *********************************************************************
## original: vanilla nimble and JAGS
## *********************************************************************
niter <- 10000
burnin <- 100
ss.ms.orig <- compareMCMCs(input1,
                           MCMCs=c('jags', 'nimble'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)


##

assignSlice <- function(MCMCconf, reorder = TRUE) {
    nodes <- MCMCconf$model$expandNodeNames(c('phi','gamma','p', 'psi1'))
    MCMCconf$removeSamplers(nodes)
    for(node in nodes) {
        MCMCconf$addSampler(target = node, type = "slice")
    }
    ## The next section was an experiment that doesn't seem to make
    ## much difference.
    if(reorder) {
        model <- MCMCconf$model
        samplers <- MCMCconf$getSamplers()
        targets <- unlist(lapply(samplers, `[[`, 'target'))
        origOrder <- seq_along(targets)
        newOrder <- integer()
        i_psi1 <- which(targets == "psi1")
        i_z1 <- which(targets %in% model$expandNodeNames("z[, 1]"))
        newOrder <- c(newOrder, i_psi1, i_z1)
        for(t in 1:9) {
            i_phit <- which(targets == paste0("phi[",t,"]"))
            i_gammat <- which(targets == paste0("gamma[",t,"]"))
            i_ztp1 <- which(targets %in% model$expandNodeNames(paste0("z[,",t+1,"]")))
            newOrder <- c(newOrder, i_phit, i_gammat, i_ztp1)
        }
        rest <- origOrder[!(origOrder %in% newOrder)]
        newOrder <- c(newOrder, rest)
        if(length(unique(newOrder)) != length(targets))
            stop("Problem reordering samplers.")
    }
    MCMCconf$setSamplers(newOrder)
    MCMCconf
}

source("AF_slice_faster.R")

assignAFSS <- function(MCMCconf) {
    phiNodes <- MCMCconf$model$expandNodeNames(c('phi'))
    gammaNodes <- MCMCconf$model$expandNodeNames(c('gamma'))
    MCMCconf$removeSamplers(phiNodes)
    MCMCconf$removeSamplers(gammaNodes)
    for(i in seq_along(phiNodes)) {
        MCMCconf$addSampler(target = c(phiNodes[i], gammaNodes[i]), type = 'AF_slice_faster')
    }
    MCMCconf
}

## 

MCMCdefs <- list(
    nim_slice = quote({
        assignSlice(configureMCMC(Rmodel))
    }),
    nim_AFSS = quote({
        assignAFSS(configureMCMC(Rmodel))
    })
)

ss.ms.slice <- compareMCMCs(input1,
                            MCMCs=c('nim_slice', 'nim_AFSS'),
                            MCMCdefs = MCMCdefs,
                            niter=niter,
                            burnin = burnin,
                            summary=FALSE,
                            check=FALSE)

ss.ms.slice[[1]]$efficiency
ss.ms.slice[[1]]$summary
ss.ms.orig[[1]]$efficiency
