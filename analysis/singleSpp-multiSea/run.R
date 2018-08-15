## setwd("~/Dropbox/occupancy")
rm(list=ls())
setwd("analysis/singleSpp-multiSea")
source('src/initialize.R')
options(mc.cores=15)
set.seed(444)

## MCMC sampler options
MCMCs <- c('nimble', 'jags', 'jags_like_nimble', 'RW_block', 'AFSS_block')
MCMCdefs <- c('nimble', 'jags', MCMCdefs.slice,  MCMCdefs.RW.block,
              MCMCdefs.AFSS.block)
names(MCMCdefs) <- MCMCs

## TRUE for model integrating over latent states
filter.opts <- c(TRUE)
## true for model including hyper paramters for year effects on phi,
## gamma and p
hyper.param.opts <- c(FALSE)

for(h in hyper.param.opts){
    for(ff in filter.opts) {
        hyper.param <- h
        ffilter <- ff
        if(ffilter) MCMCs[1:3]
        runAllModels(ffilter=ffilter,
                     hyper.param=hyper.param,
                     niter=niter,
                     burnin=burnin)
    }
}

