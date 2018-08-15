## setwd("~/Dropbox/occupancy")
rm(list=ls())
setwd("analysis/singleSpp-multiSea")
source('src/initialize.R')
options(mc.cores=15)
set.seed(444)

## MCMC sampler options
cust.MCMCs <- c('nimble', 'jags_like_nimble', 'jags', 'RW_block', 'AFSS_block')
MCMC.defs <- c('nimble', MCMCdefs.slice,  'jags', MCMCdefs.RW.block,
               MCMCdefs.AFSS.block)
names(MCMC.defs) <- cust.MCMCs

## TRUE for model integrating over latent states
filter.opts <- c(FALSE, TRUE)
## true for model including hyper paramters for year effects on phi,
## gamma and p
hyper.param.opts <- c(TRUE, FALSE)

for(h in hyper.param.opts){
    for(ff in filter.opts) {
        hyper.param <- h
        ffilter <- ff
        if(ffilter){
            these.MCMCs <- cust.MCMCs[1:2]
        } else{
            these.MCMCs <- cust.MCMCs
        }
        runAllModels(ffilter=ffilter,
                     hyper.param=hyper.param,
                     niter=niter,
                     burnin=burnin,
                     MCMCs=these.MCMCs,
                     MCMCdefs=MCMC.defs)
    }
}

