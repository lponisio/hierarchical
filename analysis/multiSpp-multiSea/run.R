## setwd('~/Dropbox/occupancy')
rm(list=ls())
setwd('analysis/multiSpp-multiSea')
source("src/initialize.R")

## MCMC sampler options
cust.MCMCs <- c('nimble', 'jags_like_nimble', 'jags', 'RW_block', 'AFSS_block')
MCMC.defs <- c('nimble', MCMCdefs.slice,  'jags', MCMCdefs.RW.block,
               MCMCdefs.AFSS.block)
names(MCMC.defs) <- cust.MCMCs

## FALSE for model integrating over latent states
latent.opts <- c(TRUE, FALSE)
## true for model including hyper paramters for year effects on phi,
## gamma and p
hyper.param.opts <- c(TRUE, FALSE)

for(h in hyper.param.opts){
    for(l in latent.opts) {
        if(l){ ## latent states
            these.MCMCs <- cust.MCMCs
        } else{ ## no latent states, but yes hyper param
            these.MCMCs <- cust.MCMCs[c(1,2,4,5)]
        }

        runAllModels(latent=l,
                     hyper.param=h,
                     niter=niter,
                     burnin=burnin,
                     MCMCs=these.MCMCs,
                     MCMCdefs=MCMC.defs)
    }
}
