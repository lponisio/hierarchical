## setwd('~/Dropbox/occupancy')
rm(list=ls())
setwd('analysis/multiSpp-singleSea')
source('src/initialize.R')

## MCMC sampler options
cust.MCMCs <- c('nimble', 'RW_block1', 'RW_block2',
                'AFSS_block1', 'AFSS_block2', 'jags')
MCMC.defs <- c('nimble', MCMCdefs.RW.block1,
               MCMCdefs.RW.block2, MCMCdefs.AFSS.block1,
               MCMCdefs.AFSS.block2, 'jags')
names(MCMC.defs) <- cust.MCMCs

## TRUE for model integrating over latent states
latent.opts <- c(TRUE, FALSE)
## true for model including hyper paramters for year effects on phi,
## gamma and p
hyper.param.opts <- c(TRUE, FALSE)

for(h in hyper.param.opts){
    for(ff in latent.opts) {
        if(ff & h){ ## latent states and hyper param
            these.MCMCs <- cust.MCMCs
        } else if (ff & !h){ ## latent states, but yes hyper param
            these.MCMCs <- cust.MCMCs[c(1,6)]
        } else if(!ff & h){
            these.MCMCs <- cust.MCMCs[1:5]
        } else if(!ff & !h){
            these.MCMCs <- cust.MCMCs[1]
        }

        runAllModels(latent=ff,
                     hyper.param=h,
                     niter=niter,
                     burnin=burnin,
                     MCMCs=these.MCMCs,
                     MCMCdefs=MCMC.defs)
    }
}

