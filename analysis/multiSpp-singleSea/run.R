## setwd('~/Dropbox/occupancy')
rm(list=ls())
setwd('analysis/multiSpp-singleSea')

source('src/initialize.R')
latent <- c(FALSE)
hyper.param <- c(TRUE)


## MCMC sampler options
cust.MCMCs <- c('nimble', 'jags')
MCMC.defs <- c('nimble',  'jags')
names(MCMC.defs) <- cust.MCMCs

## TRUE for model integrating over latent states
latent.opts <- c(TRUE, FALSE)
## true for model including hyper paramters for year effects on phi,
## gamma and p
hyper.param.opts <- c(TRUE, FALSE)

for(h in hyper.param.opts){
    for(ff in latent.opts) {
        if(ff){ ## filtering
            these.MCMCs <- cust.MCMCs[1]
        } else if(!h){ ## no filtering or hyper param
            these.MCMCs <- cust.MCMCs
        } else{ ## no filtering, but yes hyper param
            these.MCMCs <- cust.MCMCs
        }

        runAllModels(latent=ff,
                     hyper.param=h,
                     niter=niter,
                     burnin=burnin,
                     MCMCs=these.MCMCs,
                     MCMCdefs=MCMC.defs)
    }
}

