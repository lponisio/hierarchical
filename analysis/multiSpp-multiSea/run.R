## setwd('~/Dropbox/occupancy')
rm(list=ls())
setwd('analysis/multiSpp-multiSea')
source("src/initialize.R")

## MCMC sampler options
cust.MCMCs <- c('nimble', 'jags')
MCMC.defs <- c('nimble', 'jags')
names(MCMC.defs) <- cust.MCMCs


## FALSE for model integrating over latent states
latent.opts <- c(FALSE)
## true for model including hyper paramters for year effects on phi,
## gamma and p
hyper.param.opts <- c(FALSE)

for(h in hyper.param.opts){
    for(l in latent.opts) {
        if(l & h){ ## latent states and hyper param
            these.MCMCs <- cust.MCMCs
        } else if (l & !h){ ## latent states but no hyper param
            these.MCMCs <- cust.MCMCs
        } else if(!l & h){ ## no latent states, but yes hyper param
            these.MCMCs <- cust.MCMCs[1]
        } else if(!l & !h){ ## no latent states, no hyper param
            these.MCMCs <- cust.MCMCs[1]
        }

        runAllModels(latent=l,
                     hyper.param=h,
                     niter=niter,
                     burnin=burnin,
                     MCMCs=these.MCMCs,
                     MCMCdefs=MCMC.defs)
    }
}
