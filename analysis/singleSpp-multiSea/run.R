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
filter.opts <- c(TRUE, FALSE)
## true for model including hyper paramters for year effects on phi,
## gamma and p
hyper.param.opts <- c(TRUE, FALSE)
## easy or difficult to detect
mus.p <- c(1.5, 0.2)

for(mu.p in mus.p){
    for(h in hyper.param.opts){
        for(ff in filter.opts) {
            if(ff){ ## filtering
                these.MCMCs <- cust.MCMCs[1:2]
            } else if(!h){ ## no filtering or hyper param
                these.MCMCs <- cust.MCMCs[1:3]
            } else{ ## no filtering, but yes hyper param
                these.MCMCs <- cust.MCMCs
            }

            runAllModels(ffilter=ff,
                         hyper.param=h,
                         niter=niter,
                         burnin=burnin,
                         MCMCs=these.MCMCs,
                         MCMCdefs=MCMC.defs,
                         mu.p=mu.p)
        }
    }
}
