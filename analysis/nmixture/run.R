## setwd("~/Dropbox/occupancy")
rm(list=ls())
setwd("analysis/nmixture")
source('src/initialize.R')

cust.MCMCs <- c('nimble', 'jags', 'RW_block', 'jags_like_nimble','AFSS_block')
MCMC.defs <- c('nimble', 'jags', MCMCdefs.RW.block, MCMCdefs.slice,
               MCMCdefs.AFSS.block)
names(MCMC.defs) <- cust.MCMCs

## TRUE for model integrating over latent states
latent.opts <- c(FALSE)
## true for model including random effects
hyper.param.opts <- c(FALSE)

for(h in hyper.param.opts){
    for(ff in latent.opts) {
        if(ff){ ## latent states and hyper param
            these.MCMCs <- cust.MCMCs
        }else{
            these.MCMCs <- cust.MCMCs[-2]
        }
        runAllModels(latent=ff,
                     hyper.param=h,
                     niter=niter,
                     burnin=burnin,
                     MCMCs=these.MCMCs,
                     MCMCdefs=MCMC.defs)
    }
}
