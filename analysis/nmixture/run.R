## setwd("~/Dropbox/hierarchical")
rm(list=ls())
setwd("analysis/nmixture")
source('src/initialize.R')

## all sampler options
## cust.MCMCs <- c('nimble', 'jags', 'block_RW', 'jags_like_nimble','block_AFSS')
## MCMC.defs <- c('nimble', 'jags', MCMCdefs.RW.block, MCMCdefs.slice,
##                MCMCdefs.AFSS.block)

## perry, perhaps try first with just nimble
run.models <- TRUE
cust.MCMCs <- c('nimble')
MCMC.defs <- c('nimble')

names(MCMC.defs) <- cust.MCMCs

## FALSE for model integrating over latent states
latent.opts <- c(TRUE)
## true for model including random effects
hyper.param.opts <- c(FALSE)

if(run.models){
    for(h in hyper.param.opts){
        for(ff in latent.opts) {
            if(ff){ ## latent states and hyper param
                these.MCMCs <- cust.MCMCs
            }else{
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
}


effsHP <- getEffFUN("hyperparamTRUE", save.dir,  summary="efficiency",
                    make.plot=make.comp.plots)

effsNoHP <- getEffFUN("hyperparamFALSE", save.dir,
                      summary="efficiency",
                      make.plot=make.comp.plots)

pdf.f(plotEffNmixture, file.path(save.dir, "../figures/comparisons/nmixture.pdf"),
      height=6, width=7)
