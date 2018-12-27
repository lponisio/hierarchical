## setwd('~/Dropbox/occupancy')
rm(list=ls())
setwd('analysis/multiSpp-multiSea')
source("src/initialize.R")

## MCMC sampler options
cust.MCMCs <- c('nimble', 'jags','block_RW',
                'block_AFSS', 'jags_like_nimble')
MCMC.defs <- c('nimble', 'jags', MCMCdefs.RW.block,
               MCMCdefs.AFSS.block, MCMCdefs.slice )
names(MCMC.defs) <- cust.MCMCs



## FALSE for model integrating over latent states
latent.opts <- c(TRUE, FALSE)
## true for model including hyper paramters for year effects on phi,
## gamma and p
hyper.param.opts <- c(FALSE)
for(h in hyper.param.opts){
    for(l in latent.opts) {
        if(l){ ## latent states
            these.MCMCs <- cust.MCMCs
        } else{ ## no latent states, but yes hyper param, cannot use JAGS
            these.MCMCs <- cust.MCMCs[-2]
        }

        runAllModels(latent=l,
                     hyper.param=h,
                     niter=niter,
                     burnin=burnin,
                     MCMCs=these.MCMCs,
                     MCMCdefs=MCMC.defs)
    }
}

effsHP <- getEffFUN("hyperparamTRUE", save.dir,  summary="efficiency")
effsNoHP <- getEffFUN("hyperparamFALSE", save.dir,  summary="efficiency")

pdf.f(plotEffMSMS, file.path(save.dir, "../../../figures/MSMS.pdf"),
      height=6, width=7)






