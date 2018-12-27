## setwd('~/Dropbox/occupancy')
rm(list=ls())
setwd('analysis/multiSpp-singleSea')
source('src/initialize.R')
options(mc.cores=15)

cust.MCMCs <- c('nimble', 'block_RW',
                'block_AFSS', 'jags_like_nimble', 'jags')
MCMC.defs <- c('nimble', MCMCdefs.RW.block,
               MCMCdefs.AFSS.block, MCMCdefs.slice, 'jags')


## TRUE for model integrating over latent states
latent.opts <- c(TRUE, FALSE)
## true for model including hyper paramters for year effects on phi,
## gamma and p
hyper.param.opts <- c(TRUE, FALSE)

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

effsHP <- getEffFUN("hyperparamTRUE", save.dir,  summary="efficiency")
effsNoHP <- getEffFUN("hyperparamFALSE", save.dir,  summary="efficiency")


pdf.f(plotEffMSSS, file.path(save.dir, "../../../figures/MSSS.pdf"),
      height=6, width=7)


