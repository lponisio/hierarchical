## setwd('~/Dropbox/hierarchical')
rm(list=ls())
setwd('analysis/multiSpp-singleSea')
source('src/initialize.R')
cust.MCMCs <- c('nimble', 'block_RW',
                'block_AFSS', 'jags_like_nimble', 'jags')
MCMC.defs <- c('nimble', MCMCdefs.RW.block,
               MCMCdefs.AFSS.block, MCMCdefs.slice, 'jags')


## TRUE for model integrating over latent states
latent.opts <- c(TRUE, FALSE)
## true for model including hyper paramters for year effects on phi,
## gamma and p
hyper.param.opts <- c(TRUE, FALSE)

if(run.models){
for(h in hyper.param.opts){
    for(ff in latent.opts) {
        if(ff){ ## latent states and hyper param
            these.MCMCs <- cust.MCMCs
        }else{
            these.MCMCs <- cust.MCMCs[-5]
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
                      adj.xlab=2,
                      make.plot=make.comp.plots)

pdf.f(plotEffMSSS, file.path(save.dir, "../figures/comparisons/MSSS.pdf"),
      height=6, width=7)


