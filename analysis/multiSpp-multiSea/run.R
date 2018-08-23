## setwd('~/Dropbox/occupancy')
rm(list=ls())
setwd('analysis/multiSpp-multiSea')
source("src/initialize.R")

## ## MCMC sampler options
## cust.MCMCs <- c('nimble', 'jags', 'RW_block', 'AFSS_block')
## MCMC.defs <- c('nimble',  'jags', MCMCdefs.RW.block,
##                MCMCdefs.AFSS.block)
## names(MCMC.defs) <- cust.MCMCs



## MCMC sampler options
cust.MCMCs <- c('jags', 'nimble')
MCMC.defs <- c('jags', 'nimble')
names(MCMC.defs) <- cust.MCMCs


## FALSE for model integrating over latent states
latent.opts <- c(TRUE, FALSE)
## true for model including hyper paramters for year effects on phi,
## gamma and p
hyper.param.opts <- c(TRUE, FALSE)

## for(h in hyper.param.opts){
##     for(l in latent.opts) {
##         if(l){ ## latent states
##             these.MCMCs <- cust.MCMCs
##         } else{ ## no latent states, but yes hyper param
##             these.MCMCs <- cust.MCMCs[c(1)]
##         }

##         runAllModels(latent=l,
##                      hyper.param=h,
##                      niter=niter,
##                      burnin=burnin,
##                      MCMCs=these.MCMCs,
##                      MCMCdefs=MCMC.defs)
##     }
## }


effsHP <- getEffFUN("hyperparamTRUE", save.dir,  summary="efficiency")
effsNoHP <- getEffFUN("hyperparamFALSE", save.dir,  summary="efficiency")

meanHP <- getEffFUN("hyperparamTRUE", save.dir,  summary="mean")
meanNoHP <- getEffFUN("hyperparamFALSE", save.dir,  summary="mean")


pdf.f(plotEffMSMS, file.path(save.dir, "../../../figures/MSMS.pdf"),
      height=6, width=7)


