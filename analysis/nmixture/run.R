## setwd("~/Dropbox/occupancy")
rm(list=ls())
setwd("analysis/nmixture")
source('src/initialize.R')

## cust.MCMCs <- c('nimble', 'jags', 'block_RW', 'jags_like_nimble','block_AFSS')
## MCMC.defs <- c('nimble', 'jags', MCMCdefs.RW.block, MCMCdefs.slice,
##                MCMCdefs.AFSS.block)


## cust.MCMCs <- c('nimble', 'jags')
## MCMC.defs <- c('nimble', 'jags')

## names(MCMC.defs) <- cust.MCMCs

## ## FALSE for model integrating over latent states
## latent.opts <- c(FALSE)
## ## true for model including random effects
## hyper.param.opts <- c(TRUE, FALSE)

## for(h in hyper.param.opts){
##     for(ff in latent.opts) {
##         if(ff){ ## latent states and hyper param
##             these.MCMCs <- cust.MCMCs
##         }else{
##             these.MCMCs <- cust.MCMCs[-2]
##         }
##         runAllModels(latent=ff,
##                      hyper.param=h,
##                      niter=niter,
##                      burnin=burnin,
##                      MCMCs=these.MCMCs,
##                      MCMCdefs=MCMC.defs)
##     }
## }


effsHP <- getEffFUN("hyperparamTRUE", save.dir,  summary="efficiency")

effsNoHP <- getEffFUN("hyperparamFALSE", save.dir,  summary="efficiency")

pdf.f(plotEffNmixture, file.path(save.dir, "../../../figures/nmixture.pdf"),
      height=6, width=7)
