## setwd("~/Dropbox/hierarchical")
rm(list=ls())
setwd("analysis/nmixture")
source('src/initialize.R')

## all sampler options
## cust.MCMCs <- c('nimble', 'jags', 'block_RW', 'jags_like_nimble','block_AFSS')
## MCMC.defs <- c('nimble', 'jags', MCMCdefs.RW.block, MCMCdefs.slice,
##                MCMCdefs.AFSS.block)

cust.MCMCs <- c('nimble')
MCMC.defs <- c('nimble')

make.comp.plots <- TRUE
run.models <- TRUE

names(MCMC.defs) <- cust.MCMCs

## FALSE for model integrating over latent states
latent.opts <- c(FALSE)
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

## ## Perry's code to inspect a single result
## list.files(save.dir)
## iFile <- 3 ## choose the desired file
## load(file.path(save.dir, list.files(save.dir)[iFile]))

effsHP <- getEffFUN("hyperparamTRUE", save.dir,  summary="efficiency",
                    make.plot=make.comp.plots)

effsNoHP <- getEffFUN("hyperparamFALSE", save.dir,
                      summary="efficiency",
                      make.plot=make.comp.plots)

pdf.f(plotEffNmixture, file.path(save.dir,
                                 "../figures/comparisons/nmixture.pdf"),
      height=6, width=7)


mapply(function(a, b)
                   all(a == b),
                   a=latent$data,
                   b=nolatent$data,
                   SIMPLIFY=FALSE)



mapply(function(a, b)
                   all(a == b),
                   a=latent$inits[-1],
                   b=nolatent$inits,
                   SIMPLIFY=FALSE)



mapply(function(a, b)
                   all(a == b),
                   a=latent$constants,
                   b=nolatent$constants,
                   SIMPLIFY=FALSE)
