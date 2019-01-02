## setwd("~/Dropbox/occupancy")
rm(list=ls())
setwd("analysis/singleSpp-multiSea")
source('src/initialize.R')
## MCMC sampler options
cust.MCMCs <- c('nimble', 'jags', 'block_RW', 'jags_like_nimble','block_RW')
MCMC.defs <- c('nimble', 'jags', MCMCdefs.RW.block, MCMCdefs.slice,
               MCMCdefs.AFSS.block)
names(MCMC.defs) <- cust.MCMCs


## FALSE for model integrating over latent states
latent.opts <- c(TRUE, FALSE)
## true for model including hyper paramters for year effects on phi,
## gamma and p
hyper.param.opts <- c(TRUE, FALSE)
## easy or difficult to detect
mus.p <- c(1, -1)

if(run.models){
    for(mu.p in mus.p){
        data <- genDynamicOccData(mu.p=mu.p,
                                  psi1=psi1,
                                  sigma.p=sigma.p,
                                  mu.phi=mu.phi,
                                  sigma.phi=sigma.phi,
                                  mu.gamma=mu.gamma,
                                  sigma.gamma=sigma.gamma)
        for(h in hyper.param.opts){
            for(l in latent.opts) {
                if(!l){ ## filtering, no latent states
                    these.MCMCs <- cust.MCMCs[-2]
                } else{ ## no filtering, either year or no hyper param
                    these.MCMCs <- cust.MCMCs
                }
                runAllModels(latent=l,
                             hyper.param=h,
                             niter=niter,
                             burnin=burnin,
                             MCMCs=these.MCMCs,
                             MCMCdefs=MCMC.defs,
                             mu.p=mu.p,
                             data=data)
            }
        }
    }
}


effsHP <- getEffFUN("hyperparamTRUE", save.dir,  summary="efficiency",
                    make.plot=make.comp.plots)

effsNoHP <- getEffFUN("hyperparamFALSE", save.dir,
                      summary="efficiency",
                      make.plot=make.comp.plots)

pdf.f(plotEffSSMS, file.path(save.dir, "../../../figures/SSMS.pdf"),
      height=6, width=7)


