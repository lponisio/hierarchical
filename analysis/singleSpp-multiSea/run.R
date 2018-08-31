## setwd("~/Dropbox/occupancy")
rm(list=ls())
setwd("analysis/singleSpp-multiSea")
source('src/initialize.R')
options(mc.cores=15)


## MCMC sampler options
cust.MCMCs <- c('nimble', 'jags', 'RW_block', 'jags_like_nimble','AFSS_block')
MCMC.defs <- c('nimble', 'jags', MCMCdefs.RW.block, MCMCdefs.slice,
               MCMCdefs.AFSS.block)

names(MCMC.defs) <- cust.MCMCs

## FALSE for model integrating over latent states
latent.opts <- c(TRUE, FALSE)
## true for model including hyper paramters for year effects on phi,
## gamma and p
hyper.param.opts <- c(TRUE, FALSE)
## easy or difficult to detect
mus.p <- c(1.5, 0.2)

for(mu.p in mus.p){
    for(h in hyper.param.opts){
        for(l in latent.opts) {
            if(!l){ ## filtering, no latent states
                these.MCMCs <- cust.MCMCs[-2]
            } else if(!h){ ## no filtering or hyper param
                these.MCMCs <- cust.MCMCs
            } else{ ## no filtering, but yes hyper param
                these.MCMCs <- cust.MCMCs
            }

            runAllModels(latent=l,
                         hyper.param=h,
                         niter=niter,
                         burnin=burnin,
                         MCMCs=these.MCMCs,
                         MCMCdefs=MCMC.defs,
                         mu.p=mu.p,
                         psi1=psi1,
                         sigma.p=sigma.p,
                         mu.phi=mu.phi,
                         sigma.phi=sigma.phi,
                         mu.gamma=mu.gamma,
                         sigma.gamma=sigma.gamma)
        }
    }
}


effsHP <- getEffFUN("hyperparamTRUE", save.dir,  summary="efficiency")
effsNoHP <- getEffFUN("hyperparamFALSE", save.dir,  summary="efficiency")

meanHP <- getEffFUN("hyperparamTRUE", save.dir,  summary="mean")
mu.p <- NA
meanHP <- rbind(meanHP, sim=lapply(colnames(meanHP), get))

meanNoHP <- getEffFUN("hyperparamFALSE", save.dir,  summary="mean")
meanNoHP <- rbind(meanNoHP, sim=lapply(colnames(meanNoHP), get))

pdf.f(plotEffSSMS, file.path(save.dir, "../../../figures/SSMS.pdf"),
      height=6, width=7)

