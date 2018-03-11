## setwd('~/Dropbox/occupancy')
rm(list=ls())
library(nimble)
library(coda)
setwd('analysis/multiSpp-singleSea')
source('../all/plotting.R')
source('src/initialize.R')
source('src/plotPosteriorGroups.R')

## original model jags and nimble
load(file=file.path(save.dir, "orig.Rdata"))

## vectorized, likelihood for latent state, derived quantity
## calculation
load(file=file.path(save.dir, "filter.Rdata"))

## cross level sampler
load(file=file.path(save.dir, "crosslevel.Rdata"))

## subsampling
load(file=file.path(save.dir, "subsamp.Rdata"))

ms.ss.orig[[1]] <- rename_MCMC_comparison_method(c('nimble', 'jags'),
                                                 c('NIMBLE-latent',
                                                   'JAGS-latent'),
                                              comparison=ms.ss.orig[[1]])

ms.ss.filter[[1]] <- rename_MCMC_comparison_method(c('nimble'),
                                                 c('NIMBLE-filter'),
                                              comparison=ms.ss.filter[[1]])

ms.ss.crosslevel[[1]] <- rename_MCMC_comparison_method('nimbleCrosslevel',
                                                       'NIMBLE-cross-level',
                                              comparison=ms.ss.crosslevel[[1]])

ms.ss.subsamp[[1]] <- rename_MCMC_comparison_method('nimbleSubsamp',
                                                    'NIMBLE-subsample',
                                              comparison=ms.ss.subsamp[[1]])

## compare mcmcs
ms.ss.occ.all <- combine_MCMC_comparison_results(ms.ss.orig[[1]],
                                                 ms.ss.filter[[1]],
                                                 ms.ss.crosslevel[[1]],
                                                 ms.ss.subsamp[[1]],
                                                 name = "ms.ss")
save(ms.ss.occ.all, file=file.path(save.dir, "combined.Rdata"))

make_MCMC_comparison_pages(ms.ss.occ.all,
                           dir=file.path(save.dir, "../figures/comparisons"))


## look at samples

checkChains(ms.ss.occ.all$ms.ss$samples,
            f.path = file.path(save.dir,
                               "../figures/chains/%s.pdf"))

## ****************************************
## custom figs
## ****************************************

by.param <- apply(ms.ss.occ.all$ms.ss$samples, c(1,2), effectiveSize)/
  ms.ss.occ.all$ms.ss$timing
by.config <- ms.ss.occ.all$ms.ss$efficiency

source('../all/plotting.R')
plotEffSize(by.config, by.param, f.path= file.path(save.dir,
              "../figures/comparisons/%s%s.pdf"), "MultiSpp-SingleSea",
            at=0.3, adj1=0.03, adj2=0.1)

(ms.ss.occ.all$ms.ss$efficiency$mean["NIMBLE-latent"] -
 ms.ss.occ.all$ms.ss$efficiency$mean["JAGS-latent"])/
    max(ms.ss.occ.all$ms.ss$efficiency$mean["NIMBLE-latent"])

(ms.ss.occ.all$ms.ss$efficiency$mean["NIMBLE-filter"] -
 ms.ss.occ.all$ms.ss$efficiency$mean["JAGS-latent"])/
    max(ms.ss.occ.all$ms.ss$efficiency$mean["NIMBLE-filter"])


## ****************************************
## param groups
## ****************************************
params <- dimnames(ms.ss.occ.all$ms.ss$summary)[[3]]

groups <- list()
i <- 1
while(length(params) > 0){
    id <- agrep(params[1], params, max.distance = 0.2)
    groups[[i]] <- params[id]
    params <- params[-id]
    i <- i + 1
}

pdf.f(plotPosteriorGroups, file=file.path(save.dir, "../figures/allparams.pdf"),
      height=10, width=8.5 )

