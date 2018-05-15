## setwd('~/Dropbox/occupancy')
rm(list=ls())
library(nimble)
library(coda)
setwd('analysis/nmixture')
source('../all/plotting.R')
source('src/initialize.R')

## original model jags and nimble
load(file=file.path(save.dir, "orig.Rdata"))

## vectorized, likelihood for latent state, derived quantity
## calculation
load(file=file.path(save.dir, "grouped.Rdata"))

nmixture.orig[[1]] <- rename_MCMC_comparison_method(c('nimble', 'jags'),
                                                 c('NIMBLE',
                                                   'JAGS'),
                                              comparison=nmixture.orig[[1]])

nmixture.grouped[[1]] <- rename_MCMC_comparison_method(c('nimble',
                                                 'AFSS', 'RWB5'),
                                                 c('NIMBLE-grouped',
                                                   'NIMBLE-AFSS',
                                                   'NIMBLE-blocking'),
                                              comparison=nmixture.grouped[[1]])

## compare mcmcs
nmixture.occ.all <- combine_MCMC_comparison_results(nmixture.orig[[1]],
                                                 nmixture.grouped[[1]],
                                                 name = "nmixture")
save(nmixture.occ.all, file=file.path(save.dir, "combined.Rdata"))


make_MCMC_comparison_pages(nmixture.occ.all,
                           dir=file.path(save.dir, "../figures/comparisons"))


## look at samples

checkChains(nmixture.occ.all$nmixture$samples,
            f.path = file.path(save.dir,
                               "../figures/chains/%s.pdf"))

## ****************************************
## custom figs
## ****************************************

by.param <- apply(nmixture.occ.all$nmixture$samples, c(1,2), effectiveSize)/
  nmixture.occ.all$nmixture$timing
by.config <- nmixture.occ.all$nmixture$efficiency

source('../all/plotting.R')
plotEffSize(by.config, by.param, f.path= file.path(save.dir,
              "../figures/comparisons/%s%s.pdf"), "MultiSpp-SingleSea",
            at=0.3, adj1=0.03, adj2=0.1)

(nmixture.occ.all$nmixture$efficiency$mean["NIMBLE-latent"] -
 nmixture.occ.all$nmixture$efficiency$mean["JAGS-latent"])/
    max(nmixture.occ.all$nmixture$efficiency$mean["NIMBLE-latent"])

(nmixture.occ.all$nmixture$efficiency$mean["NIMBLE-filter"] -
 nmixture.occ.all$nmixture$efficiency$mean["JAGS-latent"])/
    max(nmixture.occ.all$nmixture$efficiency$mean["NIMBLE-filter"])

