rm(list=ls())
library(nimble)
library(coda)
setwd('~/Dropbox')
setwd('occupancy/analysis/multiSpp-singleSea')
source('../all/plotting.R')
source('src/initialize.R')

thinChains <- function(chains, nthin){
    apply(chains, c(1,2), function(x){
        x[!seq(nthin, length(x))] <- NULL
    })
}

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
                                                 ## ms.ss.crosslevel[[1]],
                                                 ## ms.ss.subsamp[[1]],
                                                 name = "ms.ss")
save(ms.ss.occ.all, file=file.path(save.dir, "combined.Rdata"))

make_MCMC_comparison_pages(ms.ss.occ.all,
                           dir=file.path(save.dir, "../figures/comparisons"))


## look at samples

checkChains(ms.ss.occ.all$ms.ss$samples,
            f.path = file.path(save.dir,

## paramter groups
params <- dimnames(ms.ss.occ.all$ms.ss$summary)[[3]]
groups <- list()
i <- 1
while(length(params) > 0){
    id <- agrep(params[1], params, max.distance = 0.2)
    groups[[i]] <- params[id]
    params <- params[-id]
    i <- i + 1
}

f <- function(){
    layout(matrix(1:3, ncol=1))
    par(oma=c(3,1,1,1))
    cols <- rainbow(dim(ms.ss.occ.all$ms.ss$summary)[1])
    for(group in groups){
        if(length(group) > 1){
            for(i in 1:dim(ms.ss.occ.all$ms.ss$summary)[1]){
                this.samp <- ms.ss.occ.all$ms.ss$summary[i,,group]
                xs <- jitter(1:dim(this.samp)[2])
                if(i == 1){
                    plot(x=xs, y=this.samp["mean",], pch=16,
                         col=cols[i],
                         ylim=range(c(ms.ss.occ.all$ms.ss$summary[,'CI95_upp',
                                                                  group],
                                      ms.ss.occ.all$ms.ss$summary[,'CI95_low',
                                                                  group])),
                         xaxt="n",
                         ylab="Estimate",
                         xlab="")

                    axis(1, at=1:dim(this.samp)[2],
                         labels=FALSE)
                    text(x=1:dim(this.samp)[2],
                         y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
                         labels=group, srt=45, adj=1, xpd=TRUE)
                } else{
                    points(x=xs,
                           y=this.samp["mean",], pch=16, col=cols[i])
                }
                arrows(y1=this.samp['CI95_upp',],
                       y0=this.samp['CI95_low',],
                       x0=xs,
                       code=0, angle=90, length=0.02, lwd=1, col=cols[i])
            }
        }

        legend("topright", legend=dimnames(ms.ss.occ.all$ms.ss$summary)[[1]],
               pch=16, col=cols)
    }
}

pdf.f(f, file=file.path(save.dir, "../figures/allparams.pdf"),
      height=10, width=8.5 )

## ****************************************
## custom figs
## ****************************************

by.param <- apply(ms.ss.occ.all$ms.ss$samples, c(1,2), effectiveSize)/
  ms.ss.occ.all$ms.ss$timing
by.config <- ms.ss.occ.all$ms.ss$efficiency

source('../all/plotting.R')
plotEffSize(by.config, by.param, f.path= file.path(save.dir,
              "../figures/comparisons/%s%s.pdf"), "MultiSpp-SingleSea",
            at=0.4, adj1=0.03, adj2=0.1)

(ms.ss.occ.all$ms.ss$efficiency$mean["NIMBLE-latent"] -
 ms.ss.occ.all$ms.ss$efficiency$mean["JAGS-latent"])/
    max(ms.ss.occ.all$ms.ss$efficiency$mean["NIMBLE-latent"])

(ms.ss.occ.all$ms.ss$efficiency$mean["NIMBLE-filter"] -
 ms.ss.occ.all$ms.ss$efficiency$mean["JAGS-latent"])/
    max(ms.ss.occ.all$ms.ss$efficiency$mean["NIMBLE-filter"])

