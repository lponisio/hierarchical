args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
    run.models <- FALSE
    make.comp.plots <- FALSE
    mcmc.scale <- 2e2
} else{
    run.models <- args[1]
    make.comp.plots <- args[2]
    mcmc.scale <- as.numeric(args[3])
}

library(nimble)
library(igraph)
library(parallel)
source("src/dynamicOcc.R")
source("src/setup.R")
source("../all/plotting.R")
source("../all/misc.R")
source("src/plotting.R")
source('src/models.R')
source('src/customSamplerSpec.R')

save.dir <-  "../../../hierarchical_saved/singleSpp-multiSea/saved"

## MCMC settings
burnin <- 1e1*mcmc.scale
niter <- (1e3)*mcmc.scale

set.seed(444)
psi1 <- 0.2
sigma.p <- 0.02
mu.phi <- rnorm(1)
sigma.phi <- 0.02
mu.gamma <- rnorm(1)
sigma.gamma <- 0.02


runAllMCMC <- function(i, input1, niter, burnin, latent,
                       hyper.param, MCMCdefs, mu.p){
    print(sprintf("hyperparam%s_latent%s_sampler%s_mup%s",
                  hyper.param,
                  latent, i, mu.p))

    if(i == 'nimble' | i == 'jags'){
        ss.ms.samples <- compareMCMCs(input1,
                                      MCMCs=i,
                                      niter=niter,
                                      burnin = burnin,
                                      summary=FALSE,
                                      check=FALSE)
    } else{
        ss.ms.samples <- compareMCMCs(input1,
                                      MCMCs=i,
                                      MCMCdefs = MCMCdefs[i],
                                      niter=niter,
                                      burnin = burnin,
                                      summary=FALSE,
                                      check=FALSE)

    }
    save(ss.ms.samples,
         file=file.path(save.dir,
                        sprintf("hyperparam%s_latent%s_sampler%s_mup%s.Rdata",
                                hyper.param,
                                latent, i, mu.p)))
}



runAllModels <- function(latent, hyper.param,
                         niter,
                         burnin, MCMCs,
                         MCMCdefs,
                         mu.p,
                         data){

    model.input <- prepModDataOcc(data, include.zs=latent)
    if(!hyper.param){
        to.drop <- c("sigma.phi", "sigma.gamma", "sigma.p", "p",
                     "phi", "gamma")
        model.input$inits[to.drop] <- NULL
    }else{
        ## with new inits priors logit trick
        ## to.drop <- c("mu.phi", "mu.gamma", "mu.p")
        ## model.input$inits[to.drop] <- NULL
    }
    ss.ms.occ <- makeModel(latent, hyper.param)
    input1 <- c(code=ss.ms.occ,
                model.input)
    lapply(MCMCs, runAllMCMC, input1, niter, burnin,  latent,
           hyper.param, MCMCdefs, mu.p)
}
