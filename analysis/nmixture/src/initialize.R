args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
    run.models <- FALSE
    make.comp.plots <- FALSE
} else{
    run.models <- args[1]
    make.comp.plots <- args[2]
}

library(nimble)
library(igraph)
source("src/dNmixture.R")
source("src/models.R")
source("src/customSamplerSpec.R")
source("src/setup.R")
source("../all/plotting.R")
source("src/plotting.R")

save.dir <-  "../../../hierarchical_saved/nmixture/saved"

## mcmc settings
scale <- 1
burnin <- 1e2*scale
niter <- (1e3)*scale


runAllMCMC <- function(i, input1, niter, burnin, latent,
                       hyper.param, MCMCdefs){
    print(sprintf("hyperparam%s_latent%s_sampler%s",
                  hyper.param,
                  latent, i))

    if(i == 'nimble' | i == 'jags'){
        nmixture.samples <- compareMCMCs(input1,
                                         MCMCs=i,
                                         niter=niter,
                                         burnin = burnin,
                                         summary=FALSE,
                                         check=FALSE)
    } else{
        nmixture.samples <- compareMCMCs(input1,
                                         MCMCs=i,
                                         MCMCdefs = MCMCdefs[i],
                                         niter=niter,
                                         burnin = burnin,
                                         summary=FALSE,
                                         check=FALSE)

    }
    save(nmixture.samples,
         file=file.path(save.dir,
                        sprintf("hyperparam%s_latent%s_sampler%s.Rdata",
                                hyper.param,
                                latent, i)))
}



runAllModels <- function(latent, hyper.param, niter, burnin,
                         MCMCs, MCMCdefs){
    model.input <- prepNmixtureData(latent=latent,
                                    hyper.param=hyper.param)

    nmixture <- makeModel(latent, hyper.param)
    input1 <- c(code=nmixture,
                model.input)
    lapply(MCMCs, runAllMCMC, input1, niter, burnin,  latent,
           hyper.param, MCMCdefs)
}
