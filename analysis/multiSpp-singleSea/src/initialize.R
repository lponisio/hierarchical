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
library(reshape)
source("../all/plotting.R")
source("src/plotting.R")
source("src/setup.R")
source("src/multispeciesOcc.R")
source("src/models.R")
source("src/customSamplerSpec.R")
source("../all/misc.R")

save.dir <-  "../../../hierarchical_saved/multiSpp-singleSea/saved"

survey.data <- read.csv("data/occupancy_data.csv")
species.groups <- read.csv("data/species_groups.csv")
survey.dates <- read.csv("data/survey_dates.csv")
habitat <- read.csv("data/habitat.csv")

## mcmc settings
burnin <- 1e2*mcmc.scale
niter <- (1e3)*mcmc.scale

runAllMCMC <- function(i, input1, niter, burnin, latent,
                       hyper.param, MCMCdefs){
    print(sprintf("hyperparam%s_latent%s_sampler%s",
                  hyper.param,
                  latent, i))
    if(i == 'nimble' | i == 'jags'){
        ms.ss.samples <- compareMCMCs(input1,
                                      MCMCs=i,
                                      niter=niter,
                                      burnin = burnin,
                                      summary=FALSE,
                                      check=FALSE)
    } else{
        ms.ss.samples <- compareMCMCs(input1,
                                      MCMCs=i,
                                      MCMCdefs = MCMCdefs[i],
                                      niter=niter,
                                      burnin = burnin,
                                      summary=FALSE,
                                      check=FALSE)

    }
    save(ms.ss.samples,
         file=file.path(save.dir,
                        sprintf("hyperparam%s_latent%s_sampler%s.Rdata",
                                hyper.param,
                                latent, i)))
}



runAllModels <- function(latent, hyper.param, niter, burnin,
                         MCMCs, MCMCdefs){
    model.input <- prepMutiSpData(survey.data,
                                  survey.dates,
                                  species.groups,
                                  habitat,
                                  n.zeroes =0, ## don't augment data
                                  remove.zs=!latent,
                                  hyper.param=hyper.param)
    ms.ss.occ <- makeModel(latent, hyper.param)
    input1 <- c(code=ms.ss.occ,
                model.input)
    lapply(MCMCs, runAllMCMC, input1, niter, burnin,  latent,
           hyper.param, MCMCdefs)
}
