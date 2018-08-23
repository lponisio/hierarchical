library(devtools)
## install_github("nimble-dev/nimble",
##                ref = "devel",
##                subdir = "packages/nimble")

library(nimble)
library(parallel)
library(igraph)
library(reshape)
source("../all/plotting.R")
source("src/plotting.R")
source("src/setup.R")
source("src/multispeciesOcc.R")
source("src/models.R")
source("src/customSamplerSpec.R")

dir.create(file.path("../../../occupancy_saved/saved/multiSpp-singleSea/saved"),
           showWarnings = FALSE)
save.dir <-  "../../../occupancy_saved/saved/multiSpp-singleSea/saved"

survey.data <- read.csv("data/occupancy_data.csv")
species.groups <- read.csv("data/species_groups.csv")
survey.dates <- read.csv("data/survey_dates.csv")
habitat <- read.csv("data/habitat.csv")

## mcmc settings
scale <- 1e2
burnin <- 1e2*scale
niter <- (1e3)*scale


logit <- function(x) {
    log(x/(1 - x))
}

expit <- function(x) {
    exp(x)/(1 + exp(x))
}



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
    save(ms.ss.samples, file=file.path(save.dir,
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
    mclapply(MCMCs, runAllMCMC, input1, niter, burnin,  latent,
           hyper.param, MCMCdefs)
}
