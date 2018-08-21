library(devtools)
## install_github("nimble-dev/nimble",
##                ref = "devel",
##                subdir = "packages/nimble")


## install_github("nimble-dev/nimble",
##                ref = "calcPriorFirst",
##                subdir = "packages/nimble")



library(nimble)
library(igraph)
library(parallel)
source("src/dynamicOcc.R")
source("src/setup.R")
source("../all/plotting.R")

source('src/models.R')
source('src/customSamplerSpec.R')

dir.create(file.path("../../../occupancy_saved/saved/singleSpp-multiSea/saved"),
           showWarnings = FALSE)
save.dir <-  "../../../occupancy_saved/saved/singleSpp-multiSea/saved"

## MCMC settings
scale <- 1e1
burnin <- 1e1*scale
niter <- (1e3)*scale

logit <- function(x) {
    log(x/(1 - x))
}

expit <- function(x) {
    exp(x)/(1 + exp(x))
}



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
    save(ss.ms.samples, file=file.path(save.dir,
                                       sprintf("hyperparam%s_latent%s_sampler%s_mup%s.Rdata",
                                               hyper.param,
                                               latent, i, mu.p)))
}



runAllModels <- function(latent, hyper.param, niter, burnin, MCMCs,
                         MCMCdefs, mu.p){
        data <- genDynamicOccData(mu.p=mu.p)
        model.input <- prepModDataOcc(data, include.zs=latent)
        if(!hyper.param){
            to.drop <- c("sigma.phi", "sigma.gamma", "sigma.p", "p", "phi", "gamma")
            model.input$inits[to.drop] <- NULL
        }
        ss.ms.occ <- makeModel(latent, hyper.param)
        input1 <- c(code=ss.ms.occ,
                    model.input)
        lapply(MCMCs, runAllMCMC, input1, niter, burnin,  latent,
               hyper.param, MCMCdefs, mu.p)
    }
