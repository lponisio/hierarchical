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
scale <- 1
burnin <- 1e1*scale
niter <- (1e3)*scale

logit <- function(x) {
    log(x/(1 - x))
}

expit <- function(x) {
    exp(x)/(1 + exp(x))
}



runAllMCMC <- function(i, input1, niter, burnin, ffilter,
                       hyper.param, MCMCdefs, mu.p){
    print(sprintf("hyperparam%s_filter%s_sampler%s_mup%s",
                                               hyper.param,
                                               ffilter, i, mu.p))
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
                                       sprintf("hyperparam%s_filter%s_sampler%s_mup%s.Rdata",
                                               hyper.param,
                                               ffilter, i, mu.p)))
}



runAllModels <- function(ffilter, hyper.param, niter, burnin, MCMCs,
    MCMCdefs, mu.p){
    data <- genDynamicOccData(mu.p=mu.p)
    model.input <- prepModDataOcc(data, include.zs=!ffilter)
    if(!hyper.param){
        to.drop <- c("sigma.phi", "sigma.gamma", "sigma.p", "p", "phi", "gamma")
        model.input$inits[to.drop] <- NULL
    }
    ss.ms.occ <- makeModel(ffilter, hyper.param)
    input1 <- c(code=ss.ms.occ,
                model.input)
    lapply(MCMCs, runAllMCMC, input1, niter, burnin,  ffilter,
           hyper.param, MCMCdefs, mu.p)
}
