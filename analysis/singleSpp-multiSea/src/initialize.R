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


runAllMCMC <- function(i, input1, niter, burnin, ffilter,  hyper.param){
    print(i)
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
                                       sprintf("hyperparam%s_filter%s_sampler%s.Rdata",
                                               hyper.param,
                                               ffilter, i)))
}



runAllModels <- function(ffilter, hyper.param, niter, burnin){
    data <- genDynamicOccData()
    model.input <- prepModDataOcc(data, include.zs=!ffilter)
    if(!hyper.param){
        to.drop <- c("sigma.phi", "sigma.gamma", "sigma.p", "p", "phi", "gamma")
        model.input$inits[to.drop] <- NULL
    }
    ss.ms.occ <- makeModel(ffilter, hyper.param)
    input1 <- c(code=ss.ms.occ,
                model.input)
    lapply(MCMCs, runAllMCMC, input1, niter, burnin,  ffilter,  hyper.param)
}
