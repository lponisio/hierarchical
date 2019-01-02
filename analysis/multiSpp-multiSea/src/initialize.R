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
source("../all/plotting.R")
source("src/plotting.R")
source("src/models.R")
source("src/customSamplerSpec.R")
source("src/dynamicOcc.R")

save.dir <-  "../../../hierarchical_saved/multiSpp-multiSea/saved"

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
        ms.ms.samples <- compareMCMCs(input1,
                                      MCMCs=i,
                                      niter=niter,
                                      burnin = burnin,
                                      summary=FALSE,
                                      check=FALSE)
    } else{
        ms.ms.samples <- compareMCMCs(input1,
                                      MCMCs=i,
                                      MCMCdefs = MCMCdefs[i],
                                      niter=niter,
                                      burnin = burnin,
                                      summary=FALSE,
                                      check=FALSE)

    }
    save(ms.ms.samples,
         file=file.path(save.dir,
                        sprintf("hyperparam%s_latent%s_sampler%s.Rdata",
                                hyper.param,
                                latent, i)))
}



runAllModels <- function(latent, hyper.param, niter, burnin,
                         MCMCs, MCMCdefs){
    load("data/all-5-0-2500-350.Rdata")
    if(!latent){
        model.input$data$Z <- NULL
        model.input$inits$Z <- NULL
        ## We do not want any X element equal to NA or they will not be
        ## considered data and will be sampled.
        model.input$data$X[ is.na(model.input$data$X) ] <- -1000
    }
    if(!hyper.param){
        model.input$inits <- lapply(model.input$inits, function(x){
            if(length(x) == 49){
                x <- x[1]
            } else{
                x <- x
            }
            return(x)
        })
        model.input$inits[grepl("sigma",
                                names(model.input$inits))] <- NULL
        model.input$inits[grepl("mu",
                                names(model.input$inits))] <- NULL
    }
    ms.ms.occ <- makeModel(latent, hyper.param)
    input1 <- c(code=ms.ms.occ,
                model.input)
    lapply(MCMCs, runAllMCMC, input1, niter, burnin,  latent,
           hyper.param, MCMCdefs)
}


