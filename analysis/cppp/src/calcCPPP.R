library(parallel)
library(coda)

pppFuncVirtual <- nimbleFunctionVirtual(
    run = function(N = integer(0)) returnType(double(0))
)

virtualDiscFunction <- nimbleFunctionVirtual(
    run = function() returnType(double(0))
)


likeDiscFuncGenerator <- nimbleFunction(
    setup = function(model, ...){},
    run = function(){
        output <- -calculate(model)
        returnType(double(0))
        return(output)
    },
    contains = virtualDiscFunction
)

pppFunc <- nimbleFunction(
    setup = function(model,
                     dataNames,
                     paramNames,
                     MCMCIter,
                     thin,
                     burnIn,
                     discFunc,
                     nBootstrapSDReps,
                     discFunctionArgs){
        paramDependencies <- model$getDependencies(paramNames)
        dataDependencies <- model$getDependencies(dataNames)
        discFunction <- nimbleFunctionList(virtualDiscFunction)
        discFunction[[1]] <- discFunc(model, discFunctionArgs)
        l <- ceiling(min(1000, (MCMCIter/thin - burnIn)/20)) ##length of each
        ##block, ensures
        ##it's not too big
        ## total number of blocks available
        q <- (MCMCIter/thin - burnIn) - l + 1
        ##to sample from
        ## number of blocks to use for ppp
        h <- ceiling((MCMCIter/thin - burnIn)/l)
        ##function calculation
        blockMatrix <- matrix(0, nrow = h*l, ncol = length(paramNames))
        ## approximately of size m (number of mc
        ## samples)
    },
    run = function(MCMCOutput = double(2), nPPPCalcIters = integer(0), useBlockMatrix = logical(0)){
        output <- numeric(nPPPCalcIters)
        origDataValues <- values(model, dataNames)
        for(i in 1:nPPPCalcIters){
            randNum <- ceiling(runif(1, 0, (MCMCIter)/thin - burnIn - 1 ))
            values(model, dataNames) <<- origDataValues
            if(useBlockMatrix == TRUE){
                values(model, paramNames) <<- blockMatrix[randNum, ]
            }
            else{
                values(model, paramNames) <<- MCMCOutput[burnIn + randNum, ]
            }
            calculate(model, paramDependencies)
            calculate(model, dataDependencies)
            obsDeviance <- discFunction[[1]]$run()
            simulate(model, dataNames, includeData = TRUE)
            calculate(model, paramDependencies)
            calculate(model, dataDependencies)
            simDeviance <- discFunction[[1]]$run()
            if(simDeviance >= obsDeviance) output[i] <- 1
            else output[i] <- 0
        }
        out <- mean(output)
        returnType(double(0))
        return(out)
    },
    methods = list(
        getSD = function(MCMCOutput = double(2), nPPPCalcIters=double(0)){
            blockpppValues <- numeric(nBootstrapSDReps) ## block estimates of ppp
            for(r in 1:nBootstrapSDReps){
                for(i in 1:h){
                    randNum <- runif(1,0,1)
                    ##random starting index for blocks (post burn-in)
                    randIndex <- ceiling(randNum*q)
                    for(j in 1:l){
                        ## fill in blockMV with chosen blocks
                        blockMatrix[(i - 1)*l + j,] <<- MCMCOutput[burnIn + randIndex - 1 + j,]
                    }
                }
                ##as per Caffo, calculate both Q functions using the same
                ##samples from the latent variables
                blockpppValues[r]  <- run(MCMCOutput, nPPPCalcIters, TRUE)
            }
            blockpppValuesSD <- sd(blockpppValues)
            returnType(double())
            return(blockpppValuesSD)
        }
    ))




## calculate proportion of deviances that are greater or equal to the
## observed value
calcCPPP <- function(MCMCIter,
                     burnIn,
                     nPPPCalcIters,
                     C.pppFunc,
                     C.mcmc = NULL,
                     samples = NULL,
                     returnSamples){
    if(is.null(samples)){
        C.mcmc$run(MCMCIter, reset = TRUE, simulateAll = TRUE)
        samples <- as.matrix(C.mcmc$mvSamples)
    }
    pre.pp <- C.pppFunc$run(samples, nPPPCalcIters, FALSE)
    bootSD <- C.pppFunc$getSD(samples, nPPPCalcIters)
    if(!is.finite(pre.pp))    pre.pp <- NA
    return(list(pre.pp = pre.pp,
                bootSD = bootSD,
                samples = if(returnSamples == T) samples else NULL))
}


generateCPPP <-  function(R.model,
                          origMCMCOutput,
                          mcmcCreator = NULL,
                          cpppMCMCIter = NULL,
                          thin = 1,
                          dataNames, ## names of the data
                          nPPPCalcIters,## number of samples from posterior
                          nSimPPPVals, ## number of simulated PPP values
                          burnInProp, ## proportion of mcmc to drop
                          discFuncGenerator = NULL,
                          returnSamples = FALSE,
                          nBootstrapSDReps = 200, ## number of bootstrap samples
                          nCores = 1,
                          discFunctionArgs = list()){
    if(is.null(discFuncGenerator)){
        message("No discrepancy function specified, will use negative model likelihood.")
        discFuncGenerator <- likeDiscFuncGenerator
    }
    if(!inherits(R.model, "RmodelBaseClass")){
        stop("R.model is not an Rmodel")
    }

    if(is(R.model$CobjectInterface, "uninitializedField")){
        compileNimble(R.model)
    }
    if(burnInProp >= 1 | burnInProp < 0){
        stop("burnInProp needs to be between 0 and 1")
    }

    origMCMCIter <- nrow(origMCMCOutput)*thin

    if(origMCMCIter <= 1){
        stop("origMCMCOutput has no rows!")
    }
    if(is.null(cpppMCMCIter)){
        cpppMCMCIter <- origMCMCIter
    }
    else if(cpppMCMCIter < 2){
        stop("cpppMCMCIter must be an integer > 1")
    }

    burnIn <- ceiling(burnInProp*(cpppMCMCIter/thin))

    if(nPPPCalcIters > cpppMCMCIter){
        stop("number of samples from posterior must be < number of MCMC iterations")
    }

    testDataNames <- try(R.model[[dataNames]], silent=TRUE)
    if(inherits(testDataNames, "try-error")){
        stop(paste("dataNames", dataNames,
                   "is not the name of the data in model"))
    } else{
        test2DataNames <- all(R.model$expandNodeNames(dataNames) %in%
                              R.model$getNodeNames(dataOnly=TRUE))
        if(test2DataNames == FALSE){
            stop(paste("dataNames", dataNames,
                       "is not the name of the data in model"))
        }
    }

    paramNames <- colnames(origMCMCOutput)
    testParamNames <- lapply(paramNames, function(x){
        test.this.param <- try(R.model[[x]], silent=TRUE)
        if(inherits(test.this.param, "try-error")){
            stop(paste("paramNames", x,
                       "are not parameters in model"))
        }
    })
    test2ParamNames <- all(R.model$expandNodeNames(paramNames) %in%
                           R.model$getNodeNames(includeData=FALSE,
                                                stochOnly=TRUE))
    if(test2ParamNames == FALSE){
        stop(paste("paramNames", paramNames,
                   "are not parameters in model"))
    }


    ## keep track of the real data
    origData <- nimble:::values(R.model, dataNames)

    ## sample posterior, simulate data from sample
    paramDependencies <- R.model$getDependencies(paramNames)

    modelpppFunc <- pppFunc(R.model,
                            dataNames,
                            paramNames,
                            origMCMCIter,
                            thin,
                            burnIn,
                            discFunc = discFuncGenerator,
                            nBootstrapSDReps=nBootstrapSDReps,
                            discFunctionArgs)
    C.pppFunc <- compileNimble(modelpppFunc,
                               project = R.model)

    ## calculate deviances
    obs.cppp <- calcCPPP(origMCMCIter,
                         burnIn,
                         nPPPCalcIters=nPPPCalcIters,
                         C.pppFunc,
                         samples = origMCMCOutput,
                         returnSamples = returnSamples)


    cpppParallelFunction <- function(coreNum, nBootIter){
        R.newModel <- R.model$newModel()
        C.newModel <- compileNimble(R.newModel)
        if(is.null(mcmcCreator)){
            newMCMC.spec <-   configureMCMC(R.newModel,
                                            print=FALSE,
                                            nthin=thin)
            newMCMC <- buildMCMC(newMCMC.spec)
        }
        else{
            newMCMC <- mcmcCreator(R.newModel)
        }

        modelpppFunc <- pppFunc(R.newModel,
                                dataNames,
                                paramNames,
                                cpppMCMCIter,
                                thin,
                                burnIn,
                                discFunc = discFuncGenerator,
                                nBootstrapSDReps=nBootstrapSDReps,
                                discFunctionArgs)

        C.newMcmc <- compileNimble(newMCMC, project = R.newModel)
        C.newPppFunc <- compileNimble(modelpppFunc,
                                      project = R.newModel)
        out <- list()
        for(iteration in 1:nBootIter){
            message(paste("refitting data iteration", iteration, "for core number", coreNum))
            simulate(C.newModel,  includeData =  TRUE)
            out[[iteration]] <- calcCPPP(cpppMCMCIter,
                                         burnIn,
                                         nPPPCalcIters,
                                         C.pppFunc = C.newPppFunc,
                                         C.mcmc = C.newMcmc,
                                         returnSamples = returnSamples)
        }
        return(out)
    }


    ## simulate the ppp values
    sim.ppp.output <- mclapply(X = 1:nCores, FUN = cpppParallelFunction,
                               nBootIter =  ceiling(nSimPPPVals/nCores))

    ## extract simulated ppp and boot SDs
    sim.ppp <- sapply(sapply(sim.ppp.output,
                             function(x) x), function(x)  x$pre.pp)
    sim.vars <- (sapply(sapply(sim.ppp.output, function(x) x),
                        function(x) x$bootSD))^2

    ## approximate the distbution of observed ppp and simulated ppp
    diff.ppp <- obs.cppp$pre.pp - sim.ppp
    diff.vars.ppp <- obs.cppp$bootSD^2 + sim.vars


    ## calculate the average prob that simulated values are less than
    ## the observed
    cppp <- mean(pnorm(0, diff.ppp, diff.vars.ppp),
                 na.rm=TRUE)

    ## extract chains and diagnostics
    if(returnSamples){
        sim.samples <- lapply(sapply(sim.ppp.output, function(x) x),
                              function(x) x$samples)
    }
    else sim.samples <- NULL


    ## output real data and model
    nimble:::values(R.model, dataNames) <- origData

    out <- list(cppp=cppp,
                obs.ppp=c(estimate=obs.cppp$pre.pp,
                          bootVar=(obs.cppp$bootSD)^2),
                sim.cpp.dist=cbind(esimate=sim.ppp,
                                   bootVar=sim.vars),
                samples=sim.samples)
    return(out)
}
