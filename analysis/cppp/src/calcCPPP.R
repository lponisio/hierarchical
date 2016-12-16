library(parallel)
library(coda)

## asymptotic variance is calculated using the moving-block
## bootstrap method of "Markov Chain Monte Carlo in Statistical Mechanics"
## by Mignani & Rosa, 2001 (p. 350)
## model: the C model
## C.pppFunc: the function that calculates the ppp
## fixedNodes: (may not need)
## sampledNodes: parameters that we want to sample
## mvBlock: (may not need)
## mvSqmplew: output from MCMC previously run
## burn in is the burn in
## numReps: the number of iternations of the moving block bootstrap

## nsamps : number of mcmc iterations
## theta : 
## calc_asympSD = nimbleFunction(
##   setup = function(model, sampledNodes, mvSample,
##     mvBlock, burnIn = 0, numReps=200,
##     dataNames,
##     MCMCIter,
##     thin,
##     averageParams,
##     discFunc,
##     ...){

##    PPPFunction <- pppFunc(R.model, dataNames, sampledNodes,
##                               mvBlock,
##                               MCMCIter,
##                               thin,
##                               burnIn,
##                               averageParams,
##                               discFunc=discFunc,
##                               ...)

##   ## pppFunction <- nimbleFunctionList(pppFuncVirtual)
##   ## pppFunction[[1]] <- bootstrapPPPFunc 

##   },

##   run = function(nsamps = double(0)){
##     blockpppValues <- numeric(numReps) ## block estimates of ppp
##     l <- ceiling(min(1000, (nsamps - burnIn)/20)) ##length of each
##     ##block, ensures
##     ##it's not too big
##     q <- (nsamps - burnIn) - l + 1 ##total number of blocks available
##     ##to sample from
##     h <- ceiling((nsamps - burnIn)/l) ##number of blocks to use for ppp
##     ##function calculation
##     resize(mvBlock, h*l) ##size our model value object to be
##     ##approximately of size m (number of mc
##     ##samples)

##     for(r in 1:numReps){
##       for(i in 1:h){
##         randNum <- runif(1,0,1)
##         randIndex <- ceiling(randNum*q) ##random starting index for
##         ##blocks (post burn-in)
##         for(j in 1:l){
##           nimCopy(mvSample, mvBlock, sampledNodes, sampledNodes,
##                   burnIn + randIndex-1+j,  (i-1)*l+j) ##fill in mvBlock with chosen blocks
##         }
##       }
##       ##as per Caffo, calculate both Q functions using the same
##       ##samples from the latent variables
##       blockpppValues[r]  <- PPPFunction$run(nsamps) 
##     }

##     blockpppValuesSD <- sd(blockpppValues)
##     returnType(double())
##     return(blockpppValuesSD)
##   },
##   where = getLoadingNamespace()
##   )

pppFuncVirtual <- nimbleFunctionVirtual(
  run = function(N = integer(0)) returnType(double(0))
  )

virtualDiscFunction <- nimbleFunctionVirtual(
  run = function() returnType(double(0))
  )

likeDiscFuncGenerator <- nimbleFunction(
  setup = function(model, ...){},
  run = function(){
    output <- calculate(model)
    returnType(double(0))
    return(output)
  },
  contains = virtualDiscFunction
  )


maxDiscFuncGenerator <- nimbleFunction(
  setup = function(model, ...){
    params <- list(...)
    dataNames <- params[[1]]
  },
  run = function(){
    output <- max(values(model, dataNames))
    returnType(double(0))
    return(output)
  },
  contains = virtualDiscFunction
  )


pppFunc <- nimbleFunction(
  setup = function(
    model,
    dataNames,
    paramNames,
    mcmcMV,
    MCMCIter,
    thin,
    burnIn,
    averageParams,
    discFunc,
    numReps,
    mvBlock,
    ...){
    
    paramDependencies <- model$getDependencies(paramNames)
    discFunction <- nimbleFunctionList(virtualDiscFunction)
    discFunction[[1]] <- discFunc(model,...) 
  },
  
  run = function(N = integer(0)){
    output <- numeric(N)
    if(averageParams == 0){
      discMean <- discFunction[[1]]$run()
    }
    else{
      mcmcSamps <-floor(MCMCIter/thin - burnIn)
      discMean <- 0
      for(i in 1:mcmcSamps){
        copy(mcmcMV, model, paramNames, paramNames, row = i + burnIn)
        discMean <- discMean + discFunction[[1]]$run()

      }
      discMean <- discMean / mcmcSamps
      
    }
    if(is.nan(discMean)) return(NA)
    for(i in 1:N){
      randNum <- ceiling(runif(1, 0, (MCMCIter)/thin - burnIn - 1 ))
      nimCopy(mcmcMV, model, paramNames, row = burnIn + randNum)
      calculate(model, paramDependencies)
      simulate(model, dataNames, includeData = TRUE)
      deviance <- discFunction[[1]]$run()
      if(deviance >= discMean) output[i] <- 1
      else output[i] <- 0
    }
    out <- mean(output)
    returnType(double(0))
    return(out)
  },
  methods = list(
    getSD = function(nsamps=double(0)){
      blockpppValues <- numeric(numReps) ## block estimates of ppp
      l <- ceiling(min(1000, (nsamps - burnIn)/20)) ##length of each
      ##block, ensures
      ##it's not too big
      q <- (nsamps - burnIn) - l + 1 ##total number of blocks available
      ##to sample from
      h <- ceiling((nsamps - burnIn)/l) ##number of blocks to use for ppp
      ##function calculation
      resize(mvBlock, h*l) ##size our model value object to be
      ##approximately of size m (number of mc
      ##samples)

      for(r in 1:numReps){
        for(i in 1:h){
          randNum <- runif(1,0,1)
          randIndex <- ceiling(randNum*q) ##random starting index for
          ##blocks (post burn-in)
          for(j in 1:l){
            nimCopy(mcmcMV, mvBlock, paramNames, paramNames,
                    burnIn + randIndex-1+j,  (i-1)*l+j) ##fill in mvBlock with chosen blocks
          }
        }
        ##as per Caffo, calculate both Q functions using the same
        ##samples from the latent variables
        blockpppValues[r]  <- run(nsamps) 
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
                     NSamp,
                     C.pppFunc,
                     cppp.C.mcmc,
                     firstRun,
                     runUntilConverged = NULL,
                     maxIter = NULL,
                     convStep= NULL){
  convergeTest <- NA
  samples <- NA
  if(firstRun  == 0){
    cppp.C.mcmc$run(MCMCIter)
    samples <- mcmc(as.matrix(cppp.C.mcmc$mvSamples)[-c(1:burnIn),])

    if(runUntilConverged == 1){
      ## use Geweke diagnostic to see if MCMC has converged
      convergeTest <- geweke.diag(samples)$z
      ## any na values get set to a very large z value so the MCMC
      ## will continue to be run. note that nan values correspond to
      ## posterior samples that are constant, indicating a lack of
      ## convergence
      convergeTest[!is.finite(convergeTest)] <- 10 
      zVal <- 1.96
      while(any(abs(convergeTest) > zVal)){
        cppp.C.mcmc$run(MCMCIter*convStep, reset = FALSE)
        this.niter <- MCMCIter + convStep*MCMCIter
        samples <- mcmc(as.matrix(cppp.C.mcmc$mvSamples)[-c(1:burnIn),])
        convergeTest <- geweke.diag(samples)$z
        convergeTest[!is.finite(convergeTest)] <- 10 
        if(this.niter > maxIter)
          convergeTest <- 0
      }
    }
  }
  pre.pp <- C.pppFunc$run(NSamp)
  pppSD <- C.pppFunc$getSD(NSamp)
  if(!is.finite(pre.pp))    pre.pp <- NA
  return(list(pre.pp = pre.pp,
              pppSD = pppSD,
              samples = samples,
              converge.stat = convergeTest))    
}


generateCPPP <-  function(R.model,
                          orig.C.model,
                          orig.C.mcmc,
                          orig.mcmc,
                          dataNames, ## names of the data column
                          paramNames, ## vector of parameters to monitor
                          NSamp,## number of samples from posterior
                          NPDist, ## number of simulated PPP values
                          burnInProp, ## proportion of mcmc to drop
                          discFuncGenerator,
                          averageParams,
                          returnChains = TRUE,
                          runUntilConverged = FALSE,
                          maxIter = 1*10^4,
                          convStep = 0.5,
                          nRepBoot,
                          ...){
  if(!inherits(R.model, "RmodelBaseClass")){
    stop("R.model is not an Rmodel")
  }
  if(!inherits(orig.C.model, "CmodelBaseClass")){
    stop("orig.C.model is not an Cmodel")
  }
  if(burnInProp >= 1 | burnInProp < 0){
    stop("burnInProp needs to be between 0 and 1")
  }

  thin <- orig.C.mcmc$thin
  MCMCIter <- nrow(as.matrix(orig.C.mcmc$mvSamples))*thin

  if(MCMCIter <= 1){
    stop("MCMC must be run on C model before assessment")
  }

  burnIn <- ceiling(burnInProp*(MCMCIter/thin))

  if(NSamp > MCMCIter){
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
  origData <- nimble:::values(orig.C.model, dataNames)

  ## sample posterior, simulate data from sample 
  paramDependencies <- orig.C.model$getDependencies(paramNames)
  mcmcMV <- orig.mcmc$mvSamples

  mvBlock <- modelValues(R.model)

  modelpppFunc <- pppFunc(R.model,
                          dataNames,
                          paramNames,
                          mcmcMV,
                          MCMCIter,
                          thin,
                          burnIn,
                          averageParams,
                          discFunc = discFuncGenerator,
                          numReps=nRepBoot,
                          mvBlock=mvBlock)
  C.pppFunc <- compileNimble(modelpppFunc,
                             project = R.model)

  ## calculate deviances
  obs.cppp <- calcCPPP(MCMCIter,
                       burnIn,
                       NSamp,
                       C.pppFunc,
                       orig.C.mcmc,
                       firstRun = 1)

  ## refits model with sampled data, reruns, enter inner loop,
  ## calculates distbution of PPPs
  simPppDist <- function(iteration){
    message(paste("refitting data iteration", iteration))
    simulate(orig.C.model,  includeData =  TRUE)
    out <- calcCPPP(MCMCIter,
                    burnIn,
                    NSamp,
                    C.pppFunc,
                    orig.C.mcmc,
                    firstRun = 0,
                    runUntilConverged=runUntilConverged,
                    maxIter = maxIter,
                    convStep=convStep)
    return(out)
  }

  sim.cppp <- mclapply(1:NPDist, simPppDist)
  sim.ppp <- sapply(sim.cppp, function(x) x$pre.pp)
  sim.SDs <- sapply(sim.cppp, function(x) x$pppSD)
  sim.CI <- cbind(sim.ppp + sim.SDs*1.96, sim.ppp - sim.SDs*1.96)

  obs.CI <- cbind(obs.cppp$pre.pp + obs.cppp$pppSD*1.96,
              obs.cppp$pre.pp - obs.cppp$pppSD*1.96)

  lims <- function(f1, f2){
    mapply(f1, apply(obs.CI, 1, f2), apply(sim.CI, 1, f2))
  }
  
  overlap <- lims(min, max) - lims(max, min)
  ## CI do not overlap if the above is negative
  overlap[overlap < 0] <- 0

  ## for the simulated ppp that do not overlap but are greater than
  ## the oberved value
  overlap[max(obs.CI) < apply(sim.CI, 1, min)] <- 1

  ## simulate cppp
  sim.sim.cppp <- sapply(1:NSamp, function(x){
    mean(rbinom(length(overlap), size=1, overlap))
  })
                         
  sim.samples <- lapply(sim.cppp, function(x) x$samples)
  chain.diag <- lapply(sim.cppp, function(x) x$converge.stat)

  nimble:::values(orig.C.model, dataNames) <- origData 

  out <- list(cppp=quantile(sim.sim.cppp, c(0.025, 0.5, 0.975)),
              obs.ppp=c(estimate=obs.cppp$pre.pp,
                bootSD=obs.cppp$pppSD),
              sim.cpp.dist=cbind(esimate=sim.ppp,
                bootSD=sim.SDs),
              chain.diagnostics= chain.diag,
              samples=sim.samples)
  if(!returnChains){
    out$samples <- NULL
  }
  return(out)
}

