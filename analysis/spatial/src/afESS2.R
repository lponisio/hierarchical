getParVirtual <- nimbleFunctionVirtual(
  methods = list(
    scalarGet = function(){returnType(double(0))}
  )
)

doPars2 <- nimbleFunction(
  contains = getParVirtual,
  setup = function(parName, mvSaved) {
  },
  methods = list(
    scalarGet = function(){
      paramVal <- mvSaved[parName, 1][1]
      returnType(double(0))
      return(paramVal)
    }
  ), where = getLoadingNamespace())


sampler_AFSS_to_RW_block <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    AFSS_sampler <- sampler_AF_slice(model, mvSaved, target = target, control = control$AF_sliceControl)
    RW_block_sampler <- sampler_RW_block(model, mvSaved, target = target, control = control$RWcontrol)
    nESSSamps <- control$AF_sliceControl$factorAdaptInterval
    RWAdaptInterval <- control$RWcontrol$adaptInterval
    nAFSSIters <- max(control$nAFSSIters, control$factorAdaptInterval)
    essThreshold <- control$essThreshold
    numParamVars <- control$numParamVars
    numESSAdaptations <- control$numESSAdaptations
    timeSwitch        <- control$timeSwitch
    timesRan      <- 0
    nESSReps <- 200
    ## numParamVars <- length(target)
    l <- ceiling(min(1000, nESSSamps/20)) #length of each block, ensures it's not too big
    q <- nESSSamps - l + 1 #total number of blocks available to sample from
    h <- ceiling(nESSSamps/l) #number of blocks to use for q function calculation
    doVarList <- nimbleFunctionList(getParVirtual)
    for(i in 1:numParamVars){
      doVarList[[i]] <- doPars2(target[i], mvSaved)
    }
    storeSamples <- matrix(0, nrow = numParamVars, ncol = nESSSamps)
    essMeanSD   <- rep(0, numParamVars)
    essSD <- rep(0, numParamVars)
    ESS <-  rep(0, numParamVars)
    useRWSampler <- 0
    timesAdapted <- 0
    essThreshold <- essThreshold*nESSSamps
    meanMatrix <-  matrix(0, nrow = numParamVars, ncol = q)
  },
  run = function() { ## Now there is no indicator variable, so check if the target node is exactly
    timesRan <<- timesRan + 1
    if(timeSwitch == 1){
      if(timesRan > nAFSSIters){
        RW_block_sampler$run()
      }
      else{
        AFSS_sampler$run()
        if(timesRan == nAFSSIters){
          copyAdaptiveParams()
        }
      }
    }
    else{
      if(useRWSampler == 1){
        RW_block_sampler$run()
      }
      else{
        AFSS_sampler$run()
        calcESSandAdapt()
      }
    }
  },
  methods = list(
    copyAdaptiveParams = function(){
      RW_block_sampler$propCov <<- AFSS_sampler$empirCov
    },
    calcESSandAdapt = function(){
      for(par in 1:numParamVars){
        storeSamples[par, timesRan] <<- doVarList[[par]]$scalarGet()
      }
      if(timesRan == nESSSamps){
        essCalcs <- matrix(0, nrow = numParamVars, ncol = nESSReps)
        timesAdapted <<- timesAdapted + 1
        for(par in 1:numParamVars){
          for(r in 1:nESSReps){
            for(i in 1:h){
              randNum <- rbeta(1,1,1)
              randIndex <- ceiling(randNum*q) #random starting index for blocks (post burn-in)
              for(j in 1:l){
                essCalcs[par, r] <- essCalcs[par, r] + storeSamples[par, randIndex + j - 1]
              }
            }
            essCalcs[par, r]  <-  essCalcs[par, r]/(h*l)
          }
          essMeanSD[par] <<- sd(essCalcs[par, ])
          essSD[par] <<- sd(storeSamples[par,])
          ESS[par] <<- (essSD[par]/essMeanSD[par])^2
        }
        minESS <- min(ESS)
        if(minESS > essThreshold | timesAdapted > numESSAdaptations){
          copyAdaptiveParams()
          useRWSampler <<- 1
        }
        timesRan <<- 0
      }
    },
    reset = function() {
      RW_block_sampler$reset()
      AFSS_sampler$reset()
      timesRan <<- 0
      useRWSampler <<- 0
      timesAdapted <<- 0
    }
  ))
