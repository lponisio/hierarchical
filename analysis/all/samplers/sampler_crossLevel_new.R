binarySampler_baseClass <- nimbleFunctionVirtual(
  run = function() {},
  methods = list(
    getLogProbProposal = function() {returnType(double())},
    getCurrentLogProb = function() {returnType(double())},
    getLogProbLastProposal = function() {returnType(double())},
    resetValue = function(){},
    reset = function(){}
  ))


sampler_binary_new <- nimbleFunction(
  contains = binarySampler_baseClass,
  setup = function(model, mvSaved, target, control) {
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes  <- model$getDependencies(target)
    lastProb <- 0
    currentProb <- 0
    proposalProb <- 0
    lastValue <- 0
    ## checks
    if(length(targetAsScalar) > 1)  stop('cannot use binary sampler on more than one target node')
    if(!model$isBinary(target))     stop('can only use binary sampler on discrete 0/1 (binary) nodes')
  },
  run = function() {
    lastProb <<- exp(getLogProb(model, calcNodes))
    lastValue <<- model[[target]]
    model[[target]] <<- 1 - model[[target]]
    proposalProb <<- exp(calculate(model, calcNodes))
    if(!is.nan(proposalProb) & runif(1,0,1) < proposalProb/(lastProb+proposalProb)){
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
      currentProb <<- proposalProb
    }
    else{
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
      currentProb <<- lastProb
    }
  },
  methods = list(
    getLogProbProposal = function(){
      returnType(double(0))
      logProb <- log(proposalProb)
      return(logProb)
    },
    getCurrentLogProb = function(){
      returnType(double(0))
      logProb <- log(currentProb)
      return(logProb)
    },
    getLogProbLastProposal = function(){
      returnType(double(0))
      logProb <- log(lastProb)
      return(logProb)
    },
    resetValue = function(){
      model[[target]] <<- lastValue
      calculate(model, calcNodes)
    },
    reset = function() { }
  )
)

sampler_crossLevelBinary <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## control list extraction
    adaptive       <- control$adaptive
    adaptScaleOnly <- control$adaptScaleOnly
    adaptInterval  <- control$adaptInterval
    scale          <- control$scale
    propCov        <- control$propCov
    ## node list generation
    target       <- model$expandNodeNames(target)
    lowNodes     <- model$getDependencies(target, self = FALSE, stochOnly = TRUE, includeData = FALSE)
    lowCalcNodes <- model$getDependencies(lowNodes)
    calcNodes    <- model$getDependencies(c(target, lowNodes))
    ## nested function and function list definitions
    mvInternal <- modelValues(model)
    RWblockControl <- list(adaptive = adaptive, adaptScaleOnly = adaptScaleOnly, adaptInterval = adaptInterval, scale = scale, propCov = propCov)
    topRWblockSamplerFunction <- sampler_RW_block(model, mvInternal, target, RWblockControl)
    # lsf <- sampler_binary_new(model, mvSaved, lowNodes[1], control = list())
    lowSamplerFunctions <- nimbleFunctionList(binarySampler_baseClass)
    for(iLN in seq_along(lowNodes)) {
      lowNode <- lowNodes[iLN]
      lowSamplerFunctions[[iLN]] <- sampler_binary_new(model, mvSaved, lowNode, control = list())
    }
    my_setAndCalculateTop <- setAndCalculate(model, target)
    my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
  },
  run = function() {
    modelLP0 <- getLogProb(model, calcNodes) ## 
    for(iSF in seq_along(lowSamplerFunctions))
      { modelLP0 <- modelLP0 + lowSamplerFunctions[[iSF]]$getCurrentLogProb() } ## 1

    propValueVector <- topRWblockSamplerFunction$generateProposalVector() ## propose target'
    modelLP1 <- my_setAndCalculateTop$run(propValueVector)  ## 5
    propLP1 <- 0
    for(iSF in seq_along(lowSamplerFunctions))
      { lowSamplerFunctions[[iSF]]$run() ## run lower-level samplers
         propLP1 <- propLP1 + lowSamplerFunctions[[iSF]]$getLogProbProposal()## calculate 3
         modelLP1 <- modelLP1 +  lowSamplerFunctions[[iSF]]$getCurrentLogProb() } ## calculate 4

    nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    nimCopy(from = mvSaved, to = model, row = 1, nodes = lowNodes, logProb = TRUE)

    propLP0 <- 0
    for(iSF in seq_along(lowSamplerFunctions)){
      propLP0 <- propLP0 + lowSamplerFunctions[[iSF]]$getLogProbLastProposal()
    } ## calculate 6?!

    jump <- my_decideAndJump$run(modelLP1, modelLP0, propLP1, propLP0)
    if(adaptive)     topRWblockSamplerFunction$adaptiveProcedure(jump)
    if(!jump){
      for(iSF in seq_along(lowSamplerFunctions)){
        lowSamplerFunctions[[iSF]]$resetValue()
      }
    }
  },
  methods = list(
    reset = function() {
      topRWblockSamplerFunction$reset()
      for(iSF in seq_along(lowSamplerFunctions)) {
        lowSamplerFunctions[[iSF]]$reset()
      }
    }
  )
)


