binarySampler_baseClass <- nimbleFunctionVirtual(
  run = function() {returnType(double())},
  methods = list(
    ## getLogProbProposal = function() {returnType(double())},
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
    currentProb <- model$calculate(target)
    proposalProb <- 0
    lastValue <- 0
    ## checks
    if(length(targetAsScalar) > 1)  stop('cannot use binary sampler on more than one target node')
    if(!model$isBinary(target))     stop('can only use binary sampler on discrete 0/1 (binary) nodes')
  },
  run = function() {
    lastProb <<- exp(getLogProb(model, calcNodes))
    lastValue <<- model[[target]]
    model[[target]] <<- 1 - lastValue
    proposalProb <<- exp(calculate(model, calcNodes))
    if(!is.nan(proposalProb) & runif(1,0,1) < proposalProb/(lastProb+proposalProb)){
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
      currentProb <<- proposalProb
    }
    else{
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
      currentProb <<- lastProb
    }
    return(log(proposalProb))
    returnType(double(0))
  },
  methods = list(
    ## getLogProbProposal = function(){
    ##   returnType(double(0))
    ##   logProb <- log(proposalProb)
    ##   return(logProb)
    ## },
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
    lowNodes     <- model$getDependencies(target, self = FALSE,
                                          stochOnly = TRUE,
                                          includeData = FALSE)
    lowCalcNodes <- model$getDependencies(lowNodes)
    targetCalcNodes    <- model$getDependencies(c(target))#, lowNodes))  This may be incorrect, let's try without lowNodes
    targetCalcNodes <- targetCalcNodes[!(targetCalcNodes %in% lowNodes)]
    ## nested function and function list definitions
    mvInternal <- modelValues(model)

    RWblockControl <- list(adaptive = adaptive,
                           adaptScaleOnly = adaptScaleOnly,
                           adaptInterval = adaptInterval,
                           scale = scale, propCov = propCov)
    topRWblockSamplerFunction <- sampler_RW_block(model, mvInternal,
                                                  target, RWblockControl)
    lowSamplerFunctions <- nimbleFunctionList(binarySampler_baseClass)
    for(iLN in seq_along(lowNodes)) {
      lowNode <- lowNodes[iLN]
      lowSamplerFunctions[[iLN]] <- sampler_binary_new(model, mvSaved,
                                                       lowNode, control = list())
    }
    my_setAndCalculateTop <- setAndCalculate(model, target)
    my_decideAndJump <- decideAndJump(model, mvSaved, targetCalcNodes)
  },
  run = function() {
    modelLP0 <- getLogProb(model, targetCalcNodes) ##
    ## latent node probabilities 
    for(iSF in seq_along(lowSamplerFunctions)){
        modelLP0 <- modelLP0 +
                           lowSamplerFunctions[[iSF]]$getCurrentLogProb()
      } ## 1
    propValueVector <- topRWblockSamplerFunction$generateProposalVector() ## propose target'
    modelLP1 <- my_setAndCalculateTop$run(propValueVector)  ## 5
    propLP1 <- 0
    for(iSF in seq_along(lowSamplerFunctions)){
        propLP1 <- propLP1 +
          lowSamplerFunctions[[iSF]]$run()## calculate
      ## 3
        modelLP1 <- modelLP1 +
                           lowSamplerFunctions[[iSF]]$getCurrentLogProb()
    } ## calculate 4
    nimCopy(from = mvSaved, to = model, row = 1,
            nodes = targetCalcNodes, logProb = TRUE)
    nimCopy(from = mvSaved, to = model, row = 1,
            nodes = lowNodes, logProb = TRUE)

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





binary2_virtual <- nimbleFunctionVirtual(
    run = function() { },
    methods = list(
        calcBinaryPosteriorLogDensity = function() { returnType(double()) },
        reset = function() { }
    )
)

sampler_binary2 <- nimbleFunction(
    name = 'sampler_binary2',
    contains = binary2_virtual,
    setup = function(model, mvSaved, target, control) {
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
        ## checks
        if(length(targetAsScalar) > 1)  stop('cannot use binary sampler on more than one target node')
        if(!model$isBinary(target))     stop('can only use binary sampler on discrete 0/1 (binary) nodes')
    },
    run = function() {
        currentLogProb <- getLogProb(model, calcNodes)
        model[[target]] <<- 1 - model[[target]]
        otherLogProb <- calculate(model, calcNodes)
        acceptanceProb <- 1/(exp(currentLogProb - otherLogProb) + 1)
        if(!is.nan(acceptanceProb) & runif(1,0,1) < acceptanceProb)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        else
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    },
    methods = list(
        calcBinaryPosteriorLogDensity = function() {
            currentValue <- model[[target]]
            currentLogProb <- getLogProb(model, calcNodes)
            model[[target]] <<- 1 - currentValue
            otherLogProb <- calculate(model, calcNodes)
            model[[target]] <<- currentValue
            calculate(model, calcNodes)
            lp <- -log(exp(otherLogProb - currentLogProb) + 1)
            returnType(double(0))
            return(lp)
        },
        reset = function() { }
    )
)

sampler_crossLevel_binary_DT <- nimbleFunction(
    name = 'sampler_crossLevel_binary_DT',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive <- if(!is.null(control$adaptive)) control$adaptive else TRUE
        ## node list generation
        target       <- model$expandNodeNames(target)
        lowNodes     <- model$getDependencies(target, self = FALSE, stochOnly = TRUE, includeData = FALSE)
        lowCalcNodes <- model$getDependencies(lowNodes)
        calcNodes    <- model$getDependencies(c(target, lowNodes))
        ## nested function and function list definitions
        mvInternal <- modelValues(model)
        topRWblockSamplerFunction <- sampler_RW_block(model, mvInternal, target, control)
        lowBinarySamplerFunctions <- nimbleFunctionList(binary2_virtual)
        for(iLN in seq_along(lowNodes)) {
            lowNode <- lowNodes[iLN]
            if(!model$isBinary(lowNode))     stop('non-binary lowNode \'', lowNode, '\' in crossLevel_binary_DT sampler')
            lowBinarySamplerFunctions[[iLN]] <- sampler_binary2(model, mvSaved, lowNode, control)
        }
        my_setAndCalculateTop <- setAndCalculate(model, target)
        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
    },
    run = function() {
        modelLP0 <- getLogProb(model, calcNodes)
        propLP0 <- 0
        for(iSF in seq_along(lowBinarySamplerFunctions))  { propLP0 <- propLP0 + lowBinarySamplerFunctions[[iSF]]$calcBinaryPosteriorLogDensity() }
        propValueVector <- topRWblockSamplerFunction$generateProposalVector()
        topLP <- my_setAndCalculateTop$run(propValueVector)
        if(is.na(topLP))
            jump <- my_decideAndJump$run(-Inf, 0, 0, 0)
        else {
            for(iSF in seq_along(lowBinarySamplerFunctions))
                lowBinarySamplerFunctions[[iSF]]$run()
            modelLP1 <- calculate(model, calcNodes)
            propLP1 <- 0
            for(iSF in seq_along(lowBinarySamplerFunctions))
                propLP1 <- propLP1 + lowBinarySamplerFunctions[[iSF]]$calcBinaryPosteriorLogDensity()
            jump <- my_decideAndJump$run(modelLP1, modelLP0, propLP1, propLP0)
    	}
        if(adaptive)     topRWblockSamplerFunction$adaptiveProcedure(jump)
    },
    methods = list(
        reset = function() {
            topRWblockSamplerFunction$reset()
            for(iSF in seq_along(lowBinarySamplerFunctions)) {
                lowBinarySamplerFunctions[[iSF]]$reset()
            }
        }
    )
)
