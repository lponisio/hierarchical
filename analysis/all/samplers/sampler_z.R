
## *********************************************************************
## custom sampler for zs
## *********************************************************************
custom_z_sampler <- nimbleFunction(
  contains = sampler_BASE,      ## required (class inheritance)
  setup = function(model, mvSaved, target, control) { ## required arguments
    calcNodes <- model$getDependencies(target)
    ## to be used for logProbs for values of 0 and 1 respectively
    logProbs <- c(0.0, 0.0)
  },
  run = function() {
    ## get the current value
    currentValue <- model[[target]]
    ## and the related sum of logProb's
    currentLogProb <- getLogProb(model, calcNodes)
    for(i in 1:2) {   ## check both possible values
      ## for the current value, record the current sum of logProb
      if(i - 1 == currentValue)
        logProbs[i] <<- currentLogProb
      else {
        ## for the other value, put it in the model
        model[[target]] <<- i - 1
        ## and calculate the sum of logProb's
        logProbs[i] <<- calculate(model, calcNodes)
      }
    }
    ## pick a random number between 0 and 1
    u <- runif(1, 0, 1)
    ## choose 0 with this probability, 1 otherwise
    if(u < exp(logProbs[1])/sum(exp(logProbs[1:2])))
      newValue <- 0
    else
      newValue <- 1

    if(newValue != currentValue) ## update the saved model states
      copy(from = model,
           to = mvSaved, row = 1,
           nodes = calcNodes, logProb = TRUE)
    else
      copy(from = mvSaved,
           to = model, row = 1,
           nodes = calcNodes, logProb = TRUE)
  },
  ## required for the MCMC system, but it doesn't have to do anything
  methods = list(
    reset = function () {}
    )
  )
