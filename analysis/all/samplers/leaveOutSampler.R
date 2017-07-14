sampler_latentSub <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    leaveOutProportion <- control$leaveOutProportion
    allTargets <- model$expandNodeNames(target)
    latentNodeLength <- length(allTargets)
    latentSamplerList <- nimbleFunctionList(sampler_BASE)
    for(latentNode in 1:latentNodeLength){
      latentSamplerList[[latentNode]] <-
        sampler_binary(model,
                       mvSaved,
                       target =
                       allTargets[latentNode],
                       control = control$control)
    }
  },
  run = function() {
    for(latentNode in 1:latentNodeLength){
      sampleIndicator <- rbinom(1, 1, leaveOutProportion)
      if(sampleIndicator == 0)
      latentSamplerList[[latentNode]]$run()
    }
  },
  methods = list(
    reset = function() {}
    ))


