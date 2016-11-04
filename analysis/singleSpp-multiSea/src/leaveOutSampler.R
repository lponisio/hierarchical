RW_sampler_latentSubsamp <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    leaveOutProportion <- control$leaveOutProportion
    latentNodeLength <- length(model[[target]])

    latentSamplerList <- nimbleFunctionList(sampler_BASE)
    allTargets <- model$expandNodeNames(target)
    for(latentNode in 1:latentNodeLength){
      latentSamplerList[[latentNode]] <- sampler_binary(model, mvSaved, target = allTargets[latentNode], control = control$control) 
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


