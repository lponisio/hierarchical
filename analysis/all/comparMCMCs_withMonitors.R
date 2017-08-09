compareMCMCs_withMonitors <- function(modelInfo, MCMCs = c('nimble'), MCMCdefs, monitors, BUGSdir, stanDir, stanInfo, doSamplePlots = FALSE, verbose = TRUE, summary = TRUE, ...) {
  
  ## stanInfo = list(codeFile = required, dir = optional,  stanParameterRules = optional, modelName = optional, data = optional, inits = optional)
  ## or simply a character codeFile
  
  ## set up list of models from 3 possible input formats
  if(is.character(modelInfo)) {
    modelContents <- list()
    models <- modelInfo
    for (i in seq_along(modelInfo)) {
      thisBUGSdir <- if(missing(BUGSdir)) getBUGSexampleDir(modelInfo[i]) else BUGSdir 
      modelContents[[i]] <- readBUGSmodel(model = modelInfo[i],
                                          dir = thisBUGSdir,
                                          returnComponents = TRUE)
      names(modelContents[[i]])[names(modelContents[[i]])=="model"] <- "code"
    }
  } else {
    if(!is.list(modelInfo)) stop('modelInfo must be a list if it is not a vector of BUGS example names')
    if(!is.list(modelInfo[[1]])) {
      modelContents <- list(modelInfo)
      if('name' %in% names(modelInfo)) names(modelContents) <- modelInfo$name
    }
    else modelContents <- modelInfo
    
    inputNames <- names(modelContents)
    if(is.null(inputNames)) models <- paste0('model', seq_along(modelContents))
    else {
      iNamesNeeded <- which(inputNames == "")
      models <- inputNames
      models[iNamesNeeded] <- paste0('model', iNamesNeeded)
    }
  }
  
  ## At this point models is a character vector of model names and modelContents is a list with information about each model
  
  results <- list()
  
  useStan <- 'stan' %in% MCMCs
  if(useStan) {
    useStanInfo <- !is.null(stanInfo)
    if(useStanInfo) {
      if(is.list(stanInfo)) if(!is.list(stanInfo[[1]])) stanInfo <- list(stanInfo)
      if(is.character(stanInfo)) stanInfo <- as.list(stanInfo)
      if(is.null(names(stanInfo))) names(stanInfo) <- models
    }
  } else {
    stanNameMaps <- stanDataFile <- stanInitFile <- NULL
  }
  
  noConjDef <- list(noConj = quote({ configureMCMC(Rmodel, useConjugacy=FALSE) }))
  if(missing(MCMCdefs)) MCMCdefs <- noConjDef
  else MCMCdefs <- c(MCMCdefs, noConjDef)
  
  for (i in 1:length(models)){
    if(verbose) cat(paste('Working on', models[i],'\n'))
    
    if(useStan) {
      if(!useStanInfo) {
        ## no stanInfo so use the model name for the stanModelName and stanCodefile
        stanModelName <- stanCodeFile <- models[i]
      } else {
        ## yes stanInfo so look at it
        thisStanInfo <- stanInfo[[models[i]]]
        if(is.null(thisStanInfo))
          ## no stanInfo for this 
          stanModelName <- stanCodeFile <- models[i]
        else {
          if(is.character(thisStanInfo)) {
            stanModelName <- stanCodeFile <- thisStanInfo
          } else {
            stanModelName <- thisStanInfo[[ 'modelName' ]]
            if(is.null(stanModelName)) stanModelName <- models[i]
            
            stanNameMaps <- thisStanInfo[[ 'stanParameterRules' ]]
            if(is.null(stanNameMaps)) stanNameMaps <- list()
            
            usingStanDir <- thisStanInfo[[ 'dir']]
            if(is.null(usingStanDir)) usingStanDir <- stanDir
            
            stanCodeFile <- if(is.null(names(thisStanInfo))) thisStanInfo[[ 1 ]] else thisStanInfo[['codeFile']]
            if(is.null(stanCodeFile)) stanCodeFile <- stanModelName
            stanCodeFile <- file.path(usingStanDir, stanCodeFile)
            
            if(!grepl('\\.stan', stanCodeFile)) stanCodeFile <- paste0(stanCodeFile, '.stan')
            
            
            if(is.list(thisStanInfo[[ 'data' ]])) stanDataFile <- thisStanInfo[[ 'data' ]]
            else if(is.null( thisStanInfo[[ 'data' ]])) stanDataFile <- NULL
            else stanDataFile <- file.path(usingStanDir, thisStanInfo[[ 'data' ]]) ## it's ok if this is NULL
            
            if(!is.null(stanDataFile) & !is.list(stanDataFile))
              if(!grepl('\\.data\\.R', stanDataFile)) stanDataFile <- paste0(stanCodeFile, '.data.R')
            
            if(is.list(thisStanInfo[[ 'inits' ]])) stanInitFile <- thisStanInfo[[ 'inits' ]]
            else if(is.null( thisStanInfo[[ 'inits' ]])) stanInitFile <- NULL
            else stanInitFile <- file.path(usingStanDir, thisStanInfo[[ 'inits' ]]) ## it's ok if this is NULL
            
            if(!is.null(stanInitFile) & !is.list(stanInitFile))
              if(!grepl('\\.init\\.R', stanInitFile)) stanInitFile <- paste0(stanInitFile, '.init.R')
          }
        }
      }
    }
    if(is.null(modelContents[[i]][['constants']])) {
      constants <- modelContents[[i]]$data
      data <- list()
    } else {
      constants <- modelContents[[i]]$constants
      data <- modelContents[[i]]$data
    }
    suite_output <- MCMCsuite(modelContents[[i]]$code, constants = constants, data = data, inits = modelContents[[i]]$inits,
                              monitors = monitors
                              ,setSeed = FALSE
                              ,MCMCs = MCMCs, makePlot = doSamplePlots, savePlot = doSamplePlots
                              ,summaryStats=c('mean','median','sd','CI95_low','CI95_upp')
                              ,calculateEfficiency=TRUE
                              #change
                              ,MCMCdefs = MCMCdefs 
                              ,stan_model=if(useStan) stanCodeFile else ""
                              ,stanNameMaps = stanNameMaps
                              ,stan_data = stanDataFile
                              ,stan_inits = stanInitFile
                              ,  ...
    )
    if(summary) 
      results[models[i]][[1]] <- list(summary = suite_output$summary,
                                      timing = suite_output$timing, ##create_time_df(suite_output$timing, length(MCMCs)),
                                      efficiency = suite_output$efficiency)
    else
      results[models[i]][[1]] <- suite_output
  }
  return(results)
}
