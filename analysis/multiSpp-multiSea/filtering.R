## ************************************************************
## setwd('~/Dropbox/occupancy')
rm(list=ls())
setwd('analysis/multiSpp-multiSea')

source("src/initialize.R")
load("data/5-0-350.Rdata")
source("src/complete_allInt_filtering.R")

input1 <- c(code=ms.ms.occ,
            model.input)

## *****************************************************************
## vanilla nimble/jags
## *****************************************************************

ms.ms.nimble <- compareMCMCs_withMonitors(input1,
                                          MCMCs=c('nimble'),
                                          niter=niter,
                                          burnin = burnin,
                                          summary=FALSE,
                                          check=FALSE,
                                          monitors=model.input$monitors)

save(ms.ms.orig, file=file.path(save.dir, 'filtering.Rdata'))

