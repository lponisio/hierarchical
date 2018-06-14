## ************************************************************
## setwd('~/Dropbox/occupancy')
rm(list=ls())
setwd('analysis/multiSpp-multiSea')

source("src/initialize.R")
load("data/all-5-0-2500-350.Rdata")
source("src/complete_allInt_filter.R")

model.input$data$Z <- NULL
model.input$inits$Z <- NULL
## We do not want any X element equal to NA or they will not be
## considered data and will be sampled.
model.input$data$X[ is.na(model.input$data$X) ] <- -1000


input1 <- c(code=ms.ms.occ,
            model.input)

## *****************************************************************
## vanilla nimble/jags
## *****************************************************************

ms.ms.filter <- compareMCMCs_withMonitors(input1,
                                          MCMCs=c('nimble'),
                                          niter=niter,
                                          burnin = burnin,
                                          summary=FALSE,
                                          check=FALSE,
                                          monitors=model.input$monitors)

save(ms.ms.filter, file=file.path(save.dir, 'filtering.Rdata'))

