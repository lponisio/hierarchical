
library(nimble)
setwd('~/github/occupancy/analysis/multiSpp-singleSea')
source('src/initialize.R')
source('DT_code_latent.R')
source('DT_code_filter.R')
samplesList <- vector('list', 3)
summary <- vector('list', 3)

case <- 1   ## filtering
case <- 2   ## latent states
case <- 3   ## latent states w/ crossLevel

cases <- 1:3

for(case in cases) {
    if(case == 1) { filter <- TRUE }
    if(case == 2) { filter <- FALSE; CL <- FALSE }
    if(case == 3) { filter <- FALSE; CL <- TRUE }
    set.seed(0)
    n.zeroes <- 0
    model.input <- prepMutiSpData(survey.data, survey.dates, species.groups, habitat, n.zeros, monitors, remove.zs=filter)
    ##
    ## CHANGE: get correct initial values from filtering results
    load(file=file.path(save.dir, 'filter.Rdata'))
    filter_inits <- ms.ss.filter[[1]]$summary['nimble', 'median',]
    e <- new.env()
    for(n in names(filter_inits))    eval(parse(text = paste0('e$', n, ' <- filter_inits[\'', n, '\']'))[[1]])
    new_inits <- list()
    for(n in ls(e)) new_inits[[n]] <- e[[n]]
    model.input$inits <- c(model.input$inits, new_inits)
    ##model.input$inits <- c(model.input$inits,  ms.ss.filter[[1]]$summary['nimble', 'median',])   ## OLD
    rm(ms.ss.filter)
    ##
    code <- if(filter) ms.ss.occ.filter else ms.ss.occ.latent
    constants <- model.input$constants
    data <- model.input$data
    inits <- model.input$inits
    ##
    Rmodel <- nimbleModel(code, constants, data, inits)
    if(filter) cat('filtering') else cat('latent')
    if(!filter && CL) cat('with CL')
    cat('\n')
    cat('calculate(Rmodel): '); cat(calculate(Rmodel)); cat('\n')
    ## missing inits: psi, mu.psi ???   (these ones ok: p, mu.p, X)
    ##
    conf <- configureMCMC(Rmodel)
    if(!filter && CL) {
        base.names <- c('a1', 'a2', 'a3', 'a4', 'b1', 'b2', 'u.cato', 'u.fcw', 'v.cato', 'v.fcw')
        conf$removeSamplers('Z')
        conf$removeSamplers(base.names)
        ##conf$addSampler(target = base.names, type ='sampler_crossLevel_binary_DT')
        Zparents <- c('a1', 'a2', 'a3', 'a4', 'u.cato', 'u.fcw')
        Xparents <- c('b1', 'b2', 'v.cato', 'v.fcw')
        for(i in 1:58)
            conf$addSampler(target = paste0(Zparents,'[',i,']'),
                            type ='sampler_crossLevel_binary_DT')
        for(i in 1:58)
            conf$addSampler(target = paste0(Xparents,'[',i,']'),
                            type ='sampler_crossLevel_binary_DT')
    }
    ##conf$printSamplers()
    conf$addMonitors(model.input$monitors, print = FALSE)
    Rmcmc <- buildMCMC(conf)
    ##
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    ##
    set.seed(0)
    ##niter <- 2000
    niter <- 20000
    samples <- runMCMC(Cmcmc, niter)
    ##
    ##colnames(samples)
    paramName <- 'cato.occ.mean'
    print(apply(samples, 2, mean)[paramName])
    samplesPlot(samples, paramName)
    samplesList[[case]] <- samples
    summary[[case]] <- samplesSummary(samples)
}

save(samplesList, summary, file = 'DT_results.RData')

f <- function(param) {
    l <- lapply(samplesList[1:2], function(x) x[,param])
    samp <- do.call('cbind', l)
    colnames(samp) <- c('filter', 'latent')##, 'latentCL')
    samplesPlot(samp)
}

f('cato.occ.mean')
f('u.cato[1]')
f('u.cato[5]')

names <- colnames(samplesList[[1]])
nplots <- 30
for(i in 1:nplots) {
    n <- sample(names, 1)
    f(n)
}



base.names <- c("a1", "a2", "a3", "a4", "b1", "b2", "u.cato", "u.fcw", "v.cato", "v.fcw" )
Zparents <- c('a1', 'a2', 'a3', 'a4', 'u.cato', 'u.fcw')
Xparents <- c('b1', 'b2', 'v.cato', 'v.fcw')
loopnames <- Zparents
loopnames <- Xparents

for(n in loopnames) {
    print(n)
    ##print(Rmodel$expandNodeNames(n))
    print(paste0(n, '[4]'))
    print(Rmodel$getDependencies(paste0(n, '[4]'), stochOnly=TRUE, includeData=FALSE))
}


Zparents <- c('a1', 'a2', 'a3', 'a4', 'u.cato', 'u.fcw')
Xparents <- c('b1', 'b2', 'v.cato', 'v.fcw')
i <- 4
paste0(Zparents, '[', i, ']')
paste0(Xparents, '[', i, ']')



conf <- configureMCMC(Rmodel)
conf$printSamplers()

base.names <- c('a1', 'a2', 'a3', 'a4', 'b1', 'b2', 'u.cato', 'u.fcw', 'v.cato', 'v.fcw')
conf$removeSamplers('Z')
conf$removeSamplers(base.names)
conf$printSamplers()

Zparents <- c('a1', 'a2', 'a3', 'a4', 'u.cato', 'u.fcw')
Xparents <- c('b1', 'b2', 'v.cato', 'v.fcw')
for(i in 1:58) {
    conf$addSampler(target = paste0(Zparents,'[',i,']'),
                          type ='sampler_crossLevel_binary_DT')
}
for(i in 1:58) {
    conf$addSampler(target = paste0(Xparents,'[',i,']'),
                          type ='sampler_crossLevel_binary_DT')
}

conf$printSamplers()


##filtering
##calculate(Rmodel): -3636.964
##cato.occ.mean 
##    0.4207882 
## 
##latent
##calculate(Rmodel): -7721.592
##cato.occ.mean 
##     0.405873 
## 
##latentwith CL
##calculate(Rmodel): -7721.592
##cato.occ.mean 
##    0.3439535 


ms.ss.occ.all[[1]]$efficiency




