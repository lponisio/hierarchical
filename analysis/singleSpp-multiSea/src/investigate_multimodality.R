## This file contains work to investigate why the MCMC
## in some cases ends up in an incorrect posterior region.

rm(list=ls())
source('src/initialize.R')

mc <- makeModel(latent = TRUE, hyper.param = TRUE)

# set.seed(123) ## not as bad
set.seed(1234)
mu.p <- 1
data <- genDynamicOccData(mu.p=mu.p,
                          psi1=psi1,
                          sigma.p=sigma.p,
                          mu.phi=mu.phi,
                          sigma.phi=sigma.phi,
                          mu.gamma=mu.gamma,
                          sigma.gamma=sigma.gamma)

model.input <- prepModDataOcc(data, include.zs=TRUE)

## presentButNeverSeen <- data$z & (apply(data$y, c(1, 3), sum) == 0)
## There are no false absences!

## All z's were 1 in data or inits.
## I will introduce some zeros to increase chance of good mixing.

newZinits <- model.input$inits$z
newZinits[newZinits==1] <- 0
model.input$inits$z <- newZinits

m <- nimbleModel(
    mc,
    data = model.input$data,
    inits = model.input$inits,
    constants = model.input$constants
)

mcmc <- buildMCMC(m,
                  monitors2 = c('psi1',
                                'mu.p.mean','mu.p','sigma.p',
                                'mu.phi.mean','mu.phi',
                                'mu.gamma.mean','mu.gamma',
                                'sigma.phi','sigma.gamma',
                                'p','phi','gamma',
                                'z'),
                  thin2 = 10)

cm <- compileNimble(m)
cmcmc <- compileNimble(mcmc, project = m)

cmcmc$run(10000)
samples <- as.matrix(cmcmc$mvSamples)
samples2 <- as.matrix(cmcmc$mvSamples2)
colnames(samples)
colnames(samples2)
plot(samples[,'psi1'])
plot(samples[,'mu.phi.mean'])
plot(samples[,'mu.p.mean'])
plot(samples[,'mu.gamma.mean'])
plot(samples[,'sigma.p'])

## This prints the number of 1s in the state sequences.
## Almost all 1s or almost no 1s indicates almost
## no state flipping in the sampling.
for(i in 1:15) {
    matchstring <- paste0('z\\[[[:digit:]]{1,2}\\, ',i,']')
    zibool <- grepl(matchstring, colnames(samples2))
    if(sum(zibool) != 50) stop(paste0('problem for ', i))
    print(i)
    print(apply(samples2[, zibool], 2, sum))
}

## Model without latent states is below.

mcNL <- makeModel(latent = FALSE, hyper.param = TRUE)
mNL <- nimbleModel(
    mcNL,
    data = model.input$data,
    inits = model.input$inits,
    constants = model.input$constants
)

mcmcNL <- buildMCMC(mNL)

cmNL <- compileNimble(mNL)
cmcmcNL <- compileNimble(mcmcNL, project = mNL)

cmcmcNL$run(10000)
samplesNL <- as.matrix(cmcmcNL$mvSamples)
colnames(samplesNL)
plot(samplesNL[,'psi1'])
plot(samplesNL[,'mu.phi.mean'])
plot(samplesNL[,'mu.p.mean'])
plot(samplesNL[,'mu.gamma.mean'])
plot(samplesNL[,'sigma.p'])
plot(samplesNL[,'sigma.gamma'])

## 1. compare probabilities integrated over the z's from the latent-states model
cmNL$calculate()
cmVars <- cm$getVarNames()
cmNLvars <- cmNL$getVarNames()
for(v in cmVars) {
    if(v %in% cmNLvars)
        cmNL[[v]] <- cm[[v]]
}


