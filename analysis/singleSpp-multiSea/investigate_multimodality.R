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

mcmc <- buildMCMC(m)

cm <- compileNimble(m)
cmcmc <- compileNimble(mcmc, project = m)

cmcmc$run(100000)
samples <- as.matrix(cmcmc$mvSamples)
colnames(samples)
plot(samples[,'psi1'])
plot(samples[,'mu.phi.mean'])
plot(samples[,'mu.p.mean'])
plot(samples[,'mu.gamma.mean'])
plot(samples[,'sigma.p'])

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

