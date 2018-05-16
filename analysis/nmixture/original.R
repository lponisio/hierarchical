rm(list=ls())
setwd("~/Dropbox/occupancy")
setwd("analysis/nmixture")
source('src/initialize.R')
source("src/setup.R")


## This example is adapted for NIMBLE from the AHM book by Jacob Levine and Perry de Valpine
## 6.11.1 Bayesian fitting of the basic ZIP N-mixture model
## ------------------------------------------------------------------------
## From section 6.9.2 - for data creation/organization:

## Specify model in BUGS language
## "ZIPNmix.txt"
nmixture <- nimbleCode( {

    ## Specify priors
    ## zero-inflation/suitability
    phi ~ dunif(0,1)          ## proportion of suitable sites (probability of being not a structural 0)
    theta <- 1-phi            ## zero-inflation (proportion of unsuitable)
    ltheta <- logit(theta)

    ## abundance
    beta0 ~ dnorm(0, 0.1)     ## log(lambda) intercept
    for(k in 1:7){            ## Regression params in lambda
        beta[k] ~ dnorm(0, 1)
    }
    tau.lam <- pow(sd.lam, -2)
    sd.lam ~ dunif(0, 2)      ## site heterogeneity in lambda

    ## detection
    for(j in 1:3){
        alpha0[j] <- logit(mean.p[j])
        mean.p[j] ~ dunif(0, 1) ## p intercept for occasions 1-3
    }
    for(k in 1:13){           ## Regression params in p
        alpha[k] ~ dnorm(0, 1)
    }
    tau.p.site <- pow(sd.p.site, -2)
    sd.p.site ~ dunif(0, 2)   ## site heterogeneity in p
    tau.p.survey <- pow(sd.p.survey, -2)
    sd.p.survey ~ dunif(0, 2) ## site-survey heterogeneity in p

    ## ZIP model for abundance
    for (i in 1:nsite){
        a[i] ~ dbern(phi)
        eps.lam[i] ~ dnorm(0, tau.lam)       ## Random site effects in log(abundance)
        loglam[i] <- beta0 + inprod(beta[1:7], lamDM[i, 1:7]) + eps.lam[i]
        loglam.lim[i] <- min(250, max(-250, loglam[i]))  ## Stabilize log
        lam[i] <- exp(loglam.lim[i])
        mu.poisson[i] <- a[i] * lam[i]
        N[i] ~ dpois(mu.poisson[i])
    }

    ## Measurement error model
    for (i in 1:nsite){
        eps.p.site[i] ~ dnorm(0, tau.p.site) ## Random site effects in logit(p)
        for (j in 1:nrep){
            y[i,j] ~ dbin(p[i,j], N[i])
            p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
            lp.lim[i,j] <- min(250, max(-250, lp[i,j]))  ## Stabilize logit
            lp[i,j] <- alpha0[j] + alpha[1] * elev[i] + alpha[2] * elev2[i] +
                alpha[3] * date[i,j] + alpha[4] * date2[i,j] +
                alpha[5] * dur[i,j] + alpha[6] * dur2[i,j] +
                alpha[7] * elev[i] * date[i,j] + alpha[8] * elev2[i] * date[i,j] +
                alpha[9] * elev[i] * dur[i,j] + alpha[10] * elev[i] * dur2[i,j] +
                alpha[11] * elev2[i] * dur[i,j] + alpha[12] * date[i,j] * dur[i,j] +
                alpha[13] * date[i,j] * dur2[i,j] +
                eps.p.site[i] + eps.p.survey[i,j]
            eps.p.survey[i,j] ~ dnorm(0, tau.p.survey) ## Random site-survey effects
        }
    }
}
)



input1 <- c(code=nmixture,
            model.input)


## *********************************************************************
## original: vanilla nimble and JAGS
## *********************************************************************

nmixture.orig <- compareMCMCs(input1,
                           MCMCs=c('jags', 'nimble'),
                           niter=niter,
                           burnin = burnin,
                           summary=FALSE,
                           check=FALSE)



save(nmixture.orig, file=file.path(save.dir, "orig.Rdata"))


## subset data for plotting
nsite <- 1
nrep <- 2

model.input.plot <- model.input

model.input.plot$constants <- model.input$constants
model.input.plot$constants$nsite <- nsite
model.input.plot$constants$nrep <- nrep

model.input.plot$data[sapply(model.input$data, is.matrix)] <-
    lapply(model.input$data[sapply(model.input$data, is.matrix)],
           function(x) x[1:nsite, ])


model.input.plot$data[c("y", "date", "dur", "date2", "dur2")] <-
    lapply(model.input.plot$data[c("y", "date", "dur", "date2", "dur2")], function(y)y[,1:nrep])

model.input.plot$data[sapply(model.input$data, is.vector)] <-
    lapply(model.input$data[sapply(model.input$data, is.vector)],
           function(x) x[1:nsite])

nmix <-   nimbleModel(code=nmixture,
                        constants=model.input.plot$constants,
                        data=model.input.plot$data,
                        inits=model.input.plot$inits,
                        check=FALSE)


## tkid <- tkplot(nmix$graph, vertex.shape="none", vertex.label.font=2.5,

##       vertex.label.cex=.6, edge.color="gray70",  edge.width=2, layout=layout_as_tree)

## byhand <- tkplot.getcoords(tkid)


plot(nmix$graph, vertex.shape="none", vertex.label.font=2.5,

      vertex.label.cex=.6, edge.color="gray70",  edge.width=2, layout=layout_in_circle)



## full data
nmix.full <-   nimbleModel(code=nmixture,
                        constants=model.input$constants,
                        data=model.input$data,
                        inits=model.input$inits,
                        check=FALSE)

deg <- degree(nmix.full$graph, mode="all")
hist(log(deg), main="Histogram of node degree")

max(deg)
mean(deg)
median(deg)
