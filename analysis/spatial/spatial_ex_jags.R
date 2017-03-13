rm(list=ls())
library('R2jags')
library(RCurl)
load.module("msm")


## to download the data fro dropbox
dl_from_dropbox <- function(x, key) {
                        require(RCurl)
                        bin <- getBinaryURL(paste0("https://dl.dropboxusercontent.com/s/", key, "/", x),
                                            ssl.verifypeer = FALSE)
                        con <- file(x, open = "wb")
                        writeBin(bin, con)
                        close(con)
                        message(noquote(paste(x, "read into", getwd())))
                        }

 dl_from_dropbox("SimDat.Rdata", "0cj2mzgvyv0yfuo")
load('SimDat.Rdata')



## *********************************************************************
## jags model
## *********************************************************************

sink('model.jags')
cat('model{

  ## priors
  delta ~ dunif(0.1, 10)
  sigma ~ dunif(0.1, 100)
  p ~ dunif(0, 1)
  alpha ~ dnorm(0, 0.001)
  b1 ~ dnorm(0, 0.001)

  rho[1:nsite] ~ dmnorm(zeros[1:nsite],
                        D.tau[1:nsite, 1:nsite])


  ## Likelihood
  ## Ecological model for true occurrence
  for (i in 1:nsite) {
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- alpha + b1*elev[i] + rho[i]
    p.eff[i] <- z[i] * p

    ## Observation model for replicated detection/nondetection
    ## observations
    for (j in 1:nreps) {
      y[i,j] ~ dbern(p.eff[i])
    }
  }

  ## create covariance matrix based on distances (must be 1/cov for
  ## JAGS)

   prep.cov[1:nsite, 1:nsite] <- 1/(exp(delta)*mexp(D[1:nsite, 1:nsite]))
  D.cov[1:nsite, 1:nsite] <- (sigma^2)*(0.95*prep.cov[1:nsite, 1:nsite] + 0.05*DI[1:nsite, 1:nsite])

 ##   for(i in 1:nsite){
 ##   for(j in 1:nsite){
 ##     prep.cov[i, j]  <- exp(-delta*D[i, j])
 ##     D.cov[i, j] <- (sigma^2)*(0.95*prep.cov[i, j] + 0.05*DI[i, j])
 ##   }
 ## }


  D.tau[1:nsite, 1:nsite] <- inverse(D.cov[1:nsite, 1:nsite])

  }',fill = TRUE)
sink()


## *********************************************************************
## run the model
## *********************************************************************

analyse.jags <- function(d, ni, nt, nb, nc) {
    these.inits <- d$inits

    my.inits <- function() {
        these.inits
    }

    dd <- list(data=d$data, inits=my.inits, params=model.input$monitors)
    out.jags <- jags(d$data, my.inits, d$params,'model.jags', n.chains=nc,
                     n.thin=nt, n.iter=ni, n.burnin=nb, working.directory=NULL)
    reutrn(out.jags)
}


scale <- 1e3
res <- analyse.jags(model.input,
                    ni=(1e3+1e1)*scale,
                    nt=scale,
                    nb=1e1*scale,
                    nc=3)


