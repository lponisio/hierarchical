
makeModel <- function(latent, hyper.param){
    if(latent){
        if(hyper.param){
            ## latent states and hyper param
            ms.ss.occ <- nimbleCode({
                ## Define prior distributions for community-level model
                ## parameters
                cato.occ.mean ~ dunif(0,1)
                mu.ucato <- log(cato.occ.mean) - log(1-cato.occ.mean)
                sigma.ucato ~ dunif(0, 100)
                tau.ucato <-  1/(sigma.ucato*sigma.ucato)

                fcw.occ.mean ~ dunif(0,1)
                mu.ufcw <- log(fcw.occ.mean) - log(1-fcw.occ.mean)
                sigma.ufcw ~ dunif(0, 100)
                tau.ufcw <-  1/(sigma.ufcw*sigma.ufcw)

                cato.det.mean ~ dunif(0,1)
                mu.vcato <- log(cato.det.mean) - log(1-cato.det.mean)
                sigma.vcato ~ dunif(0, 100)
                tau.vcato <-  1/(sigma.vcato*sigma.vcato)

                fcw.det.mean ~ dunif(0,1)
                mu.vfcw <- log(fcw.det.mean) - log(1-fcw.det.mean)
                sigma.vfcw ~ dunif(0, 100)
                tau.vfcw <-  1/(sigma.vfcw*sigma.vfcw)

                ## random effects
                mu.a1 ~ dnorm(0, 0.001)
                sigma.a1 ~ dunif(0, 100)
                tau.a1 <-  1/(sigma.a1*sigma.a1)
                mu.a2 ~ dnorm(0, 0.001)
                sigma.a2 ~ dunif(0, 100)
                tau.a2 <-  1/(sigma.a2*sigma.a2)
                mu.a3 ~ dnorm(0, 0.001)
                sigma.a3 ~ dunif(0, 100)
                tau.a3 <-  1/(sigma.a3*sigma.a3)
                mu.a4 ~ dnorm(0, 0.001)
                sigma.a4 ~ dunif(0, 100)
                tau.a4 <-  1/(sigma.a4*sigma.a4)
                mu.b1 ~ dnorm(0, 0.001)
                sigma.b1 ~ dunif(0, 100)
                tau.b1 <-  1/(sigma.b1*sigma.b1)
                mu.b2 ~ dnorm(0, 0.001)
                sigma.b2 ~ dunif(0, 100)
                tau.b2 <-  1/(sigma.b2*sigma.b2)

                for (i in 1:(num.species)) {
                    ## Create priors for species i from the community level prior
                    ## distributions
                    u.cato[i] ~ dnorm(mu.ucato, tau.ucato)
                    u.fcw[i] ~ dnorm(mu.ufcw, tau.ufcw)
                    a1[i] ~ dnorm(mu.a1, tau.a1)
                    a2[i] ~ dnorm(mu.a2, tau.a2)
                    a3[i] ~ dnorm(mu.a3, tau.a3)
                    a4[i] ~ dnorm(mu.a4, tau.a4)
                    v.cato[i] ~ dnorm(mu.vcato, tau.vcato)
                    v.fcw[i] ~ dnorm(mu.vfcw, tau.vfcw)
                    b1[i] ~ dnorm(mu.b1, tau.b1)
                    b2[i] ~ dnorm(mu.b2, tau.b2)
                    ## Create a loop to estimate the Z matrix (true occurrence for
                    ## species i at point j).

                    for (j in 1:num.points) {
                        logit(psi[j,i]) <- u.cato[i]*(1-habitat.ind[j]) +
                            u.fcw[i]*habitat.ind[j] +
                            a1[i]*ufc.linear[j] +
                            a2[i]*ufc.quadratic[j] +
                            a3[i]*ba.linear[j] +
                            a4[i]*ba.quadratic[j]
                        mu.psi[j,i] <- psi[j,i]
                        Z[j,i] ~ dbern(mu.psi[j,i])
                        ## Create a loop to estimate detection for species i at point k
                        ## during sampling period k.
                        for (k in 1:num.reps[j]) {
                            logit(p[j,k,i]) <-  v.cato[i]*(1-habitat.ind[j]) +
                                v.fcw[i]*habitat.ind[j] +
                                b1[i]*date.linear[j,k] +
                                b2[i]*date.quadratic[j,k]
                            mu.p[j,k,i] <- p[j,k,i]*Z[j,i]
                            X[j,k,i] ~ dbern(mu.p[j,k,i])
                        }
                    }
                }
            })




        } else {
            ## latent states, no hypper param
            ms.ss.occ <- nimbleCode({
                ## Define prior distributions
                a1 ~ dnorm(0, 0.001)
                a2 ~ dnorm(0, 0.001)
                a3 ~ dnorm(0, 0.001)
                a4 ~ dnorm(0, 0.001)
                b1 ~ dnorm(0, 0.001)
                b2 ~ dnorm(0, 0.001)

                cato.occ ~ dunif(0,1)
                u.cato <- log(cato.occ) - log(1-cato.occ)

                fcw.occ ~ dunif(0,1)
                u.fcw <- log(fcw.occ) - log(1-fcw.occ)

                cato.det ~ dunif(0,1)
                v.cato <- log(cato.det) - log(1-cato.det)

                fcw.det ~ dunif(0,1)
                v.fcw <- log(fcw.det) - log(1-fcw.det)


                for (i in 1:(num.species)) {
                    ## Create a loop to estimate the Z matrix (true occurrence for
                    ## species i at point j).

                    for (j in 1:num.points) {
                        logit(psi[j,i]) <- u.cato*(1-habitat.ind[j]) +
                            u.fcw*habitat.ind[j] +
                            a1*ufc.linear[j] +
                            a2*ufc.quadratic[j] +
                            a3*ba.linear[j] +
                            a4*ba.quadratic[j]
                        mu.psi[j,i] <- psi[j,i]
                        Z[j,i] ~ dbern(mu.psi[j,i])
                        ## Create a loop to estimate detection for species i at point k
                        ## during sampling period k.
                        for (k in 1:num.reps[j]) {
                            logit(p[j,k,i]) <-  v.cato*(1-habitat.ind[j]) +
                                v.fcw*habitat.ind[j] +
                                b1*date.linear[j,k] +
                                b2*date.quadratic[j,k]
                            mu.p[j,k,i] <- p[j,k,i]*Z[j,i]
                            X[j,k,i] ~ dbern(mu.p[j,k,i])
                        }
                    }
                }
            })
        }

    } else if (!latent){
        if(hyper.param){
            ## integrate over latent states and with  hyper param
            ms.ss.occ <- nimbleCode({
                ## Define prior distributions for community-level model parameters
                cato.occ.mean ~ dunif(0,1)
                mu.ucato <- log(cato.occ.mean) - log(1-cato.occ.mean)
                sigma.ucato ~ dunif(0, 100)

                fcw.occ.mean ~ dunif(0,1)
                mu.ufcw <- log(fcw.occ.mean) - log(1-fcw.occ.mean)
                sigma.ufcw ~ dunif(0, 100)

                cato.det.mean ~ dunif(0,1)
                mu.vcato <- log(cato.det.mean) - log(1-cato.det.mean)
                sigma.vcato ~ dunif(0, 100)

                fcw.det.mean ~ dunif(0,1)
                mu.vfcw <- log(fcw.det.mean) - log(1-fcw.det.mean)
                sigma.vfcw ~ dunif(0, 100)

                ## random effects
                mu.a1 ~ dnorm(0, 0.001)
                sigma.a1 ~ dunif(0, 100)
                mu.a2 ~ dnorm(0, 0.001)
                sigma.a2 ~ dunif(0, 100)
                mu.a3 ~ dnorm(0, 0.001)
                sigma.a3 ~ dunif(0, 100)
                mu.a4 ~ dnorm(0, 0.001)
                sigma.a4 ~ dunif(0, 100)
                mu.b1 ~ dnorm(0, 0.001)
                sigma.b1 ~ dunif(0, 100)
                mu.b2 ~ dnorm(0, 0.001)
                sigma.b2 ~ dunif(0, 100)


                for (i in 1:(num.species)) {
                    ## Create priors for species i from the community level prior
                    ## distributions

                    u.cato[i] ~ dnorm(mu.ucato, sd=sigma.ucato)
                    u.fcw[i] ~ dnorm(mu.ufcw, sd=sigma.ufcw)
                    a1[i] ~ dnorm(mu.a1, sd=sigma.a1)
                    a2[i] ~ dnorm(mu.a2, sd=sigma.a2)
                    a3[i] ~ dnorm(mu.a3, sd=sigma.a3)
                    a4[i] ~ dnorm(mu.a4, sd=sigma.a4)

                    v.cato[i] ~ dnorm(mu.vcato, sd=sigma.vcato)
                    v.fcw[i] ~ dnorm(mu.vfcw, sd=sigma.vfcw)
                    b1[i] ~ dnorm(mu.b1, sd=sigma.b1)
                    b2[i] ~ dnorm(mu.b2, sd=sigma.b2)

                    ## vectorize the calculation of psi.
                    logit(psi[1:num.points,i]) <-
                        u.cato[i]*(1-habitat.ind[1:num.points]) +
                        u.fcw[i]*habitat.ind[1:num.points] +
                        a1[i]*ufc.linear[1:num.points] +
                        a2[i]*ufc.quadratic[1:num.points] +
                        a3[i]*ba.linear[1:num.points] +
                        a4[i]*ba.quadratic[1:num.points]
                    mu.psi[1:num.points,i] <- psi[1:num.points, i]

                    logit(p[1:num.points, 1:max.num.reps, i]) <-
                        (v.cato[i]*(1-habitat.ind[1:num.points]) +
                         v.fcw[i]*habitat.ind[1:num.points]) %*%
                        asRow(onesRow[1, 1:max.num.reps])+
                        b1[i]*date.linear[1:num.points,1:max.num.reps] +
                        b2[i]*date.quadratic[1:num.points,1:max.num.reps]

                    ## user defined distribution to combine the bernoulli occupancy
                    ## and detection events.  We can also make this is a single
                    ## compuation for the entire matrix of locations-x-visits, for
                    ## each species (i)

                    X[1:num.points, 1:max.num.reps, i] ~ dBernDetectionMatrix(
                        occProb = mu.psi[1:num.points,i],
                        detectionProb = p[1:num.points, 1:max.num.reps,i],
                        numReps = num.reps[1:num.points])
                }
            })

        } else{
            ## integrate over latent states and remove hyper param
            ms.ss.occ <- nimbleCode({
                ## random effects
                a1 ~ dnorm(0, 0.001)
                a2 ~ dnorm(0, 0.001)
                a3 ~ dnorm(0, 0.001)
                a4 ~ dnorm(0, 0.001)
                b1 ~ dnorm(0, 0.001)
                b2 ~ dnorm(0, 0.001)
                cato.occ ~ dunif(0,1)
                u.cato <- log(cato.occ) - log(1-cato.occ)
                fcw.occ ~ dunif(0,1)
                u.fcw <- log(fcw.occ) - log(1-fcw.occ)
                cato.det ~ dunif(0,1)
                v.cato <- log(cato.det) - log(1-cato.det)
                fcw.det ~ dunif(0,1)
                v.fcw <- log(fcw.det) - log(1-fcw.det)

                for (i in 1:(num.species)) {

                    ## vectorize the calculation of psi.
                    logit(psi[1:num.points,i]) <-
                        u.cato*(1-habitat.ind[1:num.points]) +
                        u.fcw*habitat.ind[1:num.points] +
                        a1*ufc.linear[1:num.points] +
                        a2*ufc.quadratic[1:num.points] +
                        a3*ba.linear[1:num.points] +
                        a4*ba.quadratic[1:num.points]
                    mu.psi[1:num.points,i] <- psi[1:num.points, i]

                    logit(p[1:num.points, 1:max.num.reps, i]) <-
                        (v.cato*(1-habitat.ind[1:num.points]) +
                         v.fcw*habitat.ind[1:num.points]) %*%
                        asRow(onesRow[1, 1:max.num.reps])+
                        b1*date.linear[1:num.points,1:max.num.reps] +
                        b2*date.quadratic[1:num.points,1:max.num.reps]

                    ## user defined distribution to combine the bernoulli occupancy
                    ## and detection events.  We can also make this is a single
                    ## compuation for the entire matrix of locations-x-visits, for
                    ## each species (i)

                    X[1:num.points, 1:max.num.reps, i] ~ dBernDetectionMatrix(
                        occProb = mu.psi[1:num.points,i],
                        detectionProb = p[1:num.points, 1:max.num.reps,i],
                        numReps = num.reps[1:num.points])
                }
            })

        }
    }
}
