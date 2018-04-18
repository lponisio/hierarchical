ms.ms.occ <- nimbleCode({
    ## multi-species priors
    ## detectablility
    mu.p.0     ~ dnorm(0,0.01)
    mu.p.day.1 ~ dnorm(0,0.01)
    mu.p.day.2 ~ dnorm(0,0.01)
    sigma.p.0     ~ dunif(0,10)
    sigma.p.day.1 ~ dunif(0,10)
    sigma.p.day.2 ~ dunif(0,10)
    tau.p.0 <- 1/(sigma.p.0*sigma.p.0)
    tau.p.day.1 <- 1/(sigma.p.day.1*sigma.p.day.1)
    tau.p.day.2 <- 1/(sigma.p.day.2*sigma.p.day.2)

    ## phi/gam random intercepts
    mu.phi.0  ~ dnorm(0,0.01)
    mu.gam.0  ~ dnorm(0,0.01)
    sigma.phi.0 ~ dunif(0,10)
    sigma.gam.0 ~ dunif(0,10)
    tau.phi.0 <-  1/(sigma.phi.0*sigma.phi.0)
    tau.gam.0 <-  1/(sigma.gam.0*sigma.gam.0)

    ## ## hedgerow area proximity
    mu.phi.hr.area  ~ dnorm(0,0.01)
    mu.gam.hr.area  ~ dnorm(0,0.01)
    sigma.phi.hr.area ~ dunif(0,10)
    sigma.gam.hr.area ~ dunif(0,10)
    tau.phi.hr.area <- 1/(sigma.phi.hr.area*sigma.phi.hr.area)
    tau.gam.hr.area <- 1/(sigma.gam.hr.area*sigma.gam.hr.area)

    ## ## semi nat habitat area proximity
    mu.phi.nat.area  ~ dnorm(0,0.01)
    mu.gam.nat.area  ~ dnorm(0,0.01)
    sigma.phi.nat.area ~ dunif(0,10)
    sigma.gam.nat.area ~ dunif(0,10)
    tau.phi.nat.area <- 1/(sigma.phi.nat.area*sigma.phi.nat.area)
    tau.gam.nat.area <- 1/(sigma.gam.nat.area*sigma.gam.nat.area)

    ## ## floral resource diversity
    mu.phi.fra  ~ dnorm(0,0.01)
    mu.gam.fra  ~ dnorm(0,0.01)
    sigma.phi.fra ~ dunif(0,10)
    sigma.gam.fra ~ dunif(0,10)
    tau.phi.fra <- 1/(sigma.phi.fra*sigma.phi.fra)
    tau.gam.fra <- 1/(sigma.gam.fra*sigma.gam.fra)

    ## ## diet breadth
    phi.k  ~ dnorm(0,0.01)
    gam.k  ~ dnorm(0,0.01)

    ## ## body size
    phi.B  ~ dnorm(0,0.01)
    gam.B  ~ dnorm(0,0.01)

    ## ## interaction between hedgerow proximity and floral resources
    ## ## (habitat quality)
    phi.hr.area.fra  ~ dnorm(0,0.01)
    gam.hr.area.fra  ~ dnorm(0,0.01)

    ## ## interaction between semi nat habitat proximity and floral
    ## ## resources (habitat quality)
    phi.nat.area.fra  ~ dnorm(0,0.01)
    gam.nat.area.fra  ~ dnorm(0,0.01)

    ## ## interaction between hedgerow proximity and species diet breadth
    phi.hr.area.k  ~ dnorm(0,0.01)
    gam.hr.area.k  ~ dnorm(0,0.01)

    ## interaction between semi nat habitat proximity and species diet
    ## breadth
    phi.nat.area.k  ~ dnorm(0,0.01)
    gam.nat.area.k  ~ dnorm(0,0.01)

    ## interaction between hedgerow proximity and species body size
    phi.hr.area.B  ~ dnorm(0,0.01)
    gam.hr.area.B  ~ dnorm(0,0.01)

    ## interaction between semi nat habitat proximity and body size
    phi.nat.area.B  ~ dnorm(0,0.01)
    gam.nat.area.B  ~ dnorm(0,0.01)

    ## species-specific  parameters
    for(sp in 1:nsp) {
        ## day
        p.0[sp]     ~ dnorm(mu.p.0,     tau.p.0)
        p.day.1[sp] ~ dnorm(mu.p.day.1, tau.p.day.1)
        p.day.2[sp] ~ dnorm(mu.p.day.2, tau.p.day.2)

        ## random intercept
        phi.0[sp] ~ dnorm(mu.phi.0, tau.phi.0)
        gam.0[sp] ~ dnorm(mu.gam.0, tau.gam.0)

        ## ## hedgerow area
        phi.hr.area[sp] ~ dnorm(mu.phi.hr.area,
                                tau.phi.hr.area)
        gam.hr.area[sp] ~ dnorm(mu.gam.hr.area,
                                tau.gam.hr.area)

        ## ## natural habitat
        phi.nat.area[sp] ~ dnorm(mu.phi.nat.area,
                                 tau.phi.nat.area)
        gam.nat.area[sp] ~ dnorm(mu.gam.nat.area,
                                 tau.gam.nat.area)

        ## ## fra
        phi.fra[sp] ~ dnorm(mu.phi.fra,
                            tau.phi.fra)
        gam.fra[sp] ~ dnorm(mu.gam.fra,
                            tau.gam.fra)

        for(site in 1:nsite) {
            for(yr in 1:nyear) {
                for(rep in 1:nrep[site,yr,sp]) {
                    logit(p[site,yr,rep,sp]) <-
                        p.0[sp] +
                        p.day.1[sp]*day[site,yr,rep,sp] +
                        p.day.2[sp]*day.2[site,yr,rep,sp]
                }
            }
        }
    }

    for(sp in 1:nsp) {
        for(site in 1:nsite) {
            ## start off at the average for each species, site across years
            logit(phi.site.sp.mean[site,sp]) <- mean(phi[site, 1:(nyear-1),sp])
            logit(gam.site.sp.mean[site,sp]) <- mean(gam[site, 1:(nyear-1),sp])

            psi.1[site,sp] <- gam.site.sp.mean[site,sp]/
                (1 - phi.site.sp.mean[site,sp] + gam.site.sp.mean[site,sp])

            ## occupancy in year 1
            psi[site,1,sp] <- psi.1[site,sp]
            Z[site,1,sp] ~ dbern(psi[site,1,sp])

            ## detectability in year 1
            for(rep in 1:nrep[site,1,sp]) {
                mu.p[site,1,rep,sp] <- Z[site,1,sp]*p[site,1,rep,sp]
                X[site,1,rep,sp] ~ dbern(mu.p[site,1,rep,sp])
            }

            ## occupancy in subsequent years
            for(yr in 1:(nyear-1)) {
                phi[site,yr,sp] <-
                    phi.0[sp] +
                    phi.k*k[sp] +
                    phi.B*B[sp] +
                    phi.hr.area[sp]*HRarea[site] +
                    phi.nat.area[sp]*natural[site]  +
                    phi.fra[sp]*fra[site, yr] +
                    phi.hr.area.fra*fra[site, yr]*HRarea[site]+
                    phi.nat.area.fra*fra[site, yr]*natural[site] +
                    phi.hr.area.k*k[sp]*HRarea[site] +
                    phi.nat.area.k*k[sp]*natural[site] +
                    phi.hr.area.B*B[sp]*HRarea[site] +
                    phi.nat.area.B*B[sp]*natural[site]

                gam[site,yr,sp] <-
                    gam.0[sp] +
                    gam.k*k[sp] +
                    gam.B*B[sp] +
                    gam.hr.area[sp]*HRarea[site] +
                    gam.nat.area[sp]*natural[site]  +
                    gam.fra[sp]*fra[site, yr]  +
                    gam.hr.area.fra*fra[site, yr]*HRarea[site] +
                    gam.nat.area.fra*fra[site, yr]*natural[site]+
                    gam.hr.area.k*k[sp]*HRarea[site] +
                    gam.nat.area.k*k[sp]*natural[site] +
                    gam.hr.area.B*B[sp]*HRarea[site] +
                    gam.nat.area.B*B[sp]*natural[site]

                logit(psi[site,yr+1,sp]) <-
                    Z[site,yr,sp] * phi[site,yr,sp] +
                    (1-Z[site,yr,sp]) * gam[site,yr,sp]

                Z[site,yr+1,sp] ~ dbern(psi[site,yr+1,sp])

                ## detectability in != year 1
                for(rep in 1:nrep[site,yr+1,sp]) {
                    mu.p[site,yr+1,rep,sp] <- Z[site,yr+1,sp]*p[site,yr+1,rep,sp]
                    X[site,yr+1,rep,sp] ~ dbern(mu.p[site,yr+1,rep,sp])
                }
            }
        }
    }


})

