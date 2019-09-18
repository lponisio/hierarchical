library(RColorBrewer)
library(coda)

pdf.f <- function(f, file, ...) {
    cat(sprintf('Writing %s\n', file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
}

checkChains <- function(all.mods.samps, f.path,
                        only.one="jags",
                        params=NULL, prop.plot=0.1){
    ## function to plot values of parameters as a function of
    ## iterations. Input is a MCMCcompare object
    niter <- dim(all.mods.samps)[3]
    if(is.null(params)){
        params <- names(all.mods.samps[1,,1])
    }
    mods <- names(all.mods.samps[,1,1])
    if(is.null(mods)) mods <- only.one

    lapply(mods, function(z){
        f <- function(){
            layout(matrix(1:4, ncol=2))
            apply(all.mods.samps[z,params,], 1, function(x){
                plot(x[seq(from=1, to=length(x),
                           length.out=length(x)*prop.plot)], type="l",
                     xlab = 'iteration',
                     main= params[which(apply(all.mods.samps[z,,], 1,
                                              function(y)
                         all(match(x,y))))]
                     )
            })
        }

        pdf.f(f,
              file= file.path(sprintf(f.path, z)),
              height=11, width=8.5)
    })
}



plotPointsMakeTable <- function(occ.all, adj.xlab, sim.data=FALSE){
    layout(matrix(1:2, ncol=1), heights=c(1,3))
    par(oma=c(5, 6, 0.5, 1),
        mar=c(0.5, 0, 0.5, 1), cex.axis=1.5)
    sampler.names <-
        getSamplerNames(occ.all$MCMCresults$efficiency$min)$samplers
    unique.samplers <- unique(sampler.names)
    cols <- brewer.pal(length(unique.samplers), "Dark2")
    names(cols) <- unique.samplers
    params <- colnames(occ.all$MCMCresults$summary[,'efficiency',])

    latent <-
        sapply(strsplit(dimnames(occ.all$MCMCresults$summary)[[1]],
                        "latent"), function(x) x[2])
    latent <- sapply(strsplit(latent,'_'), function(x) x[1])
    pchs <- ifelse(latent, 16, 15)

    ltys=rep(1, length(latent))

    plot(NA, ylim=c(0,1), xlim=c(0,1), yaxt= "n", xaxt="n", bty="n")
    legend("top", legend=unique.samplers,
           cex=0.7, lty=1, lwd=2,
           col=cols, bty="n", ncol=length(unique.samplers))

    legend("bottomleft", legend=c("Latent states",
                                  "Latent state integration"),
           cex=0.7, pch=c(16,15), bty="n", ncol=2)

     if(sim.data){
        mus <-
            sapply(strsplit(dimnames(occ.all$MCMCresults$summary)[[1]],
                            "mup"), function(x) x[2])
        ltys <- ifelse(mus == 1, 1, 2)
        legend("bottomright",
               legend=c("High detectability", "Low detectability"),
           cex=0.7, lty=c(1,2), bty="n", ncol=2)
    }

    plot(NA,
         ylim=range(occ.all$MCMCresults$summary[,'efficiency',]),
         xlim=c(1, length(params)),  xaxt="n",
         ylab="", xlab="", las=2, cex.axis=0.8)

    text(x=1:length(params), par('usr')[3],
         srt = 45, adj=adj.xlab,
         labels = params,
         xpd = NA,
         cex=0.7)
    mtext("ESS/second",
          2, line=3.2, cex=1.2)


    for(i in 1:dim(occ.all$MCMCresults$summary)[1]){
        sum.output <- round(t(occ.all$MCMCresults$summary[i,,]),
                            digits=4)
        sum.output <- as.data.frame(sum.output)
        sum.output$param <- rownames(sum.output)
        rownames(sum.output) <- NULL
        sum.output$geweke <-
            geweke.diag(t(occ.all$MCMCresults$samples[i,,]))$z
        sum.output$sec1000samp <- 1000/sum.output$efficiency
        sum.output$min1000samp <- (1000/sum.output$efficiency)/60
        sum.output$hrs1000samp <- ((1000/sum.output$efficiency)/60)/60
        write.table(sum.output,
                    sep=",", row.names=FALSE,
                    file=file.path(save.dir,
                    sprintf("../tables/%s.csv",
                    dimnames(occ.all$MCMCresults$summary)[[1]][i])))
        points(x=1:length(params),
               y=occ.all$MCMCresults$summary[i,'efficiency',],
               col=cols[sampler.names[i]],
               pch=pchs[i],
               ltys[i],
               xaxt="n",
               type="o")
    }
}


makeCombinedTables <- function(pattern, save.dir, patternp=NULL){
    ## function for taking the summary tables and combining by scenario
    ## for supplemental table
    list.tables <- list.files(file.path(save.dir, "../tables"),
                              pattern=pattern)

    if(!is.null(patternp)){
        list.tables <- list.tables[grepl(patternp, list.tables)]
        pattern <- paste0(pattern, patternp)
        }
    tables <- lapply(file.path(save.dir, "../tables", list.tables),
                     read.csv)
    sum.tables <- lapply(tables, function(x){
        x[,c("mean", "sd", "geweke")]
    })

    sum.tables <- do.call(cbind, sum.tables)
    sum.tables <- round(sum.tables, digits=2)

    rownames(sum.tables) <- tables[[1]]$param
    cols <- colnames(sum.tables)
    sum.tables <- rbind(cols, sum.tables)

    ## split by whether latent states were sampled
    samp.latent <- grepl("latentTRUE", list.tables)
    nosamp.latent <- grepl("latentFALSE", list.tables)
    samp.names.prep <- sapply(strsplit(list.tables, "sampler"),
                              function(x) x[2])
    samp.names <-
        sapply(strsplit(samp.names.prep, ".csv"), function(x) x[1])
    colnames(sum.tables) <- rep(samp.names, each=3)

    ## split tables by whether latent states are sampled
    sum.tab.latent <- sum.tables[, rep(samp.latent, each=3)]
    sum.tab.nolatent <- sum.tables[, rep(nosamp.latent, each=3)]

    write.table(sum.tab.latent,
                sep="&",
                file=file.path(save.dir,
                               sprintf("../latentTRUE%s.txt",
                                       pattern)))
      write.table(sum.tab.nolatent,
                sep="&",
                file=file.path(save.dir,
                               sprintf("../latentFALSE%s.txt",
                                       pattern)))


}
