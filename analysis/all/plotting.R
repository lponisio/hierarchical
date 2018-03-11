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
                plot(x[seq(from=1, to=length(x), length.out=length(x)*prop.plot)], type="l",
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


plotEffSize <- function(eff.size, eff.param,
                        f.path, name, at, adj1, adj2,
                        widths=c(4.5, 6)){
    cols <- brewer.pal(length(eff.size$mean)+1, "Greys")[-1]
    names(cols) <- names(eff.size$mean)
    f <- function(){
        layout(matrix(1:2, ncol=1))
        par(oma=c(6.5, 8, 0.5, 1),
            mar=c(0.5, 0, 2.5, 1), cex.axis=1.5)
        ## barplots
        mp1 <- barplot(eff.size$mean, names="", las=1, col=cols)
        mtext("Mean", 3, line=0.5, cex=1.2)

        mp2 <- barplot(eff.size$min, names="", las=1, col=cols)
        text(mp2, par('usr')[3] - adj1,
             srt = 45, adj = 1,
             labels = names(eff.size$mean),
             xpd = NA,
             cex=1)

        mtext("Minimum", 3, line=0.5, cex=1.2)
        mtext("Effective sample size \n per second",
              2, line=4.5, cex=1.5, at=at)
    }

    f2 <- function(){
        layout(matrix(1:2, ncol=1), heights=c(1,3))

        par(oma=c(5, 6, 0.5, 1),
            mar=c(0.5, 0, 0.5, 1), cex.axis=1.5)

        plot(NA, ylim=c(0,1), xlim=c(0,1), yaxt= "n", xaxt="n", bty="n")
        legend("top", legend=rownames(eff.param),
               pch=16, col=cols, bty="n", ncol=3, cex=0.7)

        plot(NA, ylim=log(range(eff.param)), xlim=c(1, ncol(eff.param)),
             xlab="", ylab="", xaxt="n")
        for(i in 1:nrow(eff.param)){
            points(x=1:ncol(eff.param), y=log(eff.param[i,]),
                   col=cols[i],
                   pch=16,
                   xaxt="n",
                   type="o")
        }
        text(x=1:ncol(eff.param), par('usr')[3],
             srt = 45, adj = 1 + adj2,
             labels = colnames(eff.param),
             xpd = NA,
             cex=0.8)
        mtext("Effective sample size \n per second (log)",
              2, line=3.2, cex=1.5)
    }

    f3 <- function(){
        layout(matrix(1:2, ncol=1), heights=c(1,3))

        par(oma=c(5, 6, 0.5, 1),
            mar=c(0.5, 0, 0.5, 1), cex.axis=1.5)

        diffs <- t(apply(eff.param, 1,
                         function(x) log(x) - log(eff.param["JAGS-latent",])))[-1,]

        plot(NA, ylim=c(0,1), xlim=c(0,1), yaxt= "n", xaxt="n", bty="n")
        legend("top", legend=rownames(diffs),
               cex=0.7, pch=16, col=cols[-1], bty="n", ncol=3)

        plot(NA, ylim=range(diffs), xlim=c(1, ncol(diffs)),
             xlab="", ylab="", xaxt="n")
        for(i in 1:nrow(diffs)){
            points(x=1:ncol(diffs), y=diffs[i,],
                   col=cols[i+1],
                   pch=16,
                   xaxt="n",
                   type="o")
        }
        text(x=1:ncol(diffs), par('usr')[3],
             srt = 45, adj = 1 + adj2,
             labels = colnames(diffs),
             xpd = NA,
             cex=0.8)
        mtext("Ratio of effective sample size \n per second (log)",
              2, line=3.2, cex=1.5)
        abline(h=0, lty=2)
    }
    pdf.f(f,
          file= file.path(sprintf(f.path, name, "Bar")),
          height=6, width=widths[1])
    pdf.f(f2,
          file= file.path(sprintf(f.path, name, "Points")),
          height=4, width=widths[2])
    pdf.f(f3,
          file= file.path(sprintf(f.path, name, "Diffs")),
          height=4, width=widths[2])
}

