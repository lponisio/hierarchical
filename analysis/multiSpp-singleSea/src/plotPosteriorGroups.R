
plotPosteriorGroups <- function(){
    layout(matrix(1:3, ncol=1))
    par(oma=c(3,1,1,1))
    cols <- rainbow(dim(ms.ss.occ.all$ms.ss$summary)[1])
    for(group in groups){
        if(length(group) > 1){
            for(i in 1:dim(ms.ss.occ.all$ms.ss$summary)[1]){
                this.samp <- ms.ss.occ.all$ms.ss$summary[i,,group]
                xs <- jitter(1:dim(this.samp)[2])
                if(i == 1){
                    plot(x=xs, y=this.samp["mean",], pch=16,
                         col=cols[i],
                         ylim=range(c(ms.ss.occ.all$ms.ss$summary[,'CI95_upp',
                                                                  group],
                                      ms.ss.occ.all$ms.ss$summary[,'CI95_low',
                                                                  group])),
                         xaxt="n",
                         ylab="Estimate",
                         xlab="")

                    axis(1, at=1:dim(this.samp)[2],
                         labels=FALSE)
                    text(x=1:dim(this.samp)[2],
                         y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
                         labels=group, srt=45, adj=1, xpd=TRUE)
                } else{
                    points(x=xs,
                           y=this.samp["mean",], pch=16, col=cols[i])
                }
                arrows(y1=this.samp['CI95_upp',],
                       y0=this.samp['CI95_low',],
                       x0=xs,
                       code=0, angle=90, length=0.02, lwd=1, col=cols[i])
            }
        }

        legend("topright", legend=dimnames(ms.ss.occ.all$ms.ss$summary)[[1]],
               pch=16, col=cols)
    }
}
