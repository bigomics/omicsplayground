

##x=gx;y=ngs$samples$group;bee=bar=TRUE;offx=3;sig.stars=FALSE;xoff=0;srt=60
gx.b3plot <- function(x, y, width=1, bar=TRUE, bee=TRUE, sig.stars=FALSE,
                      ymax=NULL, bee.cex=0.3, max.stars=5, srt=NULL, xoff=0,
                      names.cex=1, names=TRUE, max.points=1000, col="grey80", ...)
{
    require(beeswarm)
    ##require(sinaplot)
    stats.segments <- function(x, y, xoffset=0, lwd=2) {
        bx = boxplot(y ~ x, plot=FALSE)
        nx = length(bx$n)
        x0 = xoffset + (1:nx)
        ##segments( x0-0.2, bx$stats[3,], x0+0.2, bx$stats[3,],lwd=lwd*1.5)
        segments( x0-0.1, bx$conf[1,],  x0+0.1, bx$conf[1,],lwd=lwd)
        segments( x0-0.1, bx$conf[2,],  x0+0.1, bx$conf[2,],lwd=lwd)
        segments( x0,     bx$conf[1,],  x0,     bx$conf[2,],lwd=lwd*0.5)
    }
    y = factor(y)
    mx = tapply(x, y, median)

    sig = yc = NULL
    if(sig.stars) {
        y.levels = unique(y)
        yc = combn(y.levels,2)
        pv <- rep(NA,ncol(yc))
        i=1
        for(i in 1:ncol(yc)) {
            grp = yc[,i]
            pv[i] = t.test(x[which(y==grp[1])], x[which(y==grp[2])])$p.value
        }
        pv
        sig = c("","*","**","***")[ 1 + 1*(pv<0.05) + 1*(pv<0.01) + 1*(pv<0.001) ]
        sig
        nthree = sum(sig=="***")
        jj = which(sig!="")
        jj = jj[order(pv[jj])]  ## only top 4 ??
        jj = head(jj,max(max.stars,nthree)) ## only top 5 ??
        j=1
        yc = apply(yc, 2, as.integer)
        yc
        dd = abs(as.vector(diff(yc)))
        jj = jj[order(dd[jj])]
        sig = sig[jj]
        yc = yc[,jj,drop=FALSE]
    }

    ##dx = (max(x,na.rm=TRUE)-min(x,na.rm=TRUE))*0.11
    dx = max(x,na.rm=TRUE)*0.11
    ylim = c(xoff,max(x)*1.3)
    if(!is.null(ymax)) ylim = c(xoff,ymax)
    if(min(x) < 0) ylim = c(1.3*min(c(x,xoff)),max(x)*1.3)
    if(sig.stars) {
        if(ncol(yc) > 8) dx = dx/5
        ylim[2] = ylim[2]*1.05 + (2+NCOL(yc))*dx
    }

    ##par(mfrow=c(1,1));srt=60
    ##bx = barplot( mx-xoff, width=0.6666, space=0.5, ylim=ylim, offset=xoff, names.arg=NA)
    bx = barplot( mx, width=0.6666, space=0.5, ylim=ylim, offset=xoff,  names.arg=NA, col=col, ... )
    y0 = min(ylim) - diff(ylim)*0.08
    pos <- ifelse(srt==0, 1, 2)

    n = length(unique(y))
    if(is.null(srt)) srt <- ifelse(n>3, 30, 0)
    if(names==TRUE) {
        text( bx[,1], y0, names(mx), cex=names.cex, adj=1, srt=srt, xpd=TRUE,
             pos=pos, offset=0)
    }
    if(bee) {
        jj <- 1:length(x)
        ##jj <- which(y %in% names(which(table(y)>2)))
        j1 <- which(table(jj)==1)
        if(length(j1)) jj <- c(jj,j1)
        if(max.points>0 && length(jj)>max.points) jj <- sample(jj,max.points)
        beeswarm(x[jj] ~ y[jj], add=TRUE, at=1:n-0.33, pch=19, cex=bee.cex, col="grey20")
        ## sinaplot( x[jj] ~ y[jj], add=TRUE, pch=19, cex=bee.cex, col="grey20")
    }
    if(bar) stats.segments(y, x, xoffset=-0.333, lwd=1.2)

    if(sig.stars) {
        i=1
        for(i in 1:NCOL(yc)) {
            grp = yc[,i]
            xmax = max(x,na.rm=TRUE)*1.05 + dx*i
            j1 = grp[1] - 0.4
            j2 = grp[2] - 0.4
            segments( j1, xmax, j2, xmax, lwd=0.5)
            if(ncol(yc)<=8)
                text((j1+j2)/2, xmax, labels=sig[i], pos=1, offset=-0.33, adj=0, cex=1.4)
        }
    }
}

gx.hist <- function(gx, main="",ylim=NULL) {
    h0 <- hist(as.vector(gx), breaks=120, main=main,
               col="grey",freq=FALSE, ylim=ylim, xlab="signal")
    i = 1
    for(i in 1:ncol(gx)) {
        h1 <- hist(gx[,i], breaks=h0$breaks,plot=FALSE)
        lines( h0$mids, h1$density, col=i+1 )
    }
}
