##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

gx.barplot <- function(x, main="", cex.main=1.2, cex.names=0.85,
                       cex.legend=0.9, srt=0, xlab='', ylab='',
                       group=NULL, group.names=NULL, bar.names=NULL,
                       voff=2, legend=TRUE)
{
    ##col1 = rev(grey.colors(2))
    if(0) {
        main="";cex.main=1.2;cex.names=0.85
        cex.legend=0.9;srt=0;xlab='';ylab='';
        group=NULL;legend=TRUE
    }
    cpal = rep(RColorBrewer::brewer.pal(12,"Paired"),99)
    if(is.null(group)) group <- rep(1,ncol(x))
    group2 <- as.vector(sapply(group,rep,2))
    col1 <- cpal[(group2-1)*2 + rep(1:2,NCOL(x))]
    if(NCOL(x)==1) col1 = col1[2]
    ylim <- c(0, 1.15*max(x,na.rm=TRUE))
    str.ht <- strheight("XXX",cex=cex.names)
    str.ht <- 0.04*ylim[2]
    space <- c(0,0.4)
    is.multi = (NCOL(x)>1)
    is.multi
    if(!is.multi) space <- c(0.2,0)
    
    barplot(x, beside=TRUE, col=col1,
            main=main, cex.main=cex.main,
            xlab=xlab, ylab=ylab, ylim=ylim,
            ##names.arg = rep('',ncol(x))
            las=1, space=space,
            names.arg = rep('',length(x)),
            ## names.arg = rep(rownames(x),ncol(x)),
            mgp=c(2,0.9,0) )
    
    if(is.multi) {
        sp0 <- nrow(x)+space[2]
        tx0 <- sp0*(0:(ncol(x)-1)) + 1
        tx1 <- tx0 + 0.5*nrow(x) + space[1]
        tx  <- as.vector(mapply(tx0,tx1,FUN=c))
        tx2 <- (tx0+tx1)/2
        if(is.null(bar.names)) {
            bar.names <- rep(rownames(x),NCOL(x))
            bar.names <- toupper(substring(bar.names,1,3))
        }
        mtext(bar.names, side=1, line=0, at=tx, adj=0.5, cex=0.5, srt=0, xpd=TRUE)
        if(is.null(group.names)) {
            group.names <- colnames(x)
        }
        if(srt==0) {
            mtext(group.names, side=1, line=1.5, at=tx2,
                  adj=0.5, cex=cex.names, srt=0, xpd=TRUE)
        } else {
            text( tx2, -voff*str.ht, group.names, cex=cex.names,
                 srt=srt, xpd=TRUE, pos=2, offset=0)
        }        

    } else {

        text(1.2*(1:length(x))-0.5, -voff*str.ht, names(x), cex=cex.names,
             srt=srt, xpd=TRUE, pos=2, offset=0*str.ht)

    }
    
    if(legend) {
        nn <- names(x)
        if(is.multi) nn <- rownames(x)
        legend("topleft", legend=nn, fill=col1,
               xpd=TRUE, bty='n', inset=c(-0.0,-0.04), cex=cex.legend,
               y.intersp=0.8, x.intersp=0.3)
    }
    
}


gx.b3plot <- function(x, y,
                      first = NULL,
                      width = 1,
                      bar = TRUE,
                      bee = TRUE,
                      sig_stars = FALSE,
                      ymax = NULL,
                      bee_cex = 2,
                      max_stars = 5,
                      srt = NULL,
                      xoff = 0,
                      names_cex = 1,
                      names = TRUE,
                      max_points = 100,
                      col = "#a6cee3F0",
                      mult_x_lim = 0.11,
                      mult_y_lim = 1.3,
                      ylab = 'CPM',
                      xlab = '',
                      pch = 19,
                      ...)
{
    # browser()

    # use the boxplot function to extract whiskers --------------------
    stats_segments <- function(x, y,
                               xoffset = 0,
                               lwd = 2) {
        bx <- boxplot(y ~ x, plot = FALSE)
        nx <- length(bx$n)
        x0 <- xoffset + (1:nx)

        segments(x0 - 0.1, bx$conf[1, ], x0 + 0.1, bx$conf[1, ], lwd = lwd)
        segments(x0 - 0.1, bx$conf[2, ], x0 + 0.1, bx$conf[2, ], lwd = lwd)
        segments(x0, bx$conf[1, ], x0, bx$conf[2, ], lwd = lwd * 0.5)
    }

    # Define NAs in factor levels ---------------------------------------
    ylevel <- levels(y)
    y <- as.character(y)

    if (any(is.na(y))) {
        y[is.na(y)] <- 'NA'
        ylevel <- c(ylevel,'NA')
    }

    y <- factor(y, levels = ylevel, exclude = NULL)
    if (!is.null(first)) y <- relevel(y, ref = first)

    # Condition medians --------------------------------------------------
    mx <- tapply(x, y, median, na.rm = TRUE)

    # ??? ----------------------------------------------------------------
    sig <- yc <- NULL


    # Not sure about this part -------------------------------------------
    # does it add t-test stars? Is it ever used? -------------------------
    # if (sig_stars) {
    #     y_levels <- unique(y)
    #     yc <- combn(y_levels, 2)
    #     pv <- rep(NA, ncol(yc))
    #     i <- s1
    #     for (i in 1:ncol(yc)) {
    #         grp <- yc[, i]
    #         pv[i] <- t.test(x[which(y == grp[1])],
    #                         x[which(y == grp[2])])$p.value
    #     }
    #     pv
    #     sig <- c("", "*", "**", "***")[1 + 1 * (pv < 0.05) + 1 * (pv < 0.01) + 1 * (pv < 0.001)]
    #     sig
    #     nthree <- sum(sig == "***")
    #     jj <- which(sig != "")
    #     jj <- jj[order(pv[jj])] ## only top 4 ??
    #     jj <- Matrix::head(jj, max(max.stars, nthree)) ## only top 5 ??
    #     j <- 1
    #     yc <- apply(yc, 2, as.integer)
    #     yc
    #     dd <- abs(as.vector(diff(yc)))
    #     jj <- jj[order(dd[jj])]
    #     sig <- sig[jj]
    #     yc <- yc[, jj, drop = FALSE]
    # }

    # compute graphical paramenters ------------------------------------
    dx <- max(x, na.rm = TRUE) * mult_x_lim
    ylim <- c(xoff, max(x) * mult_y_lim)
    if (!is.null(ymax)) ylim <- c(xoff, ymax)
    if (min(x) < 0) {
        ylim  <- c(mult_y_lim * min(c(x, xoff)),
                   max(x) * mult_y_lim)
        }

    if(length(col) == 1) {
        col <- rep(col, length(mx))
        names(col) <- names(mx)
    } else {
        col <- col[match(names(mx), names(col))]
        col[is.na(col)] <- "grey90"
    }

    # plot beeswarm -----------------------------------------------------
    jj <- 1:length(x)
    if (max_points > 0 &&
          length(jj) > max_points) {
            jj <- unlist(
                    tapply(jj, y, function(i) {
                        Matrix::head(sample(i),
                                     max_points)
                      })
                    )
        }

    # svg()
    x_levels <- factor(ylevel,
                       levels = ylevel)
    boxplot(mx ~ x_levels,
            ylim = ylim,
            xlab = xlab,
            ylab = ylab,
            notch = TRUE, # make median lines smaller
            frame.plot = FALSE,
            offset = xoff)
    beeswarm::beeswarm(x[jj] ~ y[jj],
                       ylim = ylim,
                       pch = pch,
                       cex = bee_cex,
                       corral = "random",
                       col = col,
                       add = TRUE)
    # dev.off()


    # if(is.null(srt)) {
    #     nnchar <- sum(sapply(unique(y),nchar))
    #     srt <- ifelse(nnchar > 24, 30, 0)
    # }
    # 
    # pos <- ifelse(srt==0, 1, 2)
    # 
    # n = length(unique(y))
    # if(names==TRUE) {
    #     y0 = min(ylim) - diff(ylim)*0.05
    #     text( bx[,1], y0, names(mx), cex=names_cex,
    #          srt=srt, adj=ifelse(srt==0,0.5,0.965), xpd=TRUE,
    #          pos=pos, offset=0)
    # }
    # 
    # if(bar) 
    # 
    # if(sig_stars) {
    #     i=1
    #     for(i in 1:NCOL(yc)) {
    #         grp = yc[,i]
    #         xmax = max(x,na.rm=TRUE)*1.05 + dx*i
    #         j1 = grp[1] - 0.4
    #         j2 = grp[2] - 0.4
    #         segments( j1, xmax, j2, xmax, lwd=0.5)
    #         if(ncol(yc)<=8)
    #             text((j1+j2)/2, xmax, labels=sig[i], pos=1, offset=-0.33, adj=0, cex=1.4)
    #     }
    # }
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
