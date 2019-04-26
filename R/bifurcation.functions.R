##' Estimate probabil
##' @param  n.samples a number of seed samplings.
##' @param seeds a vector of seeds to use. Overwrites n.samples.
##' @return a list of pptree objects
##' @export
slide.cells <- function(r,root,leaves,wind=50,mapping=TRUE,regions=NULL,n.cores=10){
  segs <- extract.subtree(r,c(root,leaves))
  if (is.null(regions)){
    segs.prog <- intersect(extract.subtree(r,c(root,leaves[1]))$segs,extract.subtree(r,c(root,leaves[2]))$segs)
    segs.b1 <- setdiff(extract.subtree(r,c(root,leaves[1]))$segs,segs.prog)
    segs.b2 <- setdiff( segs$segs,c(segs.prog,segs.b1))
    pp.probs <- colSums(r$R)
    pps <- r$pp.info$PP[r$pp.info$seg %in% segs$segs]

    # recursive function to estimate cell probabilisties in sliding window along the tree
    region.extract <- function(pt.cur,segs.cur){
      freq <- list()
      pps.next <- pps[r$pp.info[pps,]$time >= pt.cur & r$pp.info[pps,]$seg %in% segs.cur]
      cmsm <- cumsum( (pp.probs[pps.next])[order(r$pp.info[pps.next,]$time)] )
      inds <- which(cmsm > wind)
      #is.na(inds[1])
      if ( is.na(inds[1])  ){
        if ( max(cmsm) > wind/2 ){
          if (mapping==TRUE){
            cell.probs <- apply(r$R[,pps.next],1,sum)
          }else{
            cell.probs <- as.numeric(apply(r$R[,],1, function(x) which.max(x) %in% pps.next  ))
            names(cell.probs) <- rownames(r$R)
          }
          freq <- c(freq,list(cell.probs))
          #pt.cur <- NULL; segs.cur <- NULL
          return(freq)
        }
      }else{
        pps.region <- pps.next[order(r$pp.info[pps.next,]$time)][1:inds[1]]
        if ( mapping==TRUE ){
          cell.probs <- apply(r$R[,pps.region],1,sum)
        }else{
          cell.probs <- as.numeric(apply(r$R[,],1, function(x) which.max(x) %in% pps.region  ))
          names(cell.probs) <- rownames(r$R)
        }
        freq <- c(freq,list(cell.probs))
        pt.cur <- max(r$pp.info[pps.region,]$time)
        table(r$pp.info[pps.region,]$seg)

        if ( sum(!r$pp.info[pps.region,]$seg %in% segs.prog)==0 ){
          res <- region.extract(pt.cur,segs.cur)
          return( c(freq,(res)) )
        }else if ( sum(!r$pp.info[pps.region,]$seg %in% segs.b1)==0 ){
          segs.cur <- segs.b1
          res <- region.extract(pt.cur,segs.cur)
          return( c(freq,(res)) )
        }else if ( sum(!r$pp.info[pps.region,]$seg %in% segs.b2)==0 ){
          segs.cur <- segs.b2
          res <- region.extract(pt.cur,segs.cur)
          return( c(freq,(res)) )
        }else if ( !sum(!r$pp.info[pps.region,]$seg %in% segs.prog)==0 ){
          pt.cur1 <- max(r$pp.info[pps.region,]$time[r$pp.info[pps.region,]$seg%in%segs.b1])
          segs.cur1 <- segs.b1
          pt.cur2 <- max(r$pp.info[pps.region,]$time[r$pp.info[pps.region,]$seg%in%segs.b2])
          segs.cur2 <- segs.b2
          res1 <- region.extract(pt.cur1,segs.cur1)
          res2 <- region.extract(pt.cur2,segs.cur2)
          return( c(freq,(res1),(res2)) )
        }
      }
    }
    pt.cur <- min(r$pp.info[pps,]$time)
    segs.cur <- segs$segs
    freq <- region.extract(pt.cur,segs.cur)
    #cell.probs <- res[[1]]
    #show_points(names(cell.probs)[cell.probs > 0.2])
    return(freq)
  }else{
    freq = mclapply( 1:length(regions),function(j){
      x = do.call(rbind,lapply( 1:length(r$cell.info),function(i){
        cell.pseudotime <- r$cell.info[[i]]
        x = cell.pseudotime[ cell.pseudotime$seg %in% unlist(regions[[j]][[1]]), ];
        cls = rownames(r$R)[x[order(-regions[[j]][[4]]*x$t),]$cell[ regions[[j]][[2]]:regions[[j]][[3]] ]]
        return( as.numeric(rownames(r$R)%in%cls) )
      }))
      freq = apply(x,2,mean); names(freq) = rownames(r$R)
      return(freq)
    },mc.cores = n.cores)
  }
}


##' Estimate probabil
##' @param  n.samples a number of seed samplings.
##' @param seeds a vector of seeds to use. Overwrites n.samples.
##' @return a list of pptree objects
##' @export
slide.cors <- function(freq,mat,genesetA,genesetB){
  lapply( 1:length(freq),function(j){
    fpm1 = t(apply(mat[c(genesetA,genesetB),],1,function(x) rank(x)))
    cormat = cov.wt( t(fpm1[c(genesetA,genesetB),names(freq[[j]])]),wt=freq[[j]],cor=T ); cormat$cor[is.na(cormat$cor)]=0
    corA = apply(cormat$cor[,genesetA],1,mean); corB = apply(cormat$cor[,genesetB],1,mean)
    corA[genesetA] = (corA[genesetA] - 1/length(genesetA))*length(genesetA)/(length(genesetA)-1)
    corB[genesetB] = (corB[genesetB] - 1/length(genesetB))*length(genesetB)/(length(genesetB)-1)
    return(cbind(corA,corB))
  })
}



##' Estimate probabil
##' @param  n.samples a number of seed samplings.
##' @param seeds a vector of seeds to use. Overwrites n.samples.
##' @return a list of pptree objects
##' @export
fig.cells <- function(emb,freq){
  lapply( 1:length(freq),function(j){
    fre <- rep(0,nrow(emb)); names(fre) <- rownames(emb)
    fre[names(freq[[j]])] <- freq[[j]]
    colr = colorRampPalette(c('lightgrey','black'))(100)[as.numeric(cut( fre,breaks=100))]
    return(ggplot()+geom_point( aes(emb[,1],emb[,2]),color=colr,size=0.6 )+theme_bw()+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+
             theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    )
  })
}


##' Estimate probabil
##' @param  n.samples a number of seed samplings.
##' @param seeds a vector of seeds to use. Overwrites n.samples.
##' @return a list of pptree objects
##' @export
fig.cors <- function(cors,genesetA,genesetB,val=FALSE){
  lapply( cors,function(corc){
    corA <- corc[,1]
    corB <- corc[,2]
    colr = rep("red",length(corA)); colr[ c(genesetA,genesetB) %in% genesetA ] = "blue"
    repulsion = cor(corA,(ifelse( c(genesetA,genesetB)%in%genesetA,1,2 )))+cor(corB,(ifelse( c(genesetA,genesetB)%in%genesetB,1,2 )))
    return(ggplot()+geom_point( aes(corA,corB),colour=colr,size=1)+
             #xlim(-0.3,0.3)+ylim(-0.3,0.3)+ # this is for sens-auto
             xlim(-0.25,0.25)+ylim(-0.25,0.25)+ # this is for auto-mes
             geom_vline(xintercept = 0,linetype=2)+geom_hline(yintercept = 0,linetype=2)+
             #annotate("text",x=0.2,y=0.25,label=paste(round(repulsion/2,1),sep=""),size=5,fontface=3)+ # this is for sens-auto
             annotate("text",x=0.15,y=0.2,label=paste(round(repulsion/2,1),sep=""),size=ifelse(val==TRUE,5,0),fontface=3)+ # this is for auto-mes
             theme(legend.position="none")+theme_bw()+theme(axis.title=element_blank(),axis.ticks=element_blank(),axis.text=element_blank())+
             theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    )
  })
}



##' Estimate probabil
##' @param  n.samples a number of seed samplings.
##' @param seeds a vector of seeds to use. Overwrites n.samples.
##' @return a list of pptree objects
##' @export
fig.cor.evol <- function(cors,freq,emb,genesetA,genesetB){
  fig_cells <- fig.cells(emb,freq)
  fig_cor <- fig.cors(cors,genesetA,genesetB)
  pl <- marrangeGrob( t(c(fig_cells,fig_cor)),ncol=length(fig_cells),nrow=2,
                      layout_matrix = matrix(seq_len(2*length(fig_cells)), nrow = 2, ncol = length(fig_cells),byrow=TRUE),top=NA)
  return(pl)
}


##' Estimate probabil
##' @param  n.samples a number of seed samplings.
##' @param seeds a vector of seeds to use. Overwrites n.samples.
##' @return a list of pptree objects
##' @export
synchro <- function(ppt,fpm,root,leaves,genesetA,genesetB,w,step,n.mapping,n.points=100,span.smooth = 0.75,perm=FALSE){
  cat('process path 1 of 2');cat('\n')
  subtree <- extract.subtree(ppt,c(root,leaves[1]))
  xx1 <- synchro.path(ppt,fpm,genesetA,genesetB,w,step,n.mapping = n.mapping,subtree,perm=perm);
  subtree <- extract.subtree(ppt,c(root,leaves[2]))
  cat('process path 2 of 2');cat('\n')
  xx2 <- synchro.path(ppt,fpm,genesetA,genesetB,w,step,n.mapping = n.mapping,subtree,perm=perm);
  xx2$run <- xx2$run + max(xx1$run)

  pts <- seq(min(c(xx1$t,xx2$t)),max(c(xx1$t,xx2$t)),length.out = n.points)
  xx.fit <- lapply( list(xx1,xx2),function(xx.temp){
    ftAA <- loess( xx.temp$corAA ~ xx.temp$t,span=span.smooth ); corAA <- predict(ftAA,pts)
    ftAB <- loess( xx.temp$corAB ~ xx.temp$t,span=span.smooth ); corAB <- predict(ftAB,pts)
    ftBB <- loess( xx.temp$corBB ~ xx.temp$t,span=span.smooth ); corBB <- predict(ftBB,pts)
    segs <- unlist(lapply(pts,function(pt) xx.temp$seg[which.min(abs(xx.temp$t-pt))] ))
    cols <- unlist(lapply(pts,function(pt) xx.temp$colr[which.min(abs(xx.temp$t-pt))] ))
    xx.fit <- data.frame(pts,corAA = corAA,corAB = corAB,corBB = corBB,seg = segs,col = cols)
    return(xx.fit)
  })
  xx.fit <- rbind(xx.fit[[1]],xx.fit[[2]])

  xx.smooth <- do.call(rbind,apply( unique(cbind(xx.fit$pts,xx.fit$seg)),1,function(rw){
    rv <- xx.fit[xx.fit$pts==rw[1] & xx.fit$seg==rw[2],]
    data.frame(rw[1],mean(rv[,2]),mean(rv[,3]),mean(rv[,4]),rw[2],rv[1,6])
  }))
  colnames(xx.smooth) <- colnames(xx.fit)

  xx <- rbind(xx1,xx2)

  return( list(xx.smooth,xx) )

}


##' Estimate probabil
##' @param  n.samples a number of seed samplings.
##' @param seeds a vector of seeds to use. Overwrites n.samples.
##' @return a list of pptree objects
##' @export
synchro.path = function(r,mat,genesetA,genesetB,w,step,n.mapping = length(r$cell.info),subtree,perm=FALSE){
  xx1=do.call(rbind,mclapply( 1:n.mapping,function(j){
    cell.pseudotime <- r$cell.info[[j]]
    texpr1 = mat;
    if (perm==TRUE){
      for (seg in subtree$segs){
        tloc = cell.pseudotime[cell.pseudotime$seg==seg,]
        tloc = tloc[order(tloc$t),]
        cell_order = rownames(r$R)[tloc$cell]
        winperm=min(5,length(cell_order))
        for ( i in seq(1,length(cell_order)-winperm+1,winperm) ){
          texpr1[,cell_order[i:(i+winperm-1)]] = t(apply(mat[,cell_order[i:(i+winperm-1)]],1,function(x){
            sample(x)
          }))
        }
      }
    }

    cinfo = cell.pseudotime[cell.pseudotime$seg %in% subtree$segs,]
    xx = do.call(rbind,lapply( seq(1,(nrow(cinfo)-w),step),function(i){
      #cls = rownames(r$R)[ cinfo$cell[order(cinfo$t)][i:(i+w)] ]
      cls = rownames(cinfo)[ order(cinfo$t)[i:(i+w)] ]

      cormat = cor( t(texpr1[c(genesetA,genesetB),cls]),method="spearman",use="na.or.complete" ); cormat[is.na(cormat)]=0
      corA = apply(cormat[,genesetA],1,mean); corB = apply(cormat[,genesetB],1,mean)
      corA[genesetA] = (corA[genesetA] - 1/length(genesetA))*length(genesetA)/(length(genesetA)-1)
      corB[genesetB] = (corB[genesetB] - 1/length(genesetB))*length(genesetB)/(length(genesetB)-1)
      t = mean(cinfo$t[order(cinfo$t)][i:(i+w)])
      seg = cinfo$seg[which.min(abs(cinfo$t-t))]; colr = cinfo$color[which.min(abs(cinfo$t-t))]
      c( t, (mean(corA[genesetA])-mean(corA[genesetB]))^2+(mean(corB[genesetA])-mean(corB[genesetB]))^2,
         mean(corA[genesetA]),mean(corB[genesetB]),mean(corA[genesetB]),j,seg,colr)
    }))
  },mc.cores = 1))
  xx1 = data.frame(xx1);
  colnames(xx1) = c("time","dist","corAA","corBB","corAB","run","seg","colr")
  xx1$time=as.numeric(as.character(xx1$time));
  xx1$dist=as.numeric(as.character(xx1$dist));
  xx1$corAA=as.numeric(as.character(xx1$corAA));
  xx1$corAB=as.numeric(as.character(xx1$corAB));
  xx1$corBB=as.numeric(as.character(xx1$corBB));
  xx1$run=as.numeric(as.character(xx1$run));
  xx1$seg=as.numeric(as.character(xx1$seg))
  return(xx1)
}


##' Estimate probabil
##' @param  n.samples a number of seed samplings.
##' @param seeds a vector of seeds to use. Overwrites n.samples.
##' @return a list of pptree objects
##' @export
visualize.synchro <- function(crd){
  consensus <- crd[[1]]
  indiv <- crd[[2]]

  pl_aa <- ggplot()+geom_point(aes(consensus$pts,consensus$corAA),colour=consensus$col)+
    geom_line(aes( indiv$time,indiv$corAA,group=indiv$run),colour=indiv$colr,size=0.3,alpha=I(0.1))+
    theme_bw()+theme(legend.position="none")+xlab("")+ylab("module A")+
    theme(axis.title.y=element_text(size=16),axis.text=element_text(size=20),axis.title.x=element_blank())+
    scale_y_continuous(limits=c(0,max(indiv$corAA,consensus$corAA)))+geom_hline(yintercept = 0,linetype=2)

  pl_bb <- ggplot()+geom_point(aes(consensus$pts,consensus$corBB),colour=consensus$col)+
    geom_line(aes( indiv$time,indiv$corBB,group=indiv$run),colour=indiv$colr,size=0.3,alpha=I(0.1))+
    theme_bw()+theme(legend.position="none")+xlab("")+ylab("module B")+
    theme(axis.title.y=element_text(size=16),axis.text=element_text(size=20),axis.title.x=element_blank())+
    scale_y_continuous(limits=c(0,max(indiv$corBB,consensus$corBB)))+geom_hline(yintercept = 0,linetype=2)

  pl_ab <- ggplot()+geom_point(aes(consensus$pts,consensus$corAB),colour=consensus$col)+
    geom_line(aes( indiv$time,indiv$corAB,group=indiv$run),colour=indiv$colr,size=0.3,alpha=I(0.1))+
    theme_bw()+theme(legend.position="none")+xlab("")+ylab("modules A-B")+
    theme(axis.title.y=element_text(size=16),axis.text=element_text(size=20),axis.title.x=element_blank())+
    scale_y_continuous(limits=c(min(indiv$corAB,consensus$corAB),max(indiv$corAB,consensus$corAB)))+geom_hline(yintercept = 0,linetype=2)

  grid.arrange( pl_aa,pl_bb,pl_ab,bottom=textGrob("pseudotime",gp=gpar(cex=2)),
                left=textGrob("local correlation",rot=90,gp=gpar(cex=2)),ncol=1, nrow = 3 )
}

