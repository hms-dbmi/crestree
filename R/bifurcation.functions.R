##' Assign cells in a probabilistic manner to non-intersecting windows along pseudotime
##' @param  r ppt tree
##' @param root root of progenitor branch of bifurcation
##' @param leaves leaves of derivative branches of bifurcation
##' @param wind number of cell per local pseudotime window
##' @param mapping boolean, if project cells onto tree pseudotime in a probabilistic manner. TRUE by default.
##' @param regions manually assign coordinates of windows
##' @param n.cores number of cores to use
##' @return a list of cell probabilisties of being in each local pseudotime window
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


##' Assign cells in a probabilistic manner to non-intersecting windows along pseudotime
##' @param freq list of per-window cell probabilities. Outcome of slide.cells
##' @param mat expression matrix
##' @param genesetA a vector of fate-specific genes for first post-bifurcation branch
##' @param genesetB a vector of fate-specific genes for second post-bifurcation branch
##' @return a list of local per-window matrices of average gene correlations with two fate-specific gene modules
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



##' Visualize cells assigned to each local window
##' @param emb embedding
##' @param freq cell probabilities of assignment to local windows.
##' @return ggplot2 object
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


##' Visualize patterns of local correlations in consequtive windows.
##' @param  cors a list of local average correlations with fate-specific modules. Outcome of slide.cors function.
##' @param genesetA a vector of fate-specific genes of bifurcation branch 1
##' @param genesetB a vector of fate-specific genes of bifurcation branch 2
##' @param val boolean, if show the degree of "repulsion" of fate-specific modules. FALSE by default.
##' @return ggplot object visualising pattern of local correlations in consequtive windows.
##' @export
fig.cors <- function(cors,genesetA,genesetB,val=FALSE){
  cr <- unlist(lapply(cors,function(corc) as.vector(corc)))
  cr.abs <- max(abs(cr))
  lapply( cors,function(corc){
    corA <- corc[,1]
    corB <- corc[,2]
    colr = rep("red",length(corA)); colr[ c(genesetA,genesetB) %in% genesetA ] = "blue"
    repulsion = cor(corA,(ifelse( c(genesetA,genesetB)%in%genesetA,1,2 )))+cor(corB,(ifelse( c(genesetA,genesetB)%in%genesetB,1,2 )))
    return(ggplot()+geom_point( aes(corA,corB),colour=colr,size=1)+
             #xlim(-0.3,0.3)+ylim(-0.3,0.3)+ # this is for sens-auto
             xlim(-cr.abs,cr.abs)+ylim(-cr.abs,cr.abs)+ # this is for auto-mes
             geom_vline(xintercept = 0,linetype=2)+geom_hline(yintercept = 0,linetype=2)+
             #annotate("text",x=0.2,y=0.25,label=paste(round(repulsion/2,1),sep=""),size=5,fontface=3)+ # this is for sens-auto
             annotate("text",x=cr.abs-0.1,y=cr.abs-0.05,label=paste(round(repulsion/2,1),sep=""),size=ifelse(val==TRUE,5,0),fontface=3)+ # this is for auto-mes
             theme(legend.position="none")+theme_bw()+theme(axis.title=element_blank(),axis.ticks=element_blank(),axis.text=element_blank())+
             theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    )
  })
}



##' Estimates pseudotime trends of local intra- and inter-module correlations of fates-specific modules
##' @param ppt ppt object
##' @param fpm expression matrix
##' @param root root of progenitor branch of bifurcation
##' @param leaves leaves of derivative branches of bifurcation
##' @param genesetA a vector of fate-specific genes for first post-bifurcation branch
##' @param genesetB a vector of fate-specific genes for second post-bifurcation branch
##' @param w local window, in number of cells, to estimate correlations
##' @param step step, in number of cells, between local windows
##' @param n.mapping number of probabilistic cells projection to use for estimates
##' @param n.points number of pseudotime points to extrapolate and visualize patterns
##' @param span.smooth smooth parameters of loess
##' @param perm boolean, estimate control trends for local permutations instread of real expression matrix
##' @return estimates of local correlation trend for each probabilistic cell projection and average among cell projections
##' @export
synchro <- function(ppt,fpm,root,leaves,genesetA,genesetB,w,step,n.mapping,n.points=100,span.smooth = 0.1,perm=FALSE){
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


##' Estimates pseudotime trends of local intra- and inter-module correlations of fates-specific modules along a single trajectory
##' @param ppt ppt object
##' @param mat expression matrix
##' @param genesetA a vector of fate-specific genes for first post-bifurcation branch
##' @param genesetB a vector of fate-specific genes for second post-bifurcation branch
##' @param w local window, in number of cells, to estimate correlations
##' @param step step, in number of cells, between local windows
##' @param subtree trajectory, consisting of a vector of segments comprising trajectory
##' @param perm boolean, estimate control trends for local permutations instread of real expression matrix
##' @return estimate of local correlation trends along trajectory
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


##' Visualise local correlation trends of inter- and intra- module local correlations of fate-specific modules
##' @param  outcome of synchro function
##' @return visualize combiation of ggplot2 objects
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



##' Estimates inclusion times of genes in a correlated module for a number of probabilistic cell projections
##' @param  ppt ppt tree
##' @param geneset a set of genes that form a module
##' @param nodes tips of a tree that form subtree (e.g, trajectory or fork)
##' @param expr expression matrix
##' @param alp parameter regulating stringency of inclusion event
##' @param w local window of cells along pseudotime to estimate local correlation
##' @param step shift of a window along pseudotime in number of cells
##' @param n.cells total number of cells to consider from root in progenitor branch
##' @param mappings number of probabilistic cell projections to assess
##' @param do.perm do local estimates for locally permuted expression matrix
##' @param winperm number of local cells for permutations
##' @param n.cores number of cores
##' @param permut.n number of permutations used to estiamte background local correlations in onset.est
##' @param n.cores1 number of cores to calculate permutations used to estimate background local correlations in onset.est
##' @return matrix of inclusion timing for each gene (rows) in each probabilistic cells projection (columns)
##' @export
onset =  function(ppt,geneset,nodes=NULL,expr,alp=20,w=40,step=10,n.cells=280,mappings=1,do.perm=FALSE,winperm=w,n.cores=1,permut.n=10,n.cores1=1){
  if (is.null(nodes)){nodes <- ppt$tips}
  res <- do.call(cbind,mclapply(mappings, function(perm){
    cat("mapping: ");cat(perm);cat("\n")
    res <- onset.est(ppt,perm,geneset,nodes,expr,alp=alp,w=w,step=step,do.perm=do.perm,winperm=winperm,n.cells=n.cells,
                     permut.n=permut.n,n.cores=n.cores1)
    return(res)
  },mc.cores = n.cores))
}


##' Estimates inclusion times of genes in a correlated module
##' @param ppt ppt tree
##' @param mapping id of probabilistic cells project to use
##' @param geneset a set of genes that form a module
##' @param nodes tips of a tree that form subtree (e.g, trajectory or fork)
##' @param expr expression matrix
##' @param alp parameter regulating stringency of inclusion event
##' @param w local window of cells along pseudotime to estimate local correlation
##' @param step shift of a window along pseudotime in number of cells
##' @param do.perm do local estimates for locally permuted expression matrix
##' @param winperm number of local cells for permutations
##' @param n.cells total number of cells to consider from root in progenitor branch
##' @param permut.n of local cells for permutations
##' @param n.cores number of cores to calculate permutations used to estimate background local correlations in onset.est
##' @param autocor.cut cutoff on correlation of inclusion times between sequential iterations of the algorithm to stop it
##' @param iterations maximum number of iterations of the algorithm
##' @return vector of inclusion times of genes
##' @export
onset.est <- function(ppt,mapping=1,geneset,nodes,expr,alp=20,w=40,step=10,do.perm=FALSE,winperm=30,n.cells=280,
                      permut.n=10,n.cores=1,autocor.cut=0.95,iterations=15){
  # segments of subtree
  segs <- extract.subtree(ppt,nodes)$segs

  # select first n.cells by pseudotime
  subtbl <- ppt$cell.info[[mapping]][ppt$cell.info[[mapping]]$seg %in% segs,]
  cells <- rownames(subtbl)[order(subtbl$t)]
  cells <- cells[1:min(n.cells,length(cells))]

  # initialize inclusion matrix
  logW = matrix(1,nrow=length(geneset),ncol=length(cells)/step+1); rownames(logW) = geneset

  # locally permute expression matrix if do.perm=TRUE
  expr.use <- expr
  if (do.perm==TRUE){
    cell_order = cells
    for ( i in seq(1,length(cell_order)-winperm,winperm) ){
      expr.use[,cell_order[i:(i+winperm-1)]] <- t(apply(expr[,cell_order[i:(i+winperm-1)]],1,function(x){
        sample(x)
      }))}
  }

  # estimate local Spearman gene-gene correlations (real and locally permuted control)
  corMatr <- corMatrix(geneset,cells,expr.use[,cells],step,w)
  corMatrC <- corMatrixC(geneset,cells,expr.use[,cells],step,w,permutN = permut.n,n.cores=n.cores)

  # iterative estimation of inclusion time.
  logList <- list()
  autocor <- 0
  i <- 1
  auto <- vector()
  while ( (is.na(autocor) | autocor < autocor.cut) & i <= iterations){
    #print(i)
    pMat = corEvol( corMatr,corMatrC,logW)
    logList[[i]]=switchPoint(pMat,alp)
    rownames(logList[[i]]) = geneset
    logW=logList[[i]]
    if (i>1) {
      autocor=cor(apply(logList[[i-1]],1,which.max),apply(logList[[i]],1,which.max),method="spearman")
      auto[i]=autocor
    }
    #print(paste("autocor",autocor,sep=" "))
    i=i+1
  }
  #print(paste("autocor",autocor,"interations",i,sep=" "))

  # estimate time of inclusion
  apply(logW,1,function(x) if (sum(x)==0) {NA}else {
    pos = which.max(x)
    mean(subtbl[cells,]$t[(1+pos*(step-1)):(pos*(step-1)+w)])
  })
}


##' Estimates local gene-gene correlations in a sequence of windows along pseudotime
##' @param genes a set of genes to estimate local pairwise expression correlations
##' @param cells cells ordered by pseudotime
##' @param cdnormX expression matrix
##' @param step step between sequential windows (in number of cells)
##' @param w window length along pseudotime (in number of cells)
##' @return  local gene-gene correlation matrices per each window of cells
##' @export
corMatrix = function(genes,cells,cdnormX,step,w){
  corMat = lapply( seq(1,(length(cells)-w),step),function(wind){
    cls0 = cells[wind:min(wind+w-1,length(cells))]
    cdnorm1=cdnormX[genes,cls0];
    #cdnorm1[ cdnorm1==0 ] = NA;
    cdnorm2=cdnorm1;cdnorm2[ is.na(cdnorm2)]=0;cdnorm2[ cdnorm2!=0]=1;cdnorm2=cdnorm2%*%t(cdnorm2);
    cdnorm2[cdnorm2<10]=0;cdnorm2[cdnorm2>=10]=1;
    cormat=cor( t(cdnorm1), use="pairwise.complete.obs",method="spearman" ) * cdnorm2; cormat[is.na(cormat)]=0
    #cormat = cor( t(cdnormX[genes,cls0]),use="pairwise.complete.obs",method="spearman" ); cormat[is.na(cormat)]=0
    rownames(cormat)=genes
    return(cormat)
  })
  return(corMat)
}

##' Estimates local gene-gene correlations in a sequence of windows along pseudotime for locally permuted matrix
##' @param genes a set of genes to estimate local pairwise expression correlations
##' @param cells cells ordered by pseudotime
##' @param cdnormX expression matrix
##' @param step step between sequential windows (in number of cells)
##' @param w window length along pseudotime (in number of cells)
##' @param permutN number of local permutations
##' @param n.cores number of cores to use
##' @return local multiple permuted gene-gene correlation matrices per each window of cells
##' @export
corMatrixC = function(genes,cells,cdnormX,step,w,permutN,n.cores){
  corMat = mclapply( seq(1,(length(cells)-w),step),function(wind){
    cls0 = cells[wind:min(wind+w-1,length(cells))]
    lapply(1:permutN, function(i){
      cdnorm1=t(apply( cdnormX[genes,cls0],1,function(x) sample(x) ));colnames(cdnorm1)=cls0
      #cdnorm1[ cdnorm1==0 ] = NA;
      cdnorm2=cdnorm1;cdnorm2[ is.na(cdnorm2)]=0;cdnorm2[ cdnorm2!=0]=1;cdnorm2=cdnorm2%*%t(cdnorm2);
      cdnorm2[cdnorm2<10]=0;cdnorm2[cdnorm2>=10]=1;
      cormat=cor( t(cdnorm1), use="pairwise.complete.obs",method="spearman" ) * cdnorm2; cormat[is.na(cormat)]=0
      #cormat = cor( t(cdnorm1[genes,cls0]),use="pairwise.complete.obs",method="spearman" ); cormat[is.na(cormat)]=0
      rownames(cormat)=genes
      return(cormat)})
  },mc.cores = n.cores)
  return(corMat)
}


##' Estimates switch point in correlation trend
##' @param matr matrix time-ordered scores (columns) of genes (rows)
##' @param alp parameter regulating stringency of switch event
##' @return vector of gene-specific switch points
##' @export
switchPoint = function(matr,alp){ t(apply(matr,1,function(x){
  logL = unlist(lapply(1:(length(x)+1),function(j){
    v = 0
    if ( j >= 2) v = v + 0#-sum( 0.5 - 2*x[1:(j-1)],na.rm=TRUE ) #+ 0 #+ 0.25*sum(is.na(x[1:(j-1)]))#+ sum( log(x[1:(j-1)]+0.05),na.rm=TRUE )  #-sum( 0.5 - 2*x[1:(j-1)],na.rm=TRUE ) #+ 0.25*sum(is.na(x[1:(j-1)]))
    if (j <= length(x)) v = v + sum( -alp*x[j:length(x)]+log( alp/(1-exp(-alp)) ) ,na.rm=TRUE) #- sum( 2*x[j:length(x)]-0.5,na.rm=TRUE )#+ sum( -alp*x[j:length(x)]+log( alp/(1-exp(-alp)) ) ,na.rm=TRUE) #- 0.25*sum(is.na(x[j:length(x)]))#+ (length(x)-j)*log(alp) #- sum( 2*x[j:length(x)]-0.5,na.rm=TRUE ) #- 0.25*sum(is.na(x[j:length(x)]))
    return(v)
  }))
  logL[1:length(x)][is.na(x)]=NA
  switch_point = which.max(logL)
  return( c(rep(0.00,switch_point-1),rep(1, length(x)+2-switch_point)) )
}))}



##' Estimate switch point in correlation trend
##' @param corMatr list of local gene-gene correlation matrices per window
##' @param corMatrC list of multiple control local gene-gene correlation matrices per window
##' @param logW matrix of binary step-wise indicators of gene inclusion in correlated module
##' @return estimates of empirical p-values of genes inclusion in correlated module
##' @export
corEvol = function( corMatr,corMatrC,logW){
  cor_Pmat = do.call(cbind,lapply( 1:length(corMatr),function(ipos){

    cormat = corMatr[[ipos]]
    corTrue = apply( t(t(cormat)*logW[,ipos]),1,mean);

    cor_control = do.call(cbind,lapply(1:length(corMatrC[[ipos]]),function(i){
      cormat = corMatrC[[ipos]][[i]]
      corB = apply( t(t(cormat)*logW[,ipos]),1,mean);
      return(corB)
    }))

    cor_P = unlist(lapply( 1:length(corTrue),function(i){
      if (corTrue[i]==0 | sum(logW[,ipos])<2 ) {NA}else {
        sum( (cor_control[i,]+runif(length(cor_control[i,]),0,0.01)) >= corTrue[i] )/length(cor_control[i,])}
    }))
    return(cor_P)
  } ))
  #rownames(cor_Pmat)=genes
  return(cor_Pmat)
}

##' Show cumulative probability of gene inclusion in correlated module
##' @param gn gene name to show
##' @param matSwitch matrix of gene inclusion time in different probabilistic cells projections
##' @param bins number of bins
##' @param par boolean. Use external plot parameters par() if TRUE. FALSE by default.
##' @export
show.gene.inclusion <- function(gn,matSwitch,bins=8,par=FALSE){
  if (!gn%in%rownames(matSwitch)){stop('gene is not in matrix')}
  sg=seq(0,max(matSwitch),(max(matSwitch))/bins)
  logHist = t(apply(matSwitch,1,function(x){unlist(lapply(sg,function(y) sum( x< y )/length(x) ))}))

  if (par==FALSE) {par(mar=c(4,7,3,1),mgp=c(2.5,0.5,0))}
  #gn = "Neurog2"
  lbl <- (2:ncol(logHist))/ncol(logHist)*max(matSwitch); #lbl <- round(lbl,2)
  plot( lbl,logHist[gn,2:length(sg)],lwd=6,
        type="l",yaxt="na",xlab="pseudotime",ylab="inclusion\nprobability",cex.lab=1.0,main=gn,cex.main=1.0,font.main=3)
  axis(side=2,at=c(0,0.5,1),cex.axis=1.5)
  abline( h=seq(0,1,0.5),lty = 2,lwd=2,col="grey60")
}

##' Show summary cumulative probability of genes inclusion in correlated module
##' @param matSwitch matrix of gene inclusion time in different probabilistic cells projections
##' @param gns.show genes to show. By default all genes are shown
##' @param bins number of bins
##' @param par boolean. Use external plot parameters par() if TRUE. FALSE by default.
##' @export
show.inclusion.summary <- function(matSwitch,gns.show=rownames(matSwitch),bins=8,par=FALSE){
  sg=seq(0,max(matSwitch),(max(matSwitch))/bins)
  logHist = t(apply(matSwitch,1,function(x){unlist(lapply(sg,function(y) sum( x< y )/length(x) ))}))
  swOrder = apply(matSwitch,1,mean)

  if (par==FALSE) {par(mar=c(4,0.5,1,7),mgp=c(2,0.5,0))}
  image( t(logHist[order(swOrder),2:(length(sg)-1)]),zlim=c(0,1),col=colorRampPalette(c("grey", "darkgreen"))(n = 60),
         labels=FALSE,lab.breaks=rownames(logHist[order(swOrder),]),yaxt="n",xaxt="n",xlab="pseudotime",cex.lab=1 )
  nms = rownames(logHist[order(swOrder),]); nms[!(nms %in% gns.show)] = ""
  seqq = seq(0,1,length.out=nrow(logHist))[! nms==""]
  axis( 4, at=seqq, labels=nms[! nms==""],hadj=0.0,xaxt="s",cex.axis=1,font = 3,las= 1)
  axis( 1, at=seq(0,1,length.out = 4),labels=round(seq(min(matSwitch),max(matSwitch),length.out = 4),2),hadj=0.3)
  box()
}


##' Show summary cumulative probability of genes inclusion in correlated module
##' @param alps levels of alp parameter
##' @param incl.real average inclusion time for alp levels
##' @param incl.perm average inclusion time for alp levels in permuted matrix of expression levels
##' @export
inclusion.stat <- function(alps, incl.real, incl.perm){
  par(mar=c(4,4,1,1))
  mm <- max(c(incl.real,incl.perm))
  plot(alps,mm-incl.real,pch=19,cex=1,col="red",type="l",ylim=c(0,max(c(mm-incl.real,mm-incl.perm))),xlab="alp",ylab="total early inclusion",cex.lab=1);
  points(alps,mm-incl.real,pch=19,cex=1,col="red")
  lines(alps,mm-incl.perm,pch=19,cex=1,col="blue")
  points(alps,mm-incl.perm,pch=19,cex=1,col="blue")
  legend("topright",c("real","permutations"),fill=c("red","blue"),bty="n")
}

