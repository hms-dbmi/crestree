#' @useDynLib crestree
NULL

##' Sample pptree objects using different seeds
##' @param  n.samples a number of seed samplings.
##' @param seeds a vector of seeds to use. Overwrites n.samples.
##' @return a list of pptree objects
##' @export
mppt.tree <- function( ... , n.cores=parallel::detectCores()/2,n.samples=n.cores, seed=NULL,seeds=NULL) {
  if(!is.null(seed)) {
    set.seed(seed);
  }
  # sample seeds
  if(is.null(seeds)) {
    seeds <- round(runif(n.samples,0,.Machine$integer.max))
  }
  mclapply(seeds,function(i) ppt.tree(..., seed=i),mc.cores=n.cores)
}

##' Sample pptree objects using bootstrap
##' @param X expression matrix of genes (rows) and cells (columns).
##' @param M number of principal points of pptree.
##' @param n.samples a number of seed samplings.
##' @param replace sampling with replacement (logical).
##' @return a list of pptree objects
##' @export
bootstrap.ppt <- function( ..., X, M=ncol(X),n.cores=parallel::detectCores()/2,n.samples=n.cores, seed=NULL,replace=T) {
  if(!is.null(seed)) {
    set.seed(seed);
  }
  parallel::mclapply(1:n.samples,function(i) {
    # take a bootstrap sample
    b.X <- X[,sample(1:ncol(X),M,replace=replace)];
    ppt.tree(..., X=b.X, M=M, init=b.X)
  },mc.cores=n.cores)
}

##' Calculate weighted pairwise correlations between columns of matrices A and B
##' @export
cor.mat <- function(A,B){
  A1 <- t(t(A)-colMeans(A))
  B1 <- t(t(B)-colMeans(B))
  res <- (crossprod(A1,B1))/sqrt( tcrossprod(colSums(A1^2),(colSums(B1^2))) )
  return(res)
}

##' Calculate pairwise euclidean distances between columns of matrices A and B
euclidean.mat <- function(A,B){
  x <- do.call(cbind,rep(list(colSums(A^2)),ncol(B)))
  y <- do.call(rbind,rep(list(colSums(B^2)),ncol(A)))
  suppressWarnings(res <- sqrt(x + y - 2*crossprod(A,B)))
  res[is.na(res) | is.nan(res)] <- 0
  return(res)
}

##' calculate weighted correlation between columns of a matrix and a given vector
wcr <- function(X,y,w){
  w <- w/sum(w)
  X1 <- X*w
  y1 <- y*w
  X2 <- t(t(X)-colSums(X1))
  y2 <- y - sum(y1)

  cv1 <- (y2*w)%*%X2
  cv2 <- sqrt(colSums(X2^2*w)*sum(y2^2*w))
  cvv <- cv1/cv2
  return(cvv[1,])
}

##' Reconstruction of the tree
##'
##' Using SimplePPT approach to model principal tree (pptree) of the data
##' @name ppt.tree
##' @param X gene (row) vs cell (columns) expression matrix
##' @param emb embdedding to visalize cells and principal tree together
##' @param M number of principal points to use (more than zero, no more than number of cells)
##' @param init matrix of initial gene coordinates of principal points
##' @param plot plot or not intermediate trees
##' @param lambda penalty for the tree  length, as used in SimplePPT
##' @param sigma parameter as used in SimplePPT
##' @param seed used to make initial assignment of principal points to a subset of cells
##' @param n.steps number of iteraions
##' @param metrics metrics used to calculated distances between cells or principal points. "euclidean" or "cosine"
##' @param p.power if cosine metrics used, option p.power allows to use (1-cor)^p.power (p.power=1 by default)
##' @param err.cut stop algorithm if proximity of principal points between iterations less than err.cut
##' @return pptree object
##' @export
ppt.tree <- function(X,W=NA,emb=NA,M,init=NULL,plot=TRUE,output=TRUE,lambda=1e1,sigma=0.1,seed=NULL,n.steps=50,err.cut = 5e-2,metrics="cosine",p.power=1,knn=NULL,...) {

  if ( metrics!="euclidean" & metrics!="cosine" ){ stop("metrics paramterer is nethier 'euclidean' nor 'cosine'") }
  if ( M < 0 | M > ncol(X)) { stop("M should be more than zero and less or equal than the number of cells") }
  if (!is.na(emb)){
    if ( sum(!colnames(X)%in%rownames(emb))>0 ) { stop("column names of gene expression matrix (X) are not consistent with row names of embedding (emb)") }
  }

  X <- as.matrix(X)
  wt <- TRUE
  if (is.na(W)) {
    wt <- FALSE
    W <- matrix(1,nrow=nrow(X),ncol=ncol(X))
  }else{
    W <- as.matrix(W[rownames(X),colnames(X)])
  }

  if(is.null(init)){
    if(!is.null(seed)){
      set.seed(seed);
    }
    F.mat <- X[,sample(1:ncol(X),M)]; rownames(F.mat) <- NULL; colnames(F.mat) <- NULL;
  } else {
    F.mat <- init;
  }

  # row-normalize W
  rwm <- matrix(rowSums(W),nrow=nrow(F.mat),ncol(F.mat))
  W <- W/rowSums(W)*ncol(W);


  # repeat untile convergence
  j=1; err=100;
  while(j <= n.steps & err > err.cut) {
    # calculate R
    if (metrics=="euclidean"){
      # simple correlation or column-wise weighted correlation.
      if (wt==FALSE) {
        R <- euclidean.mat(F.mat,X)^p.power
      }else{
        R <- do.call(cbind,lapply(1:ncol(X),function(i) {
          sqrt(colSums(((F.mat-X[,i])^2)*W[,i]))^p.power
        }))
      }
      R <- t(exp(-R/sigma))
    }else if(metrics=="cosine"){
      # simple correlation or column-wise weighted correlation.
      if (wt==FALSE) {
        cordist <- (1-cor.mat(F.mat,X))^p.power
      }else{
        cordist <- do.call(cbind,lapply(1:ncol(X),function(i) {
          (1-matWVCorr(F.mat,X[,i],W[,i]))^p.power
          #(1-wcr(F.mat,X[,i],W[,i]))^p.power
        }))
        colnames(cordist) <- colnames(X)
      }

      cordist <- (cordist-mean(cordist))
      R <- t(exp( -(cordist)/sigma ))
    }
    R[is.na(R) | is.nan(R)] <- 0
    if (!is.null(knn)){
      R = apply(R,2,function(x){
        x[ x < sort(x,decreasing = TRUE)[knn] ] <- 0
        x
      })
    }
    R <- R/rowSums(R)
    R[is.na(R) | is.nan(R)] <- 0

    # calculate distance between principal points
    if (metrics=="euclidean"){
      d <- euclidean.mat(F.mat,F.mat)
    }else if (metrics=="cosine"){
      if (wt==FALSE) {
        d <-  1-cor.mat(F.mat,F.mat)
      }
      else{
        d <- do.call(cbind,lapply(1:ncol(F.mat),function(i) {
          (1-matWVCorr(F.mat,F.mat[,i],rwm[,i]))^p.power
          #(1-wcr(F.mat,F.mat[,i],rwm[,i]))^p.power
        }))
      }
      d <- abs(d)^p.power*sign(d)
    }
    bt <- minimum.spanning.tree(graph.adjacency(as.matrix(d),weighted=T,mode="undirected"))

    B <- as.matrix(get.adjacency(bt))
    D <- diag(nrow(B))*rowSums(B)
    L <- D-B
    M <- L*lambda + diag(ncol(R))*colSums(R)

    old.F <- F.mat;
    #F.mat <- (X%*%R) %*% chol2inv(chol(M))
    F.mat <- t(solve( t(M),t((X*W)%*%R) ))# slightly faster, 15%
    F.mat <- as.matrix(F.mat)

    if (plot==TRUE){plotppt(list(F=F.mat,B=B,R=R,L=L,lambda=lambda,sigma=sigma),emb,...)}

    if (output==TRUE){
      cat(j,":")
      cat("\n")
      err = max(sqrt(colSums(F.mat-old.F)^2)/apply(F.mat,2,function(x)sqrt(sum(x^2))))
      cat(err,"\n")
    }
    j=j+1
  }

  if (plot==TRUE){plotppt(list(F=F.mat,B=B,R=R,L=L,lambda=lambda,sigma=sigma),emb,...)}

  g = graph.adjacency(B,mode="undirected");tips = V(g)[degree(g)==1];forks = V(g)[degree(g)>2]

  score = c( sum( t(1-cor.mat(F.mat,X))*R)/nrow(R), sigma/nrow(R)*sum(R*log(R),na.rm=T),lambda/2*sum(d*B))

  colnames(R) <- colnames(F.mat) <- rownames(B) <- colnames(B) <- as.character(1:nrow(B))
  invisible(list(score=score,F=F.mat,B=B,R=R,L=L,DT=d,lambda=lambda,sigma=sigma,n.steps=n.steps,metrics=metrics,M=M,cells=vi,tips=tips,forks=forks))
}


##' Estimate optimal sigma parameter.
##'
##' Using cross-validation criteria to select sigma parameter.
##' @param X gene (rows) vs cell (columns) expression matrix
##' @param M number of principal points in pptree modeling
##' @param n.sample number of sampled trees per each sigma
##' @param sig.lims a vector of sigma for which cross-validation estimated
##' @param metrics similarity measure. "cosine" or "euclidean"
##' @return optimal sigma parameter
##' @export
sig.explore <- function(X,W=NA,M=as.integer(ncol(X)/2),n.sample=1,sig.lims=seq(0.01,0.2,0.03),metrics="cosine",p.power = 1,plot=TRUE,err.cut=5e-1,n.steps=20,n.cores=1){
  if (is.na(X)) {stop("matrix X should be specified")}
  if (is.na(M)) {stop("number of principal points M should be specified")}
  cells <- colnames(X)
  for (i in 1:n.sample){
    cv <- do.call(rbind,mclapply(sig.lims,function(sig){
      x <- ppt.tree(X = X,W,M=M,err.cut=err.cut,metrics=metrics,n.steps=n.steps,p.power = p.power,lambda=0,sigma=sig,plot=FALSE,output=FALSE,seed=sample(100,1))
      y <- cor(X,x$F)
      apply(y,1,max)
    },mc.cores = n.cores))
    if (i==1){
      cv.tot <- cv
    }
    else{
      cv.tot <- cv.tot + cv
    }
  }
  cv.tot <- cv.tot/n.sample

  sig.opt <- sig.lims[which.max(apply(cv.tot,1,mean))]
  if (plot==TRUE){
    par(mfrow=c(1,1),mar=c(5,5,1,1))
    plot( sig.lims, apply(cv.tot,1,mean),lty=2,lwd=2,type="l",xlab="sigma",ylab="CV",cex.lab=1.5)
    points( sig.lims, apply(cv.tot,1,mean),pch=19,cex=1)
    abline(v=sig.opt,col="red",lty=2)
  }
  #return( cbind(sig.lims,apply(cv.tot,1,mean)) )
  return(sig.opt)

}

##' Explore lambda
##'
##' Explores multiple lambda and choose the optimal
##' @param X gene (rows) vs cell (columns) expression matrix
##' @param M number of principal points in pptree modeling
##' @param sigma fixed parameter sigma used in pptree modeling
##' @param emb embdedding to visalize cells and principal tree together. If emb is given than pptrees for a range of lambda are shown
##' @export
lambda.explore <- function(X=NA,M=ncol(X),sigma=0.1,emb=NA,metrics="cosine",tips.min=2,tips.max=10,base=2,lambda.init=100,err.cut=5e-3,n.steps=40,p.power=1){
  if (is.na(X)) {stop("matrix X should be specified")}
  if (is.na(M)) {stop("number of principal points M should be specified")}
  cells <- colnames(X)

  min.reached <- FALSE;max.reached <- FALSE
  lambda <- round(lambda.init)
  tr.list <- list()
  while (min.reached==FALSE | max.reached==FALSE){
    print(paste("lambda:",round(lambda,2) ))
    tr <- ppt.tree(X=X,M=M,lambda=lambda,sigma=sig,err.cut=err.cut,metrics=metrics,n.steps=n.steps,p.power = p.power,plot=FALSE,output=FALSE,seed=sample(100,1))
    tr <- setroot(tr,root=as.character(tr$tips[1]))
    tr.list[[as.character(round(lambda,1))]] <- tr#c(tr.list,tr)

    tips <- length(tr$tips);
    len <- sum(tr$pp.segments$d)
    entropy.ind <- sum(tr$pp.segments$d*log(tr$pp.segments$d))
    # add entry to the lambda.info matrix
    if (lambda == lambda.init){
      lambda.info <- matrix(c(lambda=lambda,tips=tips,length=len,entropy=entropy.ind),nrow=1,ncol=4)
      #tr.list[[as.character(lambda)]] <- tr
    }else{
      if (lambda < lambda.info[1,1]){
        lambda.info <- rbind(c(lambda=lambda,tips=tips,length=len,entropy=entropy.ind),lambda.info)
        #tr.list[[as.character(lambda)]] <- tr#c(tr,tr.list)
      }else{
        lambda.info <- rbind(lambda.info,c(lambda=lambda,tips=tips,length=len,entropy=entropy.ind))
        #tr.list[[as.character(lambda)]] <- #c(tr.list,tr)
      }
    }
    # update lambda
    if (min.reached == FALSE & tips < tips.max){
      lambda <- lambda/base
    }else if (min.reached == FALSE & tips >= tips.max){
      min.reached <- TRUE
      lambda <- lambda.info[nrow(lambda.info),1]*base
    }else if (tips <= tips.min ){# | tips >= lambda.info[nrow(lambda.info)-1,2]){
      max.reached <- TRUE
    }else{
      lambda <- lambda.info[nrow(lambda.info),1]*base
    }
  }

  ent.per.tip <- lambda.info[,4]/lambda.info[,2]
  i.opt <- which.min(ent.per.tip)
  if (!is.na(emb)){
    par(mfrow=c(2,2))
    par(mar=c(5,5,1,1))
    plot( lambda.info[,1], ent.per.tip,log="x",lty=2,lwd=2,type="l",xlab="lambda",ylab="entropy per tip",cex.lab=1.5)
    points(lambda.info[,1], ent.per.tip,pch=19,cex=1)
    abline(v=lambda.info[i.opt,1],col="red",lty=2)
    par(mar=rep(1,4))
    lamb <- lambda.info[i.opt,1]; lamb <- round(lamb,1)
    plotppt(tr.list[[as.character(lamb)]],emb,cex.tree = 0.1,lwd.tree = 3,main=paste("lambda =",lamb))
    box(col="red",lwd=3);
    lamb <- lambda.info[median(1:i.opt),1]; lamb <- round(lamb,1)
    plotppt(tr.list[[as.character(lamb)]],emb,cex.tree = 0.1,lwd.tree = 3,main=paste("lambda =",lamb))
    lamb <- lambda.info[median((i.opt+1):nrow(lambda.info)),1]; lamb <- round(lamb,1)
    plotppt(tr.list[[as.character(lamb)]],emb,cex.tree = 0.1,lwd.tree = 3,main=paste("lambda =",lamb))
  }
  return(lambda.info)
  #return(list(lambda.info[i.opt,1],lambda.info))
}


##' Visualize pptree onto embedding
##'
##' Projects pptree onto embedding (e.g. tSNE)
##' @name plotppt
##' @param r - pptree object
##' @param emb - (x,y) coordinates data frame (e.g Rtsne $Y result)
##' @param F - coordinates of principal points (optional)
##' @param gene - a gene to show expression of (optional)
##' @param mat - gene vs cell expression matrix (needed if option 'gene' is activated)
##' @param pattern.cell - numeric profile of a quantity for each cell (e.g. expression of a gene or cell cycle stage)
##' @param pattern.tree - numeric profile of a quantity for each principal point (e.g. expression of a gene or cell cycle stage)
##' @param cex.main - cex of points
##' @param cex.col - color of points
##' @param cex.title - cex of title
##' @param cex.tree - cex of principal points
##' @param tips - logical, to draw indecies of tips of the tree. Usefull before usage of cleanup.branches()
##' @export
plotppt <- function(r,emb,F=NULL, gene=NULL, main=gene, mat=NULL, pattern.cell=NULL, pattern.tree=NULL,
                    cex.col=NA, tree.col = NULL,
                    cex.main=0.5, cex.title=1,
                    cex.tree=1.5,lwd.tree=1,par=TRUE,tips=FALSE,forks=FALSE,subtree=NA,...) {
  if ( sum(!rownames(r$R)%in%rownames(emb))>0 ) { stop("cell names used for tree reconstruction are not consistent with row names of embedding (emb)") }
  if (sum(!is.na(cex.col))==0 ) {cex.col=rep("grey70",nrow(emb)); names(cex.col) <- rownames(emb)}
  vi = rownames(emb)%in%rownames(r$R); names(vi) <- rownames(emb)
  if(is.null(F)) { F <- t(t(t(emb[rownames(r$R),])%*%r$R)/colSums(r$R)) }
  if ( is.null(pattern.cell) & !is.null(gene) ){
    if (is.null(mat)) { stop("mat expression matrix should be defined together with gene parameter") }
    if (gene %in% rownames(mat) == FALSE) { stop("gene is not in mat matrix") }
    if ( sum(!rownames(r$R) %in% colnames(mat)) > 0 ) { stop("cell names used for tree reconstruction are not consistent with mat column names") }
    pattern.cell = mat[gene,rownames(r$R)]#mat[gene,rownames(r$R)]
  }

  if ( !is.null(pattern.tree) & length(pattern.tree) != ncol(r$R) ) { stop("length of pattern.tree vector is inconsistent with cell number used for tree reconstruction") }
  if ( !is.null(pattern.cell) & is.null(pattern.tree) ){
    if ( sum(!names(pattern.cell) %in% rownames(r$R)) > 0 ){ stop("pattern.cell vector should contain names for all cells used to reconstruct the tree")}
    pattern.cell <- pattern.cell[rownames(r$R)] ## is it correct?
    aggr <- colSums(r$R)
    pattern.tree <- t(r$R)%*%pattern.cell[rownames(r$R)]/aggr
    pattern.tree[aggr==0] <- NA
  }

  if (is.null(tree.col)) {tree.col = "black"}
  if( !is.null(pattern.cell) ){
    cex.col <- rep("black",nrow(emb)); names(cex.col) <- rownames(emb)
    cex.col[names(pattern.cell)] <- colorRampPalette(c("blue","gray50","red"))(1024)[round((pattern.cell-min(pattern.cell))/diff(range(pattern.cell))*1023)+1]
    #cex.col <- colorRampPalette(c("blue","gray50","red"))(1024)[round((pattern.cell-min(pattern.cell))/diff(range(pattern.cell))*1023)+1]
  }
  if ( !is.null(pattern.tree) ){
    tree.col <- colorRampPalette(c("blue","gray50","red"))(1024)[round((pattern.tree-min(pattern.tree,na.rm=T))/diff(range(pattern.tree,na.rm = T))*1023)+1]
    #r$fitting$pp.fitted[gene,]
  }

  if (!is.na(subtree)){
    #cex.col[rownames(r$cell.summary)][!r$cell.summary$seg %in% subtree$seg] <- "black"
    tree.col[!r$pp.info$seg %in% subtree$seg] <- "grey80"
    vi[vi==TRUE][rownames(r$cell.summary)][!r$cell.summary$seg %in% subtree$seg] <- FALSE
  }

  if ( sum(names(cex.col)%in%rownames(emb))==0 ) {stop('cex.col names do not match row names of emb')}

  cols <- rep("black",nrow(emb)); names(cols) <- rownames(emb)
  cols[ intersect(names(cex.col),rownames(emb)) ] <- cex.col[intersect(names(cex.col),rownames(emb))]
  if (par==TRUE) {par(mar=rep(1,4))}
  plot(emb,pch=ifelse(vi,19,1),cex=cex.main,col = adjustcolor(cols,ifelse(is.null(pattern.tree),1,0.1)),xlab=NA,ylab=NA,xaxt='n',yaxt='n',main=main,cex.main=cex.title,font.main=1)

  al <- get.edgelist(graph.adjacency(r$B>0))
  al <- matrix(as.integer(al),ncol=2)
  segments(F[1,al[,1]],F[2,al[,1]],F[1,al[,2]],F[2,al[,2]],lwd=lwd.tree)
  points(t(F),pch=21,
         col=tree.col,bg=tree.col,cex=cex.tree)

  if (tips==TRUE){
    coord = do.call(rbind,lapply(r$tips,function(tip){
      x1 = F[1,tip]; y1 = F[2,tip]
      x2 = F[1,which(r$B[tip,]>0)]; y2 = F[2,which(r$B[tip,]>0)]
      xnew = x1 + 1.5*sign(x1-x2)#(1+sign(x1-x2)/0.5)*sign(x1-x2)#alpha*(x1-x2)
      ynew = y1 + 1.5*sign(y1-y2)#xnew*(y2-y1)/(x2-x1) + (y1*x2-y2*x1)/(x2-x1)
      c(xnew,ynew)
    }))
    text((coord),col=1,cex=1,adj=c(0,0),labels=r$tips,font=2);#text(t(F[, r$tips ]),col=1,cex=1.2,adj=c(0,0),labels=r$tips);
  }
  if (forks==TRUE & length(r$forks) > 0){
    coord = do.call(rbind,lapply(r$forks,function(fork){
      x1 = F[1,fork]; y1 = F[2,fork]
      x2 = F[1,which(r$B[fork,]>0)]; y2 = F[2,which(r$B[fork,]>0)]
      xnew = x1 #+ 1.5*sign(x1-x2)#(1+sign(x1-x2)/0.5)*sign(x1-x2)#alpha*(x1-x2)
      ynew = y1 #+ 1.5*sign(y1-y2)#xnew*(y2-y1)/(x2-x1) + (y1*x2-y2*x1)/(x2-x1)
      c(xnew,ynew)
    }))
    text((coord),col=1,cex=1,adj=c(0,0),labels=r$forks,font=2);#text(t(F[, r$tips ]),col=1,cex=1.2,adj=c(0,0),labels=r$tips);
  }
  #legend(x="bottomright",legend=c(paste("lambda=",r$lambda[1],sep=""),paste("sigma=",r$sigma[1],sep="")))
}


##' Visualize list of pptree objects onto embedding
##'
##' Projects pptree objects onto embedding (e.g. tSNE)
##' @param rl list of pptree objects (as calculated using bootstrap.tree or mppt.tree)
##' @param emb (x,y) coordinates data frame (e.g Rtsne $Y result)
##' @param cols vector of colors for cells in emb.
##' @export
plotpptl <- function(rl,emb, cols=adjustcolor(1,alpha=0.3),alpha=1, lwd =1, ...) {
  par(mfrow=c(1,1), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 0.8);
  plot(emb,col=cols,cex=1,pch=19,xlab="",ylab="", ...)
  lapply(rl,function(r) {
    F <- t(t(t(emb[rownames(r$R),])%*%r$R)/colSums(r$R))
    al <- get.edgelist(graph.adjacency(r$B>0))
    al <- matrix(as.integer(al),ncol=2)
    #points( t(F),col=adjustcolor(cols,alpha=0.1),lwd=1,cex=0.2 )
    segments(F[1,al[,1]],F[2,al[,1]],F[1,al[,2]],F[2,al[,2]],lwd=lwd,col=adjustcolor("black",alpha))
  })
  #legend(x="bottomright",legend=c(paste("lambda=",rl[[1]]$lambda[1],sep=""),paste("sigma=",rl[[1]]$sigma[1],sep="")))
}


##' Remove spurious branches of pptree
##' @param r ppt.tree result
##' @param tips.number select and retain only fixed number of tips (tips.number) that explain the most cell-cell variation.
##' @param tips.remove vector of tips indices to remove
##' @param min.branch.length remove all branches with length less or equal than min.branch.length principal points
##' @return modified ppt.tree object with cleaned up structure
##' @export
cleanup.branches <- function(r,tips.remove=NULL,min.branch.length=3) {
  #colnames(r$F) <- NULL; colnames(r$B) <- rownames(r$B) <- NULL;
  repeat {
    g <- graph.adjacency(r$B>0,mode="undirected")
    leaves <- V(g)[degree(g)==1]
    branches <- V(g)[degree(g)>2]
    bd <-shortest.paths(g,v=leaves,to=branches)
    ivi <- which(apply(bd,1,min)<=min.branch.length)
    ivi <- unique( c(ivi, which( leaves %in% tips.remove) ) )
    if(length(ivi)==0) { break }
    toremove <- c();
    for(x in ivi) {
      bdp <- get.shortest.paths(g,leaves[x],to=branches[which.min(bd[x,])])
      toremove <- c(toremove,bdp$vpath[[1]][-length(bdp$vpath[[1]])])
    }

    # remove from the graph (B)
    r$B <- r$B[-toremove,-toremove]
    # remove from F
    r$F <- r$F[,-toremove];
    # remove from lRu
    r$lRu <- r$lRu[,-toremove]
    # remove from R and renormalize
    r$R <- r$R[,-toremove];
    r$R <- r$R/rowSums(r$R);
  }
  colnames(r$F) <- colnames(r$B) <- rownames(r$B) <- as.character(1:nrow(r$B));

  g = graph.adjacency(r$B,mode="undirected");r$tips = V(g)[degree(g)==1];r$forks = V(g)[degree(g)>2]

  r
}


##' Orient the tree by setting up the root
##'
##' Assign root, pseudotime and segment to each principal point of the tree
##' @param r pptree object
##' @param root root principal point (plotppt(tips=TRUE,..) can be used to visualize candidate tips for a root)
##' @return modified ppt.tree object with new fields r$pp.info (estimated pseudotime and branch of principal points), r$pp.segments (segments information), r$root (root id).
##' @export
setroot <- function(r,root=NULL,plot=TRUE) {
  if (is.null(root)) { stop("Assign correct root number") }
  if ( ! root %in% r$tips ) {stop("Root should be one of the tree tips")}
  # calculate time of each PP
  if (r$metrics=="euclidean"){d <- 1e-6+euclidean.mat(r$F,r$F)
  }else if (r$metrics=="cosine"){
    d <-  abs( 1e-2 + 1-cor.mat(r$F,r$F))
  }
  g <- graph.adjacency(r$B*d,weighted=T,mode="undirected")
  pp.info <- data.frame( cbind( V(g),as.double(shortest.paths(g,root,V(g))),rep(0,length(V(g))) ));
  colnames(pp.info)=c("PP","time","seg")

  # infer all segments (and put in segs) of the tree
  nodes <- V(g)[ degree(g)!=2 ]
  pp.segs = data.frame(n=numeric(),from=character(),to=character(),d=numeric())
  for (i in 1:(length(nodes)-1) ){
    for (j in (i+1):length(nodes)){
      node1 = nodes[i];node2=nodes[j];
      path12 = unlist(get.shortest.paths(g,from=as.character(node1),to=as.character(node2)))
      if ( sum(nodes %in% path12) == 2  ) {
        from = node1$name;to=node2$name
        if ( !is.null(root)){
          path_root = shortest.paths(g,root,c(node1,node2))
          from = colnames(path_root)[which.min(path_root)]
          to = colnames(path_root)[which.max(path_root)]
        }
        pp.info[path12,]$seg = nrow(pp.segs)+1
        pp.segs=rbind(pp.segs,data.frame(n=nrow(pp.segs)+1,from=from,to=to,d=shortest.paths(g,as.character(node1),as.character(node2))[1]))
      }}}
  pp.segs$color=rainbow(nrow(pp.segs))
  pp.info$color=pp.segs$color[pp.info$seg]

  r$pp.segments <- pp.segs;
  r$root <- root;
  r$pp.info <- pp.info
  r
}



##' Project cells onto the principal tree
##' @param r pptree object
##' @param emb if not NULL than cell branch assignment and color code of branches are shown
##' @param n.mapping number of probabilistic mapping of cells onto the tree to use. If n.mapping=1 then likelihood cell mapping is used.
##' @return modified pptree object with new fields r$cell.summary, r$cell.info and r$img.list. r$cell.summary contains information about cells projected onto the tree, including pseudotime and branch.
##' @export
project.cells.onto.ppt <- function(r,emb=NULL,n.mapping=1) {
  if (is.null(r$root)) { stop("Assign root first") }

  g <- graph.adjacency(r$B,weighted=TRUE,mode="undirected")

  df.list <- pblapply(1:n.mapping,function(nm){
    #print(paste("mapping",nm))
    # assign nearest principal point for each cell
    if (nm > 1){
      rrm = apply(r$R,1,function(v){sample(1:length(v),size=1,prob=v/sum(v))})
    }else{
      rrm <- apply(r$R,1,which.max)
    }

    # idenfity edge onto which each cell lies
    df <- do.call(rbind,lapply(1:ncol(r$R),function(v) {
      vcells <- which(rrm==v);
      if(length(vcells)>0) {
        # determine which edge the cells belong to neighboring PPs
        nv <- as.integer(neighborhood(g,1,nodes=c(v))[[1]])
        nvd <- shortest.paths(g,v,nv)
        spi <- apply(r$R[vcells,nv[-1],drop=FALSE],1,which.max)+1
        ndf <- data.frame(cell=vcells,v0=v,v1=nv[spi],d=nvd[spi])

        p0 <- r$R[vcells,v]
        p1 <- unlist(lapply(1:length(vcells),function(i) r$R[vcells[i],ndf$v1[i]] ))
        alpha <- runif(length(vcells))
        f <- abs( (sqrt(alpha*p1^2+(1-alpha)*p0^2)-p0)/(p1-p0) )
        ndf$t <- r$pp.info[ndf$v0,]$time+(r$pp.info[ndf$v1,]$time-r$pp.info[ndf$v0,]$time)*alpha
        ndf$seg <- r$pp.info[ndf$v0,]$seg
        ndf$color <- r$pp.info[ndf$v0,]$color

        ndf
      } else {
        return(NULL);
      }
    }))
    df$edge <- apply(df,1,function(x) paste(sort(as.numeric(x[c(2,3)])),collapse="|"))
    df <- df[order(df$t,decreasing=FALSE),]

    ### assign data from ndf table of z.ensemble1
    #ndf <- z.ensemble1[[nm]]$ndf[,1:5]
    #ndf[,6:8] <-  z.ensemble1[[nm]]$cell.pseudotime[match(z.ensemble1[[nm]]$ndf$cell,z.ensemble1[[nm]]$cell.pseudotime$cell),2:4]
    #colnames(ndf)[6] <- "t"
    #rownames(ndf) <- nc.cells[ndf$cell]
    #df <- ndf
    #df <- df[order(df$t,decreasing=FALSE),]

    return(df)
  })

  # generate graph of cells and PPs for each mapping
  img.list <- pblapply(df.list,function(df){
    img <- g#graph.adjacency(r$B,weighted=TRUE,mode="undirected")
    img <- set.vertex.attribute(img,"type",value="pp")
    for(e in unique(df$edge)){
      ii <- which(df$edge==e);
      vc <- as.integer(strsplit(e,'\\|')[[1]]);
      imin <- which.min(r$pp.info$time[vc])
      #print(imin)
      #imin <- 1
      #print(c(imin,3-imin))
      # insert the cells
      if (imin==1){
        img <- add_vertices(img,length(ii),type="cell",name=paste('c',df[ii,]$cell,sep=''))
      }else{
        img <- add_vertices(img,length(ii),type="cell",name=paste('c',rev(df[ii,]$cell),sep=''))
      }
      tw <- 1-E(g,path=c(vc[1],vc[2]))$weight
      img <- delete_edges(img,e)
      if (imin==1){
        img <- add_edges(img,c(vc[1],rep(paste0('c',df$cell[ii]),each=2),vc[2]), weight=1-tw*diff(c(0,df$t[ii],1)) )
      }else{
        img <- add_edges(img,c(vc[1],rep(paste0('c',rev(df$cell[ii])),each=2),vc[2]), weight=1-tw*diff(c(0,df$t[ii],1)) )
      }
    }
    return(img)
  })

  if (n.mapping > 1) {
    df.sd <- apply(do.call(cbind,lapply(df.list,function(el)el[rownames(r$R),]$t)),1,sd)
  }else {df.sd <- NA}
  df.summary <- cbind(df.list[[1]],t.sd=df.sd)


  if (!is.null(emb)){
    cols <- adjustcolor(df.summary[rownames(r$R),]$color,0.2); names(cols) <- rownames(r$R)
    plotppt(r,emb,cex.col=cols, tree.col = r$pp.info$color,cex.main=0.5, cex.title=1,cex.tree=1,lwd.tree=1)
  }

  r$cell.summary <- df.summary
  r$cell.info <- df.list
  r$img.list <- img.list
  #r$mg <- mg;
  return(invisible(r))
}


##' Determine a set of genes significantly associated with the tree
##' @param r pptree object
##' @param X expressinon matrix of genes (row) vs cells (column)
##' @param fdr.cut FDR (Benjamini-Hochberg adjustment) cutoff on significance; significance if FDR < fdr.cut
##' @param A.cut cmplitude cutoff on significance; significance if A > A.cut
##' @param st.cut cutoff on stability (fraction of mappings with significant (fdr,A) pair) of association; significance, significance if A > A.cut
##' @param summary show plot of amplitude vs FDR of each gene's association. By default FALSE.
##' @param subtree restrict statistical assesment to a subtree
##' @param fdr.method a method to adjust for multiple testing. Default - Bonferroni. Alternatively, "BH" can be used.
##' @return modified pptree object with a new field r$stat.association that includes pvalue, amplitude, fdr, stability and siginificane (TRUE/FALSE) of gene associations
##' @export
test.associated.genes <- function(r,X,n.map=1,n.cores=(parallel::detectCores()/2),spline.df=3,fdr.cut=1e-4,A.cut=1,st.cut=0.8,summary=FALSE,subtree=NA,fdr.method=NULL, ...) {
  if (is.null(r$root)) {stop("assign root first")}
  if (is.null(r$cell.summary) | is.null(r$cell.info)) {stop("project cells onto the tree first")}
  X <- X[,intersect(colnames(X),rownames(r$cell.summary))]
  if (sum(!colnames(X) %in% rownames(r$cell.summary)) > 0) {stop( paste("Expression matrix X contains cells not mapped onto the tree, e.g. cell",colnames(X)[!colnames(X) %in% rownames(r$cell.summary)][1]) )}
  if (n.map < 0 | n.map > length(r$cell.info)) {stop("n.map should be more than 0 and less than number of mappings")}

  genes <- rownames(X)
  subseg <- unique(r$cell.summary$seg);
  if (!is.na(subtree)) {subseg <- subtree$segs}
  # for every gene
  gtl <- lapply(1:n.map,function(ix){
    print(paste("mapping",ix,"of",n.map))
    if (n.map==1){ inf <- r$cell.summary}else{
      inf <- r$cell.info[[ix]]
    }
    gt <- do.call(rbind,mclapply(genes,function(gene) {
      #sdf <- inf; sdf$exp <- X[gene,rownames(inf)]
      sdf <- inf[inf$seg%in%subseg,]; sdf$exp <- X[gene,rownames(sdf)]#[inf$seg%in%subseg]
      # time-based models
      mdl <- tapply(1:nrow(sdf),as.factor(sdf$seg),function(ii) {
        # TODO: adjust df according to branch length?
        m <- mgcv::gam(exp~s(t,k=spline.df),data=sdf[ii,],familly=gaussian())
        rl <- list(d=deviance(m),df=df.residual(m))
        rl$p <- predict(m);
        return(rl)
      })
      mdf <- data.frame(do.call(rbind,lapply(mdl,function(x) c(d=x$d,df=x$df))))
      # background model
      odf <- sum(mdf$df)-nrow(mdf); # correct for multiple segments
      m0 <- mgcv::gam(exp~1,data=sdf,familly=gaussian())
      if (sum(mdf$d)==0){ fstat <- 0}else{
        fstat <- (deviance(m0) - sum(mdf$d))/(df.residual(m0)-odf)/(sum(mdf$d)/odf)
      }
      pval <-  pf(fstat,df.residual(m0)-odf,odf,lower.tail = FALSE);#1-pf(fstat,df.residual(m0)-odf,odf,lower.tail = T);
      pr <- unlist(lapply(mdl,function(x) x$p))
      return(c(pval=pval,A=max(pr)-min(pr)))
    },mc.cores=n.cores,mc.preschedule=T))
    gt <- data.frame(gt); rownames(gt) <- genes
    if (is.null(fdr.method)) {
      gt$fdr <- p.adjust(gt$pval)
    }else{
      gt$fdr <- p.adjust(gt$pval,method=fdr.method)
    }
    gt
  })

  stat.association <- data.frame(cbind( apply(do.call(cbind,lapply(gtl,function(gt)gt$pval)),1,median),
                                        apply(do.call(cbind,lapply(gtl,function(gt)gt$A)),1,median),
                                        apply(do.call(cbind,lapply(gtl,function(gt)gt$fdr)),1,median),
                                        apply(do.call(cbind,lapply(gtl,function(gt) gt$fdr < fdr.cut & gt$A > A.cut )),1,sum)/length(gtl)
  ))
  rownames(stat.association) <- genes; colnames(stat.association) <- c("pval","A","fdr","st")
  stat.association$sign <- stat.association$fdr < fdr.cut & stat.association$A > A.cut & stat.association$st > st.cut

  # plot amplitude vs FDR and color genes that were idenfitied as significantly associated with the tree
  if (summary==TRUE){
    par(mfrow=c(1,1),mar=c(4.5,4.5,1,1))
    plot(stat.association$A,stat.association$fdr,xlab="Amplitude",ylab="FDR, log",log="y",pch=19,cex=0.5,
         col=adjustcolor( ifelse(stat.association$sign==TRUE,"red","black") ,0.4),cex.lab=1.5)
    legend("bottomleft", legend=c( paste("DE,",sum(stat.association$sign)), paste("non-DE,",sum(!stat.association$sign))),
           col=c("red", "black"), bty="n",pch=19,cex=1,pt.cex=1)
  }
  if (is.na(subtree)){
    r$stat.association <- stat.association
    return(r)
  }else{
    return(stat.association)
  }
}

##' Model gene expression levels as a function of tree positions.
##' @param r pptree object
##' @param X expressinon matrix of genes (rows) vs cells (columns)
##' @param n.map number of probabilistic cell-to-tree mappings to use
##' @param method method of modeling. Currently only splines with option 'ts' are supported.
##' @param knn use expression averaging among knn cells
##' @param gamma stringency of penalty.
##' @return modified pptree object with new fields r$fit.list, r$fit.summary and r$fit.pattern. r$fit.pattern contains matrix of fitted gene expression levels
##' @export
fit.associated.genes <- function(r,X,n.map=1,n.cores=parallel::detectCores()/2,method="ts",knn=1,gamma=1.5) {
  if (is.null(r$root)) {stop("assign root first")}
  if (is.null(r$cell.summary) | is.null(r$cell.info)) {stop("project cells onto the tree first")}
  X <- X[,intersect(colnames(X),rownames(r$cell.summary))]
  if (sum(!colnames(X) %in% rownames(r$cell.summary)) > 0) {stop( paste("Expression matrix X contains cells not mapped onto the tree, e.g. cell",colnames(X)[!colnames(X) %in% rownames(r$cell.summary)][1]) )}
  if (n.map < 0 | n.map > length(r$cell.info)) {stop("n.map should be more than 0 and less than number of mappings")}
  if ( is.null(r$stat.association) ) {stop("identiy significantly associated genes using test.associated.genes()")}


  #gtl <- lapply(1:n.map,function(ix){
  #  print(paste("mapping",ix,"of",n.map))
  #  if (n.map==1){ inf <- r$cell.summary}else{
  #    inf <- r$cell.info[[ix]]
  #  }

  if (method=="ts"){
    gtl <- fit.ts(r,X,n.map,n.cores,gamma,knn)
  }else if (method=="sf"){
    gtl <- t.fit.sf(r,X,n.map,n.cores,gamma)
  }else if (method=="av"){
    gtl <- t.fit.av(r,X,n.map,n.cores)
  }else{stop("please choose correct method name")}
  #})


  ft.summary <- matrix(0,nrow=nrow(gtl[[1]]),ncol=ncol(gtl[[1]])); rownames(ft.summary) <- rownames(gtl[[1]]); colnames(ft.summary) <- colnames(gtl[[1]])
  if (length(gtl)>=1){
    for (k in 1:length(gtl)){
      #indx <- unlist(lapply(1:nrow(r$cell.summary),function(i) {
      #  #ind <- rownames(r$cell.info[[k]])[r$cell.info[[k]]$seg==r$cell.summary$seg[i]]
      #  #ind[which.min(abs(r$cell.info[[k]][ind,]$t-r$cell.summary$t[i]))]
      #  ind <- rownames(r$cell.summary)[r$cell.summary$seg==r$cell.summary$seg[i]]
      #  ind[which.min(abs(r$cell.summary[ind,]$t-r$cell.summary$t[i]))]
      #}))
      ft.summary <- ft.summary + gtl[[k]]#[,indx]
    }
  }
  ft.summary <- ft.summary/length(gtl)
  #colnames(ft.summary) <- rownames(r$cell.summary)
  r$fit.list <- gtl
  r$fit.summary <- ft.summary
  r$fit.pattern <- classify.genes(r)
  print(table(r$fit.pattern))

  return(r)
}


##' Model gene expression levels as a brancing spline function of tree positions.
##' @param r pptree object
##' @param X expressinon matrix of genes (rows) vs cells (columns)
##' @param n.map number of probabilistic cell-to-tree mappings to use
##' @param knn use expression averaging among knn cells
##' @param gamma stringency of penalty.
##' @return matrix of fitted gene expression levels to the tree
##' @export
fit.ts <- function(r,X,n.map,n.cores=parallel::detectCores()/2,gamma=1.5,knn=1) {
  ix <- 1
  img = r$img.list[[ix]];
  root = r$root
  tips = r$tips[r$tips != root]
  branches.ll = do.call(rbind,lapply(tips, function(tip){
    b = get.shortest.paths(img,from=as.character(root),to=as.character(tip))$vpath[[1]]$name
    b = b[grepl("^c",b)]
    ind <- paste('c',r$cell.info[[ix]]$cell,sep="") %in% b
    cbind( ids=rownames(r$cell.info[[ix]])[ind], r$cell.info[[ix]][ind,],branch=rep( which(tips==tip),length(b)) )
  }))
  # calculate knn for each vertex along the tree
  for (v in r$pp.info$PP){img <- delete_vertices(img,as.character(v))}
  dst.tree <- distances(img,v=V(img),to=V(img));
  dst.tree <- dst.tree[ paste("c",r$cell.summary$cell,sep=""),paste("c",r$cell.summary$cell,sep="") ]
  rownames(dst.tree) <- colnames(dst.tree) <- rownames(r$cell.summary)
  dst.tree[dst.tree <= knn] <- 1; dst.tree[dst.tree > knn] <- 0

  gtl <- lapply(1:n.map,function(ix){
    print(paste("fit gene expression for mapping",ix))
    img = r$img.list[[ix]];
    root = r$root
    tips = r$tips[r$tips != root]
    branches = do.call(rbind,lapply(tips, function(tip){
      b = get.shortest.paths(img,from=as.character(root),to=as.character(tip))$vpath[[1]]$name
      b = b[grepl("^c",b)]
      ind <- paste('c',r$cell.info[[ix]]$cell,sep="") %in% b
      cbind( ids=rownames(r$cell.info[[ix]])[ind], r$cell.info[[ix]][ind,],branch=rep( which(tips==tip),length(b)) )
    }))

    #branches.ll <- branches
    genes <- intersect(rownames(X),rownames(r$stat.association)[r$stat.association$sign])
    gt <- do.call(rbind,mclapply(genes,function(gene) {
      expr.fitted <- unlist(lapply(unique(branches$branch),function(br){
        branches1 <- branches[branches$branch==br,]
        expr <- X[gene,as.character(branches1$ids)]
        #gene.fit1 = gam( expr ~ s( branches1$time,k=length(branches1$time),bs="ts"),knots=list(branches1$time) )
        tt <- branches1$t
        #tt <- 1:length(tt)
        gene.fit1 = mgcv::gam( expr ~ s(tt,bs="ts"),gamma=gamma)
        #ggplot()+geom_point(aes(tt,expr))+geom_line(aes(tt,gene.fit1$fitted.values))
        td <- data.frame(matrix(branches.ll[branches.ll$branch==br,]$t,nrow=sum(branches.ll$branch==br)));
        rownames(td) <- branches.ll[branches.ll$branch==br,]$ids; colnames(td) <- "tt"
        predict(gene.fit1,td )
      }))

      # old version - averaging along shared branches
      #for( cell in names(which(table(branches.ll$ids) > 1))){
      #  expr.fitted[branches.ll$ids==cell] <- mean(expr.fitted[branches.ll$ids==cell])
      #}

      # new version - knn smoothing, where knns are estimated along the tree.
      expr.fitted <- (dst.tree[names(expr.fitted),names(expr.fitted)] %*% expr.fitted) / (apply(dst.tree[names(expr.fitted),names(expr.fitted)],1,sum))
      expr.fitted <- expr.fitted[,1]
      return(expr.fitted[!duplicated(names(expr.fitted))])
    },mc.cores = n.cores))
    rownames(gt) <- genes
    return(gt)
  })
  return(gtl)
}



##' Classify tree-associated genes
##'
##' Tree-associated genes are classified in branch-monotonous, transiently expressed and having complex patterns.
##' @param r tree
##' @param X expressinon matrix of genes (rows) vs cell (columns)
##' @param cutoff expression in local optimum should be higher/lower than both terminal branch values by cutoff.
##' @return vector of predicted classification for fitted genes.
##' @export
classify.genes <- function(r,n.cores=parallel::detectCores()/2,cutoff=0.2) {
  if (is.null(r$fit.summary)) {stop("fit gene expression to the tree first")}
  a <- do.call(cbind,lapply(unique(r$cell.summary$seg),function(seg){
    seg.summary <- r$cell.summary[r$cell.summary$seg==seg,]
    tt <- r$fit.summary[,rownames(seg.summary)][,order(seg.summary$t)]
    # calculate number of inner local optima
    apply(tt,1,function(x) {
      res <- loc.opt(x)
      if ( sum(!is.na(res))==0 ){0}else{nrow(res)}
    })
  }))
  apply(a,1,function(v){
    if (sum(v)==0) {return("branch-monotonous")}else
      if (sum(v)==1) {return("transiently expressed")}else
        if (sum(v)>1) {return("complex patterns")}
  })
}

##' Identify all local optima for a time series data
##' @name loc.opt
##' @param series - time series data
##' @param cutoff - expression in local optimum should be on cutoff higher/lower than nearby local optima. This parameter allows to eliminate small local optimas that are likely artifacts
##' @return data frame containing type of local optima (min/max) and time index.
##' @export
loc.opt <- function(series,cutoff=0.1){
  dx <- diff(series)
  cand <- (-dx[1:(length(dx)-1)]*dx[2:length(dx)]) > 0
  # remove multiple rupture-related optima
  cand[1:(length(cand)-1)][cand[1:(length(cand)-1)]&cand[2:length(cand)]] <- FALSE
  if (sum(cand)>0){
    cand <- c(TRUE,cand,TRUE)
    ds <- diff(series[cand])
    opt.type <- unlist(lapply(1:(sum(cand)-2),function(i){
      if (ds[i] > cutoff & (-ds[i+1]) > cutoff ) {
        "max"
      }else if (ds[i] < -cutoff & (-ds[i+1]) < -cutoff ){
        "min"
      }else{
        NA
      }
    }))
    if ( sum(!is.na(opt.type))>0  ){
      opt.inf <- data.frame(cbind( opt.type[!is.na(opt.type)],as.numeric(which(cand))[2:(sum(cand)-1)][!is.na(opt.type)]),stringsAsFactors=FALSE)
      colnames(opt.inf) <- c("type","index"); opt.inf$index <- as.numeric(opt.inf$index)
      return(opt.inf)
    }
  }
  return(NA)
}


##' Visualize branching trajectories of a particular gene.
##' @param r pptree object
##' @param gene gene name
##' @param X matrix with a single row containing a gene expression levels (could be a vector of gene's expression). Columns of X reflect gene names.
##' @param cex.cell size of cells
##' @param cex.lab size of axis titles
##' @param cex.axis size of axis labels
##' @param cex.main size of title showing a gene name
##' @param lwd.t1 width of the main branching trajectory
##' @param lwd.t2 width of ensemble trajectories, typically thiner than that of main trajectory.
##' @param lwd.erbar width of error bars for uncertainty of cell pseudotime assignment
##' @param subtree visualise trajectory along a given subtree
##' @export
visualise.trajectory = function(r,gene,X,cex.cell=0.3,cex.lab=2,cex.axis=1.5,cex.main=1,lwd.erbar=0.0,lwd.t1=3,lwd.t2=0.2,switch.point=NA,subtree=NA){

  if (is.null(dim(X))){
    Xgene <- X
  }else{
    if ( gene %in% rownames(X) == FALSE ) {stop("gene is not in matrix X")}
    Xgene <- X[gene,]
  }
  Xgene <- Xgene[intersect(names(Xgene),rownames(r$cell.summary))]
  if ( sum(!names(Xgene)%in%rownames(r$cell.summary)) > 0 ) {stop("matrix/vector X does not contain some cells used to recostruct tree")}


  segs <- unique(r$cell.summary$seg)
  # restrict considered segments to subtree if given
  if (!is.na(subtree)){
    segs <- intersect(segs,subtree$seg)
  }

  par(mar=c(5,5,3,1))
  # draw cells
  ind <- r$cell.summary$seg%in%segs
  plot(r$cell.summary$t[ind],Xgene[rownames(r$cell.summary)][ind],type = "n",
       xlab="pseudotime",ylab="expression",cex.axis=cex.axis,cex.lab=cex.lab,main=gene,font.main=3,cex.main=cex.main)
  grid(5,5,lwd=1.5)
  points(r$cell.summary$t[ind],Xgene[rownames(r$cell.summary)][ind],col=adjustcolor(r$cell.summary$color[ind],0.5),pch=19,cex=cex.cell)
  # draw error bars of pseudotime uncertainty if given
  if ( sum(!is.na(r$cell.summary$t.sd))>0 ){
    segments( r$cell.summary$t[ind]-r$cell.summary$t.sd[ind], Xgene[rownames(r$cell.summary)][ind], r$cell.summary$t[ind]+r$cell.summary$t.sd[ind], y1 = Xgene[rownames(r$cell.summary)][ind],
              col=adjustcolor(r$cell.summary$color[ind],0.1),lwd=lwd.erbar)
  }
  # draw ensemble of sampled trajectries if given
  if (length(r$fit.list)>1){
    for (j in 2:length(r$fit.list)){
      for(seg in segs ){
        #ind <- r$cell.info[[j]]$seg == seg
        #t.ord <- order(r$cell.info[[j]]$t[ind])
        #lines(r$cell.info[[j]]$t[ind][t.ord],r$fit.list[[j]][gene,rownames(r$cell.info[[j]])][ind][t.ord],
        #      col=adjustcolor(r$cell.info[[j]]$color[ind][t.ord],0.4),lwd=lwd.t2)

        ind <- r$cell.summary$seg == seg
        t.ord <- order(r$cell.summary$t[ind])
        lines(r$cell.summary$t[ind][t.ord],r$fit.list[[j]][gene,rownames(r$cell.summary)][ind][t.ord],
              col=adjustcolor(r$cell.summary$color[ind][t.ord],0.4),lwd=lwd.t2)

      }
    }
  }
  # draw likelihood trajectory
  for(seg in segs ){
    ind <- r$cell.summary$seg == seg
    t.ord <- order(r$cell.summary$t[ind])
    lines(r$cell.summary$t[ind][t.ord],r$fit.summary[gene,rownames(r$cell.summary)][ind][t.ord],
          col=r$cell.summary$color[ind][t.ord],lwd=lwd.t1)
  }
  if (!is.na(switch.point)){
    abline(v=switch.point,lty=1,lwd=3,col=adjustcolor("black",0.5))
  }
  # connect boundary cells from different branches
  g <- r$img.list[[1]]
  for (seg in segs){

    ind <- r$cell.summary$seg==seg
    c2.name <- rownames(r$cell.summary[ind,])[which.min(r$cell.summary$t[ind])]
    c2 <- r$cell.summary$cell[ind][which.min(r$cell.summary$t[ind])]
    c2.seg <- r$cell.summary$seg[ind][which.min(r$cell.summary$t[ind])]

    c2.path <- names(shortest_paths(g,r$root,paste("c",c2,sep="") )$vpath[[1]])
    c2.path <- c2.path[unlist(lapply(1:length(c2.path),function(i) grepl("c",c2.path[i])))]
    c2.path <- as.numeric(unlist(lapply(strsplit(c2.path,"c"),function(x)x[2])))
    ind <- r$cell.summary$cell %in% c2.path & r$cell.summary$cell != c2 #& !(r$cell.summary$seg %in% r$cell.summary[c2.name,]$seg)

    if (sum(ind)>0){
      c1.name <- rownames(r$cell.summary[ind,])[which.max(r$cell.summary$t[ind])]
      segments(r$cell.summary[c(c1.name),]$t,r$fit.summary[gene,c(c1.name)],r$cell.summary[c(c2.name),]$t,r$fit.summary[gene,c(c2.name)],
               col=r$cell.summary[c2.name,]$color,lwd=lwd.t1)
    }
  }

}


##' Visualize clusters of genes using heatmap and consensus tree-projected pattern.
##' @param r pptree object
##' @param emb cells embedding
##' @param clust a vector of cluster numbers named by genes
##' @param n.best show n.best the most representative genes on the heatmap for each cluster
##' @param best.method use method to select the most representative genes. Current options: "pca" selects genes with the highest loading on pc1 component reconstructed using genes from a cluster, "cor" selects genes that have the highest average correlation with other genes from a cluster.
##' @param cex.gene size of gene names
##' @param cex.cell size of cells on embedding
##' @param cex.tree width of line of tree on embedding
##' @param reclust whether to reorder cells inside individual clusters on heatmap according to hierarchical clustering using Ward linkage and 1-Pearson as a distance between genes. By default is FALSE.
##' @param subtree visualize clusters for a given subtree
##' @export
visualise.clusters <-function(r,emb,clust=NA,clust.n=5,n.best=4,best.method="cor",cex.gene=1,cex.cell=0.1,cex.tree=2,subtree=NA, reclust=TRUE){


  if ( !is.na(clust) & sum(!names(clust)%in%rownames(r$fit.summary))>0) {stop( paste("Expression is not fitted for",sum(!names(clust)%in%rownames(r$fit.summary)),"genes" ))}
  if (best.method!="pca" & best.method!="cor") {stop(paste("incorrect best.method option",best.method) )}
  tseg <- unlist(lapply( unique(r$cell.summary$seg),function(seg)mean(r$cell.summary$t[r$cell.summary$seg==seg]))); names(tseg) <-  unique(r$cell.summary$seg)
  tseg <- tseg[as.character(r$cell.summary$seg)]

  gns <- rownames(ppt$fit.summary)
  if (!is.na(clust)){gns <- names(clust)}
  emat <- r$fit.summary[gns,rownames(r$cell.summary)][,order(tseg,r$cell.summary$t)]
  emat <- t(apply(emat,1,function(x) (x-mean(x))/sd(x) ))
  cols <- r$cell.summary$col[order(tseg,r$cell.summary$t)]
  subcells = TRUE; if (!is.na(subtree)){subcells <- r$cell.summary$seg[order(tseg,r$cell.summary$t)]%in%subtree$seg}

  # cluster genes if necessary
  if (is.na(clust)){
    gns <- rownames(emat)#names(clust)[clust==cln]
    dst.cor <- 1-cor(t(emat[gns,]))
    hcl <- hclust(as.dist(dst.cor),method="ward.D")
    clust <- cutree(hcl,clust.n)
  }

  k <- length(unique(clust))
  genes.show <- unlist(lapply(1:k,function(i){
    n <- n.best; if ( sum(clust==i) < n) {n <- sum(clust==i)}
    if (best.method=="pca"){
      pr <- pca(t(emat[clust==i,]),center = TRUE, scale = "uv")
      pr.best <- rep(i,n); names(pr.best) <- names(sort(pr@loadings[,1],decreasing = T))[1:n]
      return(pr.best)
    }else if (best.method=="cor"){
      cr <- cor(t(emat[clust==i,]))
      cr.best <- rep(i,n); names(cr.best) <- names(sort(apply(cr,1,mean),decreasing = TRUE))[1:n]
      return(cr.best)
    }
  }))

  nf <- layout( matrix(unlist(lapply(1:k,function(i) 5*(i-1)+c(1,2,3,1,4,5))),2*k,3, byrow=T),respect = T,width=c(1,1,0.1),heights=rep(c(0.1,1),k) )
  #layout.show(nf)
  for (cln in 1:k){
    # recluster genes inside module if necessary
    gns <- names(clust)[clust==cln]
    if (reclust==TRUE){
      dst.cor <- 1-cor(t(emat[gns,]))
      hclust.cor <- hclust(as.dist(dst.cor),method="ward.D")
      gns <- gns[hclust.cor$order]
    }

    # draw cluster-wise pattern
    par(mar=c(0.3,0.1,0.0,0.2))
    plotppt(r,emb,pattern.cell = apply(emat[clust==cln,],2,mean),cex.main=cex.cell,cex.tree = cex.tree,lwd.tree = 0.1,subtree=subtree)
    # draw color-scheme for branches
    #par(mar=c(0.0,0.2,0.1,2))
    par(mar=c(0.0,0.0,0.0,0))
    col.ind <- 1:length(unique(cols)); names(col.ind) = unique(cols)
    image( t(rbind( col.ind[cols[subcells]] )),axes=FALSE,col=(unique(cols[subcells])) )
    box()

    par(mar=c(0.0,0.0,0.0,0))
    plot(0.2,0.2,ylim=c(0.05,0.95),xlim=c(0,1),xaxt='n',yaxt='n',pch='',ylab='',xlab='',bty='n')

    #par(mar=c(0.2,0.2,0.0,2))
    par(mar=c(0.3,0.0,0.0,0))
    image( t(emat[gns,subcells]),axes=FALSE,col=colorRampPalette(c("blue","grey80","red"))(n = 60))
    #axis( 4, at=seq(0,1,length.out=sum(clust==cln)),col.axis="black", labels=gns,hadj=0.1,xaxt="s",cex.axis=1.5,font = 3,las= 1,tick=FALSE)
    box()

    gns[! gns %in% names(genes.show)[genes.show==cln] ] <- ""
    ### calculate coordinates of genes.show with QP
    coord <- which( names(clust)[clust==cln] %in% names(genes.show)[genes.show==cln] )/sum(clust==cln)
    del <- 1/(sum(genes.show==cln))#0.1
    Dmat <- diag(1,length(coord),length(coord))
    dvec <- rep(0,length(coord))
    Amat <- matrix(0,nrow= 3*length(coord)-1,ncol=length(coord)); bvec = rep(0,3*length(coord)-1)
    for (i in 1:(length(coord)-1)){Amat[i,i] <- -1; Amat[i,i+1] <- 1; bvec[i] <- del - (coord[i+1]-coord[i])}
    for (i in 1:(length(coord))){j <- i+length(coord)-1; Amat[j,i] <- 1; bvec[j] <- -coord[i]+0 }
    for (i in 1:(length(coord))){j <- i+2*length(coord)-1; Amat[j,i] <- -1; bvec[j] <- coord[i]-1}
    qp = solve.QP(Dmat, dvec, t(Amat), bvec, meq=0, factorized=FALSE)
    coord_new = qp$solution + coord
    par(mar=c(0.3,0,0,0))
    plot(0.2,0.2,ylim=c(0.0,1),xlim=c(0,1),xaxt='n',yaxt='n',pch='',ylab='',xlab='',bty='n')
    axis(side = 4, at = coord_new,lwd=0.0,lwd.ticks=0,font=3,cex.axis=cex.gene,labels=gns[gns!=""],tck=0.0,hadj=0.0,line=-0.9,las=1)
    for (i in 1:length(coord)){
      arrows( 0,coord[i],1,coord_new[i],length=0.0,lwd=0.7 )
    }
    ###
  }

}



##' Determine genes differentially upregulated after bifurcation point
##' @param r pptree object
##' @param mat expression matrix of genes (rows) and cells (columnts)
##' @param root a principal point of fork root
##' @param leaves vector of two principal points of fork leaves
##' @param genes optional set of genes to estimate association with fork
##' @param n.mapping number of probabilistic cell-to-tree projections to use for robustness
##' @param n.mapping.up number of probabilistic cell-to-tree projections to estimate the amount of upregulation relative to progenitor branch
##' @return summary statistics of size effect and p-value of association with bifurcaiton fork.
##' @export
test.fork.genes <- function(r,mat,matw=NULL,root,leaves,genes=rownames(mat),n.mapping=1,n.mapping.up=1,n.cores=parallel::detectCores()/2) {
  g <- graph.adjacency(r$B>0,mode="undirected")
  vpath = get.shortest.paths(g,root,leaves)
  interPP = intersection(vpath$vpath[[1]],vpath$vpath[[2]])
  which.max(r$pp.info[interPP,]$time)
  vpath = get.shortest.paths(g, r$pp.info[interPP,]$PP[which.max(r$pp.info[interPP,]$time)],leaves)
  cat("testing differential expression between branches ..");cat("\n")
  gtll <- lapply( 1:n.mapping,function(nm){
    cat("mapping ");cat(nm);cat("\n")
    cell.info <- r$cell.info[[nm]]
    brcells = do.call(rbind,lapply( 1:length(vpath$vpath), function(i){
      x=vpath$vpath[[i]]
      segs = as.numeric(names(table(r$pp.info[x,]$seg))[table(r$pp.info[x,]$seg)>1])
      return(cbind(cell.info[cell.info$seg %in% segs,],i))
    }))
    # for every gene
    gtl <- do.call(rbind,mclapply(genes,function(gene) {
      brcells$exp <- mat[gene,rownames(brcells)]
      if (is.null(matw)) {brcells$w = 1
      }else {brcells$w <- matw[gene,r$cells][as.integer(gsub("c","",brcells$node))]}
      # time-based models
      m <- mgcv::gam(exp ~ s(t)+s(t,by=as.factor(i))+as.factor(i),data=brcells,familly=gaussian(),weights=brcells$w)
      return( c(mean(brcells$exp[brcells$i==1])-mean(brcells$exp[brcells$i==2]) , min(summary(m)$p.pv[2]) ) )

      #m <- mgcv::gam(exp ~ s(t)+as.factor(i),data=brcells,familly=gaussian(),weights=brcells$w)
      #return( c(mean(brcells$exp[brcells$i==2])-mean(brcells$exp[brcells$i==1]) , min(summary(m)$s.pv[2:3]) ) )
    },mc.cores=n.cores,mc.preschedule=T));
    colnames(gtl) = c("effect","p"); rownames(gtl) = genes; gtl = as.data.frame(gtl)
    return(gtl)
  })

  effect = do.call(cbind,lapply(gtll,function(gtl) gtl$effect ))
  if (length(gtll) > 1) {effect <- apply(effect,1,median)}
  pval = do.call(cbind,lapply(gtll,function(gtl) gtl$p ))
  if (length(gtll) > 1) {pval <- apply(pval,1,median)}
  fdr = do.call(cbind,lapply(gtll,function(gtl) p.adjust(gtl$p,"BH") ))
  if (length(gtll) > 1) {fdr <- apply(fdr,1,median)}
  st = do.call(cbind,lapply(gtll,function(gtl) gtl$p < 5e-2 ))
  if (length(gtll) > 1) {st <- apply(st,1,mean)}
  stf = do.call(cbind,lapply(gtll,function(gtl) p.adjust(gtl$p,"BH") < 5e-2 ))
  if (length(gtll) > 1) {stf <- apply(stf,1,mean)}

  ### here add a code that estimates the amount of upregulation relative to progenitor branch.
  cat("testing upregulation in derivative relative to progenitor branch ..");cat("\n")
  # n.mapping.up
  eu <- do.call(cbind,lapply(leaves[1:2],function(leave){
    segs = extract.subtree(ppt,c(root,leave))
    posit = do.call(rbind,(mclapply(genes,function(gene){
      eu <- do.call(rbind,lapply(1:n.mapping.up,function(j){
        cells = rownames(r$cell.info[[j]])[r$cell.info[[j]]$seg %in% segs$segs]
        ft = lm( mat[gene,cells] ~ r$cell.info[[j]][cells,]$t   )
        return( c(ft$coefficients[2],summary(ft)$coefficients[2,4] ) )
      }))
      if (n.mapping.up > 1) {eu <- apply(eu,2,median)}
      return(eu)
    },mc.cores = n.cores,mc.preschedule = TRUE)))
  }))
  colnames(eu) <- c("pd1.a","pd1.p","pd2.a","pd2.p")

  res <- as.data.frame(cbind(effect = effect, p = pval, fdr = fdr, st = st,stf = stf))
  colnames(res) <- c("effect","p","fdr","st","stf")
  rownames(res) <- genes
  res <- cbind(res,eu)
  return(res)
}


##' Assign genes differentially expressed between two post-bifurcation branches
##' @param fork.de statistics on expression differences betwee post-bifurcation branches, return of test.fork.genes
##' @param stf.cut fraction of projections when gene passed fdr < 0.05
##' @param effect.b1 expression differences to call gene as differentially upregulated at branch 1
##' @param effect.b2 expression differences to call gene as differentially upregulated at branch 2
##' @param pd.a  minium expression increase at derivative compared to progenitor branches to call gene as branch-specific
##' @param pd.p p-value of expression changes of derivative compared to progenitor branches to call gene as branch-specific
##' @return table fork.de  with added column stat, which classfies genes in branch-specifc (1 or 2) and non-branch-specific (0)
##' @export
branch.specific.genes <- function(fork.de,stf.cut = 0.7, effect.b1 = 0.1,effect.b2 = 0.3, pd.a = 0, pd.p = 5e-2){
  ind <- fork.de$stf >= stf.cut & fork.de$effect  > effect.b1 & fork.de$pd1.a > pd.a & fork.de$pd1.p < pd.p
  gns1 <- rownames(fork.de)[ind]

  ind <- fork.de$stf >= stf.cut & fork.de$effect  < -effect.b2 & fork.de$pd2.a > pd.a & fork.de$pd2.p < pd.p
  gns2 <- rownames(fork.de)[ind]

  state <- rep(0,nrow(fork.de)); names(state) <- rownames(fork.de)
  state[gns1] <- 1
  state[gns2] <- 2

  return(cbind(fork.de,state))
}


##' Estimate optimum of expression and time of activation
##' @param r ppt.tree object
##' @param mat expression matrix
##' @param root root of progenitor branch of bifurcation
##' @param leaves leaves of derivative branches of bifurcation
##' @param genes genes to estimate parameters
##' @param deriv.cutoff a first passage of derivative through cutoff 'deriv.cutoff' to predict activation timing
##' @param gamma gamma parameter in gam function
##' @param n.mapping results are averaged among n.mapping number of probabilsitic cell projections
##' @param n.cores number of cores to use
##' @return per gene timing of optimum and activation
##' @export
activation.statistics <- function(r,mat,root,leave,genes=rownames(mat),deriv.cutoff = 0.015,gamma=1,n.mapping=1,n.cores=parallel::detectCores()/2){
  xx = do.call(rbind,(mclapply(genes,function(gene){
    gres <- do.call(rbind,lapply(1:n.mapping,function(i){
      segs = extract.subtree(ppt,c(root,leave))
      cell.summary <- r$cell.info[[i]]
      cells <- rownames(cell.summary)[cell.summary$seg %in% segs$segs]

      ft = gam( mat[gene,cells] ~ s(cell.summary[cells,]$t),gamma=gamma)
      ord <- order(cell.summary[cells,]$t)
      deriv.n <- ft$fitted.values[ord][-1]-ft$fitted.values[ord][-length(ord)]
      #deriv.d <- r$cell.summary[cells,]$t[-1]-r$cell.summary[cells,]$t[-length(ord)]
      deriv.d <- max(ft$fitted.values[ord]) - min(ft$fitted.values[ord])
      deriv <- deriv.n/deriv.d

      c(cell.summary[cells,]$t[which.max(ft$fitted.values)],
        min(c(cell.summary[cells,]$t[-1][ deriv > deriv.cutoff ],max(cell.summary[cells,]$t))) )
    }))
    c( median(gres[,1]),median(gres[,2]) )
  },mc.cores = n.cores,mc.preschedule = TRUE)))
  rownames(xx) <- genes
  colnames(xx) <- c("optimum","activation")
  return(xx)
}


##' Estimate optimum of expression and time of activation
##' @param r ppt.tree object
##' @param fork.de outcome of test.fork.genes function
##' @param mat expression matrix
##' @param root root of progenitor branch of bifurcation
##' @param leaves leaves of derivative branches of bifurcation
##' @param deriv.cutoff a first passage of derivative through cutoff 'deriv.cutoff' to predict activation timing
##' @param gamma gamma parameter in gam function
##' @param n.mapping results are averaged among n.mapping number of probabilsitic cell projections
##' @param n.cores number of cores to use
##' @return table fork.de with added per gene timing of optimum and activation
##' @export
activation.fork <- function(r,fork.de,mat,root,leaves,deriv.cutoff = 0.015,gamma=1,n.mapping=1,n.cores=parallel::detectCores()/2){
  cat("estimate activation patterns .. branch 1"); cat("\n")
  gg1 <- rownames(fork.de)[fork.de$state==1]
  act1 <- activation.statistics(r,mat,root,leaves[1],genes=gg1,deriv.cutoff = deriv.cutoff,gamma=gamma,n.mapping=n.mapping,n.cores=n.cores)

  cat("estimate activation patterns .. branch 2"); cat("\n")
  gg2 <- rownames(fork.de)[fork.de$state==2]
  act2 <- activation.statistics(r,fpm,root,leaves[2],genes=gg2,deriv.cutoff = deriv.cutoff,gamma=gamma,n.mapping=n.mapping,n.cores=n.cores)

  act <- cbind( rep(NA,nrow(fork.de)),rep(NA,nrow(fork.de)) );
  rownames(act) <- rownames(fork.de); colnames(act) <- colnames(act1)
  act[gg1,] <- act1
  act[gg2,] <- act2
  return( cbind(fork.de,act) )
}


##' Extract subtree of the tree
##' @param r ppt.tree object
##' @param nodes set tips or internal nodes (bifurcations) to extract subtree
##' @return list of segments comprising a subtree.
##' @export
extract.subtree = function(r,nodes){
  g <- graph.adjacency(r$B>0,mode="undirected")
  if ( sum(!nodes%in%V(g)) > 0 ) {stop(paste("the following nodes are not in the tree:",nodes[!nodes%in%V(g)],collapse = " ") )}
  if ( sum( degree(g)==2 & (V(g) %in% nodes) ) > 0 ) {stop( paste("the following nodes are nethier terminal nor fork:",nodes[nodes %in% V(g)[V(g)==2] ],collapse=" ") )}
  vpath = get.shortest.paths(g,nodes[1],nodes)
  v = c()
  for (i in 1:length(vpath$vpath)){
    v=c(v,unlist(vpath$vpath[[i]]))
  }
  v=unique(v)
  segs = r$pp.info$seg[r$pp.info$PP %in% v]
  segs = segs[segs %in% names(table(segs))[table(segs) > 1]]
  #v=v[ r$pp.info[v,]$seg %in% unique(segs)  ] #list( segs = unique(segs), pp = v )
  list( segs = unique(segs) )
}

##' Extract subtree of the tree
##' @param r ppt.tree object
##' @param nodes set tips or internal nodes (bifurcations) to extract subtree
##' @return list of segments comprising a subtree.
##' @export
fork.pt = function(r,root,leaves){
  b1 <- extract.subtree(r,c(root,leaves[1]))
  b2 <- extract.subtree(r,c(root,leaves[2]))
  segs.prog <- intersect(b1$segs,b2$segs)
  segs.b1 <- setdiff(b1$segs,segs.prog)
  segs.b2 <- setdiff(b2$segs,segs.prog)

  time.stat <- c( min(r$pp.info$time[r$pp.info$seg %in% segs.prog]),
                  max(r$pp.info$time[r$pp.info$seg %in% segs.prog]),
                  max(r$pp.info$time[r$pp.info$seg %in% segs.b1]),
                  max(r$pp.info$time[r$pp.info$seg %in% segs.b2])
  )
  names(time.stat) <- c("root","bifurcation","leave 1","leave 2")
  return(time.stat)
}

##' Decompose a number by degrees of 2.
##' @param n number
decompose <- function(n){
  base.binary = c()
  while (n > 0){
    x <- as.integer(log2(n))
    base.binary <- c(base.binary,x)
    n = n - 2^x
  }
  return(base.binary)
}




