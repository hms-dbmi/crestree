library(scde)
library(Cairo)
library(parallel)

process.dataset <- function(dat,name,batch=NULL,k=min(20,ncol(cd)/n.groups),print=FALSE,max.model.plots=50,cd=NULL,varinfo=NULL,knn=NULL,max.adj.var=5) {
  cat("processing ",name,": ");
  if(is.null(cd)) {
    cd <- dat;

    if (print==TRUE){
      CairoPNG(file=paste(name,"reads.per.cell.png",sep="."),width=350,height=350)
      par(mfrow=c(1,1), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 0.9);
      hist(colSums(dat)/1e6,col="wheat",xlab="reads per cell (M)",main="Read counts across cells")
      abline(v=1e5/1e6,lty=2,col=2)
      dev.off()
    }
    table(colSums(dat)>=min.cell.reads)

    if (print==TRUE){
      CairoPNG(file=paste(name,"reads.per.gene.png",sep="."),width=350,height=350)
      par(mfrow=c(1,1), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 0.9);
      hist(log10(rowSums(dat)+1),col="wheat",xlab="reads per gene (log10)",main="Read counts across genes")
      abline(v=log10(min.gene.reads+1),lty=2,col=2)
      dev.off()
    }

    # filter out low-gene cells
    vi <- colSums(cd)>min.cell.reads; table(vi)
    cd <- cd[,vi];

    # remove genes that don't have many reads
    vi <- rowSums(cd)>min.gene.reads; table(vi)
    cd <- cd[vi,];

    # remove genes that are not seen in a sufficient number of cells
    vi <- rowSums(cd>0)>min.gene.cells; table(vi)
    cd <- cd[vi,];
  }
  cat("proceeding with ",nrow(cd)," genes across ",ncol(cd)," cells ");

  if(is.null(knn) || is.null(varinfo)) {
    knn <- knn.error.models(cd,groups=as.factor(rep(name,ncol(cd))),k=k,n.cores=n.cores,min.count.threshold=1,min.nonfailed=min.nonfailed,verbose=0,max.model.plots=max.model.plots,min.size.entries=1e3)
    cat("models ")
    prior <- scde.expression.prior(models=knn,counts=cd,length.out=400,show.plot=F)
    if (print==TRUE) {pdf(file=paste(name,"varnorm.pdf",sep="."),height=4,width=8)}
    varinfo <- pagoda.varnorm(knn,counts=cd,trim=trim/ncol(cd),plot=T,verbose=1,prior=prior,max.adj.var=max.adj.var,weight.df.power=1,batch=batch)
    if (print==TRUE) {dev.off()}
    cat("varinfo ")
    varinfo <- pagoda.subtract.aspect(varinfo,colSums(cd[,rownames(knn)]>0))
  }

  pwpca <- pagoda.pathway.wPCA(varinfo,go.env,n.components=1,n.cores=n.cores,n.internal.shuffles=0,verbose=1)
  cat("pathways ")
  if (print==TRUE) {pdf(file=paste(name,"clvar.pdf",sep="."),height=4,width=8)}
  clpca <- pagoda.gene.clusters(varinfo,trim=(trim+5)/ncol(varinfo$mat),n.clusters=150,n.cores=n.cores,verbose=1,plot=T)
  if (print==TRUE){ dev.off()}
  cat("clusters\n")
  return(list(cd=cd,knn=knn,prior=prior,varinfo=varinfo,pwpca=pwpca,clpca=clpca,batch=batch))
}



postprocess <- function(res,name,colcols=NULL,print=FALSE,n.clusters=10,distance.threshold=0.9, z.score=3,perplexity=20,seed=0,include.aspects=T,top.aspects=50,port=1469) {
  cat("processing ", name,": ")
  pwpca <- res$pwpca; clpca <- res$clpca; cd <- res$cd; prior <- res$prior; knn <- res$knn; varinfo <- res$varinfo;

  tam <- pagoda.top.aspects(pwpca[!grepl("custom",names(pwpca))],clpca,use.oe.scale=F,z.score=z.score,adjust.scores=T)
  grcol <- colorRampPalette(c("white","black"),space="Lab")(1024);
  libsize <- colSums(cd)
  scol <- grcol[round((libsize-min(libsize))/diff(range(libsize))*(length(grcol)-1))+1]

  if(!is.null(varinfo$batch)) {
    bcol<-rainbow(length(levels(varinfo$batch)),v=0.5)[as.integer(varinfo$batch)]
  } else {
    bcol<-NULL;
  }


  x <- pagoda.cluster.cells(tam,varinfo,include.aspects=T,verbose=T,return.details=T)
  hc <- x$clustering;
  cl <- cutree(hc,n.clusters)
  clcol <- rainbow(length(unique(cl)))
  clcol <- clcol[cl]
  require(Cairo)
  if (print==TRUE) {CairoPNG(file=paste(name,"tam.png",sep="."),width=600,height=600)}
  tamr <- pagoda.reduce.loading.redundancy(tam,pwpca,clpca,plot=T,labRow=NA,labCol=NA,box=T,margins=c(0.5,0.5),distance.threshold=0.001,cell.clustering=hc,n.cores=1)
  if (print==TRUE) {dev.off()}

  if (print==TRUE) {CairoPNG(file=paste(name,"tamr.png",sep="."),width=600,height=600)}
  tamr2 <- pagoda.reduce.redundancy(tamr,distance.threshold=distance.threshold,plot=T,cell.clustering=hc,labRow=NA,labCol=NA,box=T,margins=c(0.5,0.5),top=top.aspects,trim=0,weighted.correlation=T,col.cols=rbind(scol))
  if (print==TRUE) {dev.off()}

  if(!is.null(colcols)) { # side colors were supplied
    # enforce column order
    colcols <- colcols[,colnames(varinfo$mat),drop=F]
    colcols <- rbind(clusters=clcol,depth=scol,colcols); # bind cluster colors
    if(!is.null(varinfo$batch)) { colcols <- rbind(colcols,batch=bcol) }
  } else {
    colcols<-rbind(clusters=clcol,depth=scol);
    if(!is.null(varinfo$batch)) { colcols <- rbind(colcols,batch=bcol) }
  }

  if (print==TRUE){
    CairoPNG(file=paste(name,"tamr2.png",sep="."),width=900,height=500)
    pagoda.view.aspects(tamr2,cell.clustering=hc,box=T,labCol=NA,margins=c(0.5,10),row.clustering=NA,col.cols=colcols)
    dev.off();
  }

  #pagoda.view.aspects(tamr2,cell.clustering=hclust(as.dist(1-cor(tamr$xv)),method='ward'),col.cols=rbind(scol),box=T,labCol=NA,margins=c(0.5,25),ColSideColors.unit.vsize=0.05)

  # update the pathway table to include the custom sets
  #tamr2$df <- pagoda.top.aspects(pwpca,clpca,return.table=T,z.score=-Inf);

  if(!include.aspects) {
    x <- pagoda.cluster.cells(tam,varinfo,include.aspects=F,verbose=T,return.details=T)
  }
  library(Rtsne);
  set.seed(seed);
  #tSNE.pagoda <- tsne(x$distance,k=2,perplexity=20,initial_dims=100)
  tSNE.pagoda <- Rtsne(x$distance,is_distance=T,initial_dims=100,perplexity=perplexity)
  if (print==TRUE){
    CairoPNG(file=paste(name,"rtSNE.pagoda.png",sep="."),width=350,height=350)
    par(mfrow=c(1,1), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
    plot(tSNE.pagoda$Y,col=clcol,cex=1,pch=19,xlab="PC1",ylab="PC2")
    dev.off()
  }
  emb <- tSNE.pagoda$Y; rownames(emb)<-labels(x$distance);

  app <- make.pagoda.app(tamr2,tam,varinfo,go.env,pwpca,clpca,cell.clustering=hc,title=name)#,embedding=NA)
  #show.app(app,name,browse=F,port=port)
  #saveRDS(app,file=paste(name,"app.rds",sep="."))
  #app <- readRDS(paste(name,"app.rds",sep="."));
  #show.app(app,name,port=1461,browse=F)


  return(invisible(app))
}


