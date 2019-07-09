source("pagoda_call/pagoda_call.functions.R")

### download count matrix
cd <- read.table("http://pklab.med.harvard.edu/ruslan/neural_crest/cd_p85-p87.txt",row.names = 1,header=TRUE)
cdb <- gsub("_.*","",colnames(cd)); names(cdb) <- colnames(cd);
cn <- "p85-87"

### Prepapre GO data
# Use the org.Mm.eg.db package for GO annotations
library(org.Mm.eg.db)
# Translate gene names to ids
ids <- unlist(lapply(mget(rownames(cd), org.Mm.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
# Reverse map
rids <- names(ids)
names(rids) <- ids
# Convert ids per GO category to gene names
go.env <- eapply(org.Mm.egGO2ALLEGS , function(x) as.character(na.omit(rids[x])))
go.env <- list2env(go.env)  # Convert to an environment
# Test
class(go.env)
head(ls(go.env)) # Look at gene set names

### set up parameters for PAGODA
n.cores <- 25; # number of parallel cores to use
min.cell.genes <- 3e3;
min.cell.reads <- 1e5;
min.gene.reads <- 10;
min.gene.cells <- 5;
min.nonfailed <- 8;
n.groups <- 10;
trim <- 3;

## run PAGODA processing script
res <- process.dataset(cd,"p85-87",batch=cdb)
## generate PAGODA app for object
app <- postprocess(res,"p85-87")

# estimate and subtract M-phase
m.x <- pagoda.show.pathways(unique(c("GO:0051301","GO:0007076","GO:0007049","GO:0005819")),res$varinfo,goenv=go.env,cell.clustering=app$results$hvc,margins=c(0,5),show.cell.dendrogram=T,box=T,showRowLabels=T,showPC=T,n.genes=20,cexRow=1.5,return.details=T)
varinfo.cc <- pagoda.subtract.aspect(res$varinfo,m.x$scores[,1])
# estimate and subtract S-phase
pat <- pagoda.show.pathways(c("GO:0006260","GO:0006270"),varinfo.cc,goenv=go.env,cell.clustering=app$results$hvc,margins=c(0,5),plot=T,show.cell.dendrogram=T,box=T,showRowLabels=T,showPC=T,n.genes=20,cexRow=1.5,return.details=F)
varinfo.cc <- pagoda.subtract.aspect(varinfo.cc,pat)
# re-estimate PAGODA object
res.cc <- process.dataset(dat=NULL,name=cn,cd=res$cd,batch=res$batch,varinfo=varinfo.cc,knn=res$knn)
app.cc <- postprocess(m.cdl[[cn]],cn,n.clusters=15,top.aspects=20,distance.threshold=0.8)

# The results can slightly differ from those in the paper due stochastisity of the method and updated versions of internally used packages
## final PAGODA object res.ss also can be downloaded from:
#res.load <- readRDS("http://pklab.med.harvard.edu/ruslan/neural_crest/p85-87.rds")
## PAGODA app can be downloaded from
#app.load <- readRDS("http://pklab.med.harvard.edu/ruslan/neural_crest/pagoda_rds/pagoda1_Wnt1_trunk_E95.rds")


tam <- pagoda.top.aspects(res.cc$pwpca[!grepl("custom",names(res.cc$pwpca))],res.cc$clpca,use.oe.scale=F,z.score=3,adjust.scores=T)
varinfo <- res.cc$varinfo
# get expression matrices on which to calculate the principal tree
# wgm - a smaller set of genes
gw <- tam$gw
gw <- gw[(rowSums(varinfo$matw)*varinfo$arv)[names(gw)] > 1]
mi <- match(names(gw), rownames(varinfo$mat))
#mi <- varinfo$arv>1.1
wgm <- varinfo$mat[mi, ]
wgwm <- varinfo$matw[mi, ]

# wgm2 - a largrer set of genes that includes pretty much all the expressed genes, except for those that were
#        penalized as part of batch-bias correction
mi <- varinfo$arv>0.5
wgm2 <- varinfo$mat[mi, ]
wgwm2 <- varinfo$matw[mi, ]


