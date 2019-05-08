# Analysis of bifurcation points


This vignette describes analysis of individual bifurcation points based on a reconstructed transcriptional tree. As example, it explores bifurcation point between sensory and autonomic nervous systems in neural crest. The guideline starts with tree reconstruction, identifies fate-specific genes and estimates timing of their activation, assess existence and formation of fate-biases, and predicts time of genes inclusion in fate-biased phase.

## Preliminaries: loading the libraries and neural crest data


```r
library(igraph)
library(mgcv)
library(quadprog) 
library(pcaMethods) 
library(Rcpp) 
library(inline) 
library(RcppArmadillo) 
library(crestree)
library(ggplot2); library(gridExtra); library(grid);

data(crest)
emb <- crest$emb
clcol <- crest$clcol
nc.cells <- crest$nc.cells
wgm <- crest$wgm
wgwm <- crest$wgwm # matrix of expression weights
fpm <- read.table("http://pklab.med.harvard.edu/ruslan/neural_crest/fpm.txt",header=TRUE)
fpm <- as.matrix(fpm)
genes.tree <- crest$genes.tree
```
  
## Run tree reconstruction
The detailed guide for tree reconstruction is described in https://github.com/hms-dbmi/crestree/blob/master/vignettes/tree_guide.md. The tree can be either manually reconstructed:

```r
M <- length(nc.cells); 
lambda <- 250; 
sigma <- 0.04
#ppt <- ppt.tree(X=wgm[,nc.cells], W=wgwm[,nc.cells], emb=emb, lambda=250, sigma=0.04, metrics="cosine", M=M, err.cut = 5e-3, n.steps=30, seed=1, plot=FALSE)
#ppt <- cleanup.branches(ppt,tips.remove = c(139,295))
#ppt <- setroot(ppt,355)
#ppt <- project.cells.onto.ppt(ppt,n.mapping = 100)
```

Alternatively, the tree object used in the paper can be downloaded from:

```r
ppt <- readRDS(url("http://pklab.med.harvard.edu/ruslan/neural_crest/tree_structure_full.rds"))
```

## Identification of branch-specific genes

Bifurcation point is charactarized by a progenitor and derivative branches.

```r
plotppt(ppt,emb,tips=TRUE,forks=FALSE,cex.tree = 0.2,lwd.tree = 2)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

We thus start with selection a root of progenitor branch (355) and two leaves of derivative branches (165 and 91):

```r
root <- 355
leaves <- c(165,91)
```

Here is a fork charactarizing sensory-autonomic bifurcation point:

```r
subtree <- extract.subtree(ppt,c(root,leaves))
plotppt(ppt,emb,tips=TRUE,forks=FALSE,cex.tree = 0.3,lwd.tree = 3,subtree=subtree)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

A routine `test.fork.genes` performs assessment of genes differentially expression between post-bifurcation branches:

```r
fork.de <- test.fork.genes(ppt,fpm,root=root,leaves=leaves,n.mapping = 10,n.cores=30)
```

A table `fork.de` contains summary statistics of fold change `effect`, p-value `p` and adjusted p-value `fdr`  of differential expression between branches, magnitude `pd1.a` (`pd2.a`) and p-value `pd1.p` (`pd2.p`) of expression changes from first (second) derivative branch to progenitor branch:

```r
head(fork.de[order(fork.de$p),],)
```


|        |    effect|  p| fdr|  st| stf|     pd1.a| pd1.p|      pd2.a|     pd2.p|
|:-------|---------:|--:|---:|---:|---:|---------:|-----:|----------:|---------:|
|Neurog2 | 1.9032796|  0|   0| 0.9| 0.9| 0.1012255|     0| -0.0262510| 0.0001858|
|Dll1    | 1.4404789|  0|   0| 1.0| 1.0| 0.0976255|     0|  0.0043581| 0.3941055|
|Srrm4   | 0.8383684|  0|   0| 1.0| 1.0| 0.0646809|     0|  0.0099734| 0.0002647|
|Mfng    | 0.9532576|  0|   0| 1.0| 1.0| 0.0630761|     0|  0.0034360| 0.4424417|
|Eya2    | 1.0311280|  0|   0| 0.9| 0.9| 0.0737302|     0|  0.0082446| 0.0035890|
|Pcdh8   | 1.1417276|  0|   0| 1.0| 0.9| 0.0515220|     0| -0.0193619| 0.0000199|

We next consider a gene to be preferentially expressed along the first (second) branch if it has `effect.b1` (`effect.b2`) increased expression compared to another post-bifurcation branch and significant increase (p < 0.05) relative to progenitor branch:


```r
fork.de <- branch.specific.genes(fork.de,effect.b1 = 0.1,effect.b2 = 0.3)
```

Column `state` characterizes genes that are specific to first (1), second (2), or neither (0) of derivative branches.

```r
genes.sensory <- rownames(fork.de)[fork.de$state==1]
genes.autonomic  <- rownames(fork.de)[fork.de$state==2]
```

For consistency with the original results, we also limit genes to `genes.tree` set associated with the tree:

```r
genes.sensory <- intersect(genes.sensory,genes.tree)
str(genes.sensory)
##  chr [1:98] "Rdh10" "Hes6" "Cxcr4" "Nfasc" "5730559C18Rik" "Rgs16" ...

genes.autonomic <- intersect(genes.autonomic,genes.tree)
str(genes.autonomic)
##  chr [1:122] "Serpine2" "Lrrfip1" "Cdh19" "Ralb" "Angptl1" "Pbx1" ...
```

## Classification of early and late modules

Dynamics of gene expression is reflected in timing of activation and expression optimum. Routine `activation.fork` estimates timing of optimum expression of smoothed expression and activation point as a first passage of derivative through `deriv.cutoff` cutoff:


```r
fork.de.act <- activation.fork(ppt,fork.de,fpm,root,leaves,deriv.cutoff = 0.015,n.mapping=10,n.cores=10)
```

`fork.de.act` table provides additional columns `optimum` and `activation` for genes predicted to be differentially expressed (`stat` = 1 or 2).

Branch-specific sets of genes (`genes.sensory` and `genes.autonomic`) can now be partitioned in early and late genes based on time of expression activation. A logical solution is to orient early/late genes relative to bifurcation point. Timing of root, bifurcation point and leaves are


```r
fork.pt(ppt,root,leaves)
```


| root| bifurcation|  leave 1|  leave 2|
|----:|-----------:|--------:|--------:|
|    0|     17.5719| 27.59971| 31.49521|

We use `cutoff=16.0` on timing of activation to define early and late genes:


```r
cutoff <- 16.0
```


```r
genes.sensory.late <- genes.sensory[fork.de.act[genes.sensory,]$activation > cutoff]
genes.sensory.early <- setdiff(genes.sensory,genes.sensory.late)

genes.autonomic.late <- genes.autonomic[fork.de.act[genes.autonomic,]$activation > cutoff]
genes.autonomic.early <- setdiff(genes.autonomic,genes.autonomic.late)
```

Now we can check if early/late genes modules follow co-activation or mutually-exclusive patterns:



```r
cells <- rownames(ppt$cell.summary)[ppt$cell.summary$seg %in% extract.subtree(ppt,c(root,leaves))$segs]
par(mfrow=c(1,2))
plot(t(programs[c(1,3),cells]),col=ppt$cell.summary[cells,]$color,pch=19,cex=0.5)
plot(t(programs[c(2,4),cells]),col=ppt$cell.summary[cells,]$color,pch=19,cex=0.5)
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png)


## Coordination of fate biases

Co-activation of both fate-specific programs poses a question of when cell acquire bias in favor of one or another program. For that, we look for coordinated expression of each module inside more homogeneous subpopulations. First, bifurcation fork is partitioned in non-intersecting windows of `wind` cells:

```r
freq <- slide.cells(ppt,root,leaves,wind=50,n.cores=10)
```
Visualization of group of cells assigned to each non-intersecting window:

```r
fig_cells <- fig.cells(emb,freq)
marrangeGrob( c(fig_cells),ncol=length(fig_cells),nrow=1,top=NA)
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21-1.png)

Windows can be also selected manually, below we follow selection used in the paper:

```r
regions = list( list(7,151,200,1),list(7,101,151,1),list(7,51,100,1),list(7,1,50,1),list(list(6,5,1,2),1,50, -1),list(list(6,5,1,2),51,100, -1),list(5,1,50,1),list(1,1,50,1))
```


```r
freq <- slide.cells(ppt,root,leaves,wind=50,regions=regions,n.cores=10)
fig_cells <- fig.cells(emb,freq)
marrangeGrob( c(fig_cells),ncol=length(fig_cells),nrow=1,top=NA)
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23-1.png)

Next, routine `slide.cors` estimates average correlation of each early fate-specific gene with both modules (`genes.sensory.early` and `genes.autonomic.early`) in each window of cells:

```r
cors <- slide.cors(freq,fpm,genes.sensory.early,genes.autonomic.early)
```

Now joint visualization enables tracking how genes of fate-specific modules coordinate expression during progression along pseudotime:

```r
fig_cor <- fig.cors(cors,genes.sensory.early,genes.autonomic.early)
marrangeGrob( c(fig_cells,fig_cor),ncol=length(fig_cells),nrow=2,
              layout_matrix = matrix(seq_len(2*length(fig_cells)), nrow = 2, ncol = length(fig_cells),byrow=TRUE),top=NA)
```

![plot of chunk unnamed-chunk-25](figure/unnamed-chunk-25-1.png)
  
To obtain more contrasted (and reproducible with the paper) view, a set of early genes could be further cleaned up by removing fate-specific genes having low correlation with its modules around bifurcation point:

```r
corA <- cors[[5]][,1]
genesetA <- names(which(corA[genes.sensory.early] > 0.07))

corB <- cors[[5]][,2]
genesetB <- names(which(corB[genes.autonomic.early] > 0.07))
```

Re-estimation average window-specific correlations for cleaned up sets of genes `genesetA` and `genesetB`:

```r
cors <- slide.cors(freq,fpm,genesetA,genesetB)
fig_cor <- fig.cors(cors,genesetA,genesetB)
marrangeGrob( c(fig_cells,fig_cor),ncol=length(fig_cells),nrow=2,
              layout_matrix = matrix(seq_len(2*length(fig_cells)), nrow = 2, ncol = length(fig_cells),byrow=TRUE),top=NA)
```

![plot of chunk unnamed-chunk-27](figure_bifurcation/unnamed-chunk-27-1.png)


More generally, formal trends of local coordination of fate-specific modules along branching trajectories can be estimated using `synchro` routine:


```r
w=30
step=10
crd <- synchro(ppt,fpm,root,leaves,genesetA,genesetB,w,step,n.mapping=100,n.points = 300,span.smooth = 0.1,perm=FALSE)
```
And visualized:

```r
visualize.synchro(crd)
```

![plot of chunk unnamed-chunk-29](figure_bifurcation/unnamed-chunk-29-1.png)



## Formation of correlated module

In order to infer genes orchestrating emergence of local intra-module correlations, we developed an iterative approach that estimates timing of genes inclusion in a core of correlated genes of the module. Routine `onset` estimates inclusion timing of early sensory genes `genes.sensory.early` for vector `mappings` of probabilistic cell projections:


```r
inclusion.sensory <- onset(ppt,geneset=genes.sensory.early,nodes=c(root,leaves),expr=fpm,alp=20,w=40,step=10,n.cells=280,mappings=1:100,do.perm=FALSE,winperm=30,n.cores=30,permut.n=10,n.cores1=1)
```

Average inclusion timing of each gene is

```r
head(apply(inclusion.sensory,1,mean))
##         Rdh10          Hes6         Cxcr4 5730559C18Rik          Utrn 
##      11.17453      13.00036      15.66524      12.09579      13.38760 
##          Gamt 
##      15.45684
```


As an example, predicted cumulative probability of gene *Pou4f1* inclusion is shown below:

```r
show.gene.inclusion("Pou4f1",inclusion.sensory)
```

![plot of chunk unnamed-chunk-32](figure_bifurcation/unnamed-chunk-32-1.png)

Overall summary statistics of inclusion probabilistics is visualized using `show.inclusion.summary` routine with each row of the heatmap reflecting cumulative inclusion probability of a single gene (as in example above):

```r
show.inclusion.summary(inclusion.sensory,gns.show=c("Neurog2","Pou4f1"))
```

![plot of chunk unnamed-chunk-33](figure_bifurcation/unnamed-chunk-33-1.png)

Analogous inference of inclusion timing for early atuonomic module `genes.autonomic.early`:

```r
inclusion.autonomic <- onset(ppt,geneset=genes.autonomic.early,nodes=c(root,leaves),expr=fpm,alp=20,w=40,step=10,n.cells=280,mappings=1:100,do.perm=FALSE,winperm=30,n.cores=30,permut.n=10,n.cores1=1)

show.inclusion.summary(inclusion.autonomic,gns.show = c("Mef2c","Pbx1"))
```

![plot of chunk unnamed-chunk-34](figure_bifurcation/unnamed-chunk-34-1.png)

## Finding optimal parameter to regulate false positive

For comparison, inclusion events are not detected in a control expression matrix, where expression levels were locally permuted using option `do.perm=TRUE` among `winperm=10` batches of cells along pseudotime:


```r
inclusion.sensory.control <- onset(ppt,geneset=genes.sensory.early,nodes=c(root,leaves),expr=fpm,alp=20,w=40,step=10,n.cells=280,
                 mappings=1:10,do.perm=TRUE,winperm=10,n.cores=30,permut.n=10,n.cores1=1)
show.inclusion.summary(inclusion.sensory.control,gns.show=NA)
```

![plot of chunk unnamed-chunk-35](figure_bifurcation/unnamed-chunk-35-1.png)


Parameter `alp` in routine `onset`, which was set to 20 in the analysis above, regulates stringency of inclusion point detection. Assessing summary statistics on inclusion times for a range of `alp` levels for real and permutated expression profiles enables selection of optimal `alp`. First, select a range of `alp` levels:


```r
alps <- seq(1,40,length.out = 10)
```

Then estimate mean gene inclusion time for `alps` levels:

```r
incl.real <- unlist(lapply(alps,function(alp){
  cat(alp);cat('\n')
  res.est <- onset(ppt,geneset=genes.sensory.early,nodes=c(root,leaves),expr=fpm,alp=alp,w=40,step=10,n.cells=280, mappings=1:10,do.perm=FALSE,winperm=10,n.cores=30,permut.n=10,n.cores1=1)
  return(mean(res.est))
}))
```

Output is average inclusion time for each alp level:


```r
head(incl.real)
## [1]  8.349644 10.047820 11.183495 12.011598 12.615876 12.951238
```

Estimate inclusion time for locally permuted expression levels by setting `do.perm=TRUE` and local permutation window `winperm=10`:


```r
incl.perm <- unlist(lapply( alps,function(alp){
  cat(alp);cat('\n')
  res.est <- onset(ppt,geneset=genes.sensory.early,nodes=c(root,leaves),expr=fpm,alp=alp,w=40,step=10,n.cells=280, mappings=1:10,do.perm=TRUE,winperm=30,n.cores=30,permut.n=10,n.cores1=1)
  mean(res.est)
}))
```

Comparison of early inclusion timing among all genes for real and control expression levels:


```r
inclusion.stat(alps, incl.real, incl.perm)
```

![plot of chunk unnamed-chunk-40](figure_bifurcation/unnamed-chunk-40-1.png)

From that one can see that if `alp > 10` than overall false positive rate in permutations is close to zero, while detection sensitivity of real data gradually declines.

