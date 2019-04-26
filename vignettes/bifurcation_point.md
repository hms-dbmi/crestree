---
title: "Analysis of bifurcation points"
#author: "Ruslan Soldatov"
#date: "2019-04-25"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  %\VignetteIndexEntry{Analysis of branching trajectories}
---



This vignette describes analysis of individual bifurcation points based on a reconstructed transcriptional tree. It explores bifurcation point between sensory and autonomic nervous systems in neural crest. The guideline starts with tree reconstruction, identifies fate-specific genes and estimates timing of their activation, assess existence and formation of fate-biases, and predicts time of genes inclusion in fate-biased phase.

## Preliminaries: loading the libraries and neural crest data


```r
library(igraph)
library(mgcv)
library(quadprog) 
library(pcaMethods) 
library(Rcpp) 
library(inline) 
library(RcppArmadillo) 
#library(Rfast)
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
```
  
## Run tree reconstruction
The detailed guide for tree reconstruction is described in https://github.com/hms-dbmi/crestree/blob/master/vignettes/tree_guide.md. The tree can be either manually reconstructed:

```r
M <- length(nc.cells); 
lambda <- 250; 
sigma <- 0.04
ppt <- ppt.tree(X=wgm[,nc.cells], W=wgwm[,nc.cells], emb=emb, lambda=250, sigma=0.04, metrics="cosine", M=M,
                err.cut = 5e-3, n.steps=30, seed=1, plot=FALSE)
ppt <- cleanup.branches(ppt,tips.remove = c(139,295))
ppt <- setroot(ppt,355)
ppt <- project.cells.onto.ppt(ppt,emb,n.mapping = 100)
```

![plot of chunk unnamed-chunk-146](figure/unnamed-chunk-146-1.png)
or downloaded from:

```r
ppt <- readRDS(url("http://pklab.med.harvard.edu/ruslan/neural_crest/tree_structure_full.rds"))

#ppt <- readRDS("/d0/home/solrust/NC/resource_files/tree_structure.rds")
```

## Identification of branch-specific genes

Bifurcation point is charactarized by a progenitor and derivative branches. We thus start with selection a root of progenitor branch and two leaves of derivative branches:

```r
plotppt(ppt,emb,tips=TRUE,forks=FALSE,cex.tree = 0.2,lwd.tree = 2)
```

![plot of chunk unnamed-chunk-148](figure/unnamed-chunk-148-1.png)

```r
root <- 355
leaves <- c(165,91)
```

A routine `test.fork.genes` performs assessment of genes differentially expression between post-bifurcation branches:

```r
fork.de <- test.fork.genes(ppt,fpm[,],root=root,leaves=leaves,n.mapping = 10,n.cores=30)
```

A table `fork.de` contains summary statistics of fold change `effect`, p-value `p` and adjusted p-value `fdr`  of differential expression between branches, magnitude `pd1.a` and p-value `pd1.p` of expression changes from derivative branch 1 to progenitor branch:

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

We next consider a gene to be preferentially expressed along the first/second branch if it has `effect.b1`/`effect.b2` increased expression compared to another post-bifurcation branch and significant increase (p < 0.05) relative to progenitor branch:


Column `state` charactarizes genes that are specific to first (1), second (2), or neither (0) of derivative branches.


For consistency with original results, we also limit genes to `genes.tree` set associated with the tree:
 chr [1:98] "Rdh10" "Hes6" "Cxcr4" "Nfasc" "5730559C18Rik" "Rgs16" ...
 chr [1:122] "Serpine2" "Lrrfip1" "Cdh19" "Ralb" "Angptl1" "Pbx1" ...

## Classification of early and late modules

Dynamics of a gene expression is reflected in timing of activation and optimum expression. Routine `activation.fork` estimates timing of optimum expression of smoothed expression and activation point as first passage of derivative through `deriv.cutoff` cutoff:

estimate activation patterns .. branch 1
estimate activation patterns .. branch 2

`fork.de.act` table provides additional columns `optimum` and `activation` for genes predicted to be differentially expressed (`stat` = 1 or 2).

Branch-specific sets of genes (genes.sensory,genes.autonomic) can now be partitioned in early and late genes based on time of expression activation. A logical solution is to orient early/late genes relative to bifurcation point. Timing of root, bifurcation point and leaves are
       root bifurcation     leave 1     leave 2 
    0.00000    17.57189    27.59971    31.49521 

We use `cutoff = 16.0` on timing of activation to define early and late genes:


Now we can check if early/late genes modules follow co-activation or mutually-exclusive patterns:



```r
cells <- rownames(ppt$cell.summary)[ppt$cell.summary$seg %in% extract.subtree(ppt,c(root,leaves))$segs]
par(mfrow=c(1,2))
plot(t(programs[c(1,3),cells]),col=ppt$cell.summary[cells,]$color,pch=19,cex=0.5)
plot(t(programs[c(2,4),cells]),col=ppt$cell.summary[cells,]$color,pch=19,cex=0.5)
```

![plot of chunk unnamed-chunk-159](figure/unnamed-chunk-159-1.png)


## Coordination of fate biases

Co-activation of both fate-specific programs poses a question of when cell acquire bias in favor of one or another program. For that, we look for coordinated expression of each module inside more homogeneous subpopulations. First, bifurcation fork is partitioned in non-intersecting windows of `wind` cells:

Visualization of group of cells assigned to each non-intersecting window:

```r
fig_cells <- fig.cells(emb,freq)
marrangeGrob( c(fig_cells),ncol=length(fig_cells),nrow=1,top=NA)
```

![plot of chunk unnamed-chunk-161](figure/unnamed-chunk-161-1.png)

Windows can be also selected manually, below we follow selection used in the paper:



```r
freq <- slide.cells(ppt,root,leaves,wind=50,regions=regions,n.cores=10)
fig_cells <- fig.cells(emb,freq)
marrangeGrob( c(fig_cells),ncol=length(fig_cells),nrow=1,top=NA)
```

![plot of chunk unnamed-chunk-163](figure/unnamed-chunk-163-1.png)

Routine `slide.cors` next estimates average correlation of each early fate-specific gene with both modules (genes.sensory.early and genes.autonomic.early) in each window of cells:

```r
cors <- slide.cors(freq,fpm,genes.sensory.early,genes.autonomic.early)
```

Now joint visualization enables tracking how genes of fate-specific modules coordinate expression during progression along pseudotime:

```r
fig_cor <- fig.cors(cors,genes.sensory.early,genes.autonomic.early)
marrangeGrob( c(fig_cells,fig_cor),ncol=length(fig_cells),nrow=2,
              layout_matrix = matrix(seq_len(2*length(fig_cells)), nrow = 2, ncol = length(fig_cells),byrow=TRUE),top=NA)
```

![plot of chunk unnamed-chunk-165](figure/unnamed-chunk-165-1.png)
  
To obtain more contrasted (and reproducible with the paper) view, a set of early genes could be further cleaned up by removing fate-specific genes having low correlation with its modules around bifurcation point:

```r
corA <- cors[[5]][,1]
genesetA <- names(which(corA[genes.sensory.early] > 0.07))

corB <- cors[[5]][,2]
genesetB <- names(which(corB[genes.autonomic.early] > 0.07))
```

Re-calculate average window-specific correlations for cleaned up sets of genes genesetA and genesetB:

```r
cors <- slide.cors(freq,fpm,genesetA,genesetB)
fig_cor <- fig.cors(cors,genesetA,genesetB)
marrangeGrob( c(fig_cells,fig_cor),ncol=length(fig_cells),nrow=2,
              layout_matrix = matrix(seq_len(2*length(fig_cells)), nrow = 2, ncol = length(fig_cells),byrow=TRUE),top=NA)
```

![plot of chunk unnamed-chunk-167](figure/unnamed-chunk-167-1.png)


More generally, formal trends of local coordination of fate-specific modules along branching trajectories can be estimated:


```r
w=30
step=10
crd <- synchro(ppt,fpm,root,leaves,genesetA,genesetB,w,step,n.mapping=100,n.points = 300,span.smooth = 0.1,perm=FALSE)
```
And visualized:

```r
visualize.synchro(crd)
```

![plot of chunk unnamed-chunk-169](figure/unnamed-chunk-169-1.png)



