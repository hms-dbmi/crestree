% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crestree.functions.R
\name{fit.associated.genes}
\alias{fit.associated.genes}
\title{Model gene expression levels as a function of tree positions.}
\usage{
fit.associated.genes(r, X, n.map = 1,
  n.cores = parallel::detectCores()/2, method = "ts", knn = 1,
  gamma = 1.5)
}
\arguments{
\item{r}{pptree object}

\item{X}{expressinon matrix of genes (rows) vs cells (columns)}

\item{n.map}{number of probabilistic cell-to-tree mappings to use}

\item{method}{method of modeling. Currently only splines with option 'ts' are supported.}

\item{knn}{use expression averaging among knn cells}

\item{gamma}{stringency of penalty.}
}
\value{
modified pptree object with new fields r$fit.list, r$fit.summary and r$fit.pattern. r$fit.pattern contains matrix of fitted gene expression levels
}
\description{
Model gene expression levels as a function of tree positions.
}
