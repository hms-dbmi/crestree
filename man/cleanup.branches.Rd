% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crestree.functions.R
\name{cleanup.branches}
\alias{cleanup.branches}
\title{Remove spurious branches of pptree}
\usage{
cleanup.branches(r, tips.remove = NULL, min.branch.length = 3)
}
\arguments{
\item{r}{ppt.tree result}

\item{tips.remove}{vector of tips indices to remove}

\item{min.branch.length}{remove all branches with length less or equal than min.branch.length principal points}

\item{tips.number}{select and retain only fixed number of tips (tips.number) that explain the most cell-cell variation.}
}
\value{
modified ppt.tree object with cleaned up structure
}
\description{
Remove spurious branches of pptree
}
