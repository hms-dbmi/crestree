#'  Processed neural crest single-cell data
#'
#' A dataset containing normalized expression matrices, t-SNE embedding and clusters.
#'
#' @docType data
#'
#' @usage data(crest)
#'
#' @format A list of 6 items.
#' \describe{
#'   \item{fpm}{normalized expression matrix}
#'   \item{wgm}{matrix of expression levels adjusted for mean-variance trend}
#'   \item{wgwm}{matrix of expression weights}
#'   \item{emb}{2D t-SNE embedding}
#'   \item{clcol}{cell colors of clusters}
#' \item{nc.cells}{a vector of neural crest cells. Other cells are neural tube cells.}
#' \item{motmat}{a matrix of target-TF scores}
#' }
#'
#' @keywords datasets
#'
#'
#' @source \url{web_link}
"crest"
