
#' Batch Bayesian optimization using Thompson sampling
#'
#' @param model an object of class `hetGP`; e.g., as returned by `mleHetGP`
#' @param npoints an integer representing the desired number of samples
#' @param Xgrid a matrix of locations at which the samples are drawn
#'
#' @return a matrix containing the `npoints` best locations where next batch of simulations should be run
#' @export
#'
#' @examples
TS_npoints <- function(model, npoints, Xgrid){

  ## predict and conditional sample
  pred <- stats::predict(model, Xgrid, xprime = Xgrid)

  inflate.factor <- 2
  tTS <- mvtnorm::rmvnorm(npoints * 2, pred$mean, 1/2 * (pred$cov + t(pred$cov)))

  best_ids <- unique(apply(tTS, 1, which.min))
  if(length(best_ids) > npoints) best_ids <- best_ids[1:npoints]

  return(Xgrid[best_ids, ])
}
