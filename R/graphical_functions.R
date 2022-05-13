#' The sparsity of a graph
#'
#' Finds the sparsity of a graph. The measure has a value in \eqn{[0,1]} where a value of \eqn{1} means there are edges between all nodes while a value of \eqn{0} means there are no edges in the graph.
#'
#' @param g The adjacency matrix of the graph.
#'
#'
#' @seealso \code{\link{jointGHS}} and \code{\link{confusion.matrix}}
#'
#'
#' @return Numeric sparsity.
#'
#'
#' @export
#'
#' @author Camilla Lingjaerde
#'
#'
#' @examples
#'
#' # example 1: simple example with two unrelated matrices
#' # generate data
#' set.seed(123)
#' n <- 80
#' p <- 100
#' dat <- huge::huge.generator(n = n, d = p, graph = "scale-free")
#' sparsity(dat$theta)
sparsity <- function(g) {
  if (nrow(g) != ncol(g)) stop("adjacency matrix must be square.")
  g <- g + 0
  p <- dim(g)[1]
  if (g[1, 1] == 0) { # If diagonal elements are not given, do not subtract p.
    return(sum(g != 0) / (p^2 - p))
  }
  else {
    return((sum(g != 0) - p) / (p^2 - p))
  }
}

#' Calculates the log likelihood of a precision matrix
#'
#' Finds the the profile log likelihood of a precision matrix given a sample covariance matrix, in the setting of Gaussian graphical models.
#'
#' @param sample.cov The sample covariance matrix of the observed data.
#' @param theta The estimated precision matrix.
#' @param n  The number of observations.
#'
#' @seealso \code{\link{jointGHS}}
#'
#' @return Numerical value of the multivariate Gaussian profile log likelihood of the precision matrix.
#'
#' @export
#'
#' @author Camilla Lingjaerde
#'
#'
#' @examples
#'
#' # example 1: simple example where we check the log likelihood of the true precision matrix
#' # generate data
#' set.seed(123)
#' n <- 80
#' p <- 100
#' dat <- huge::huge.generator(n = n, d = p, graph = "scale-free")
#' prec.mat <- dat$omega # true precision matrix
#' gaussianloglik(var(dat$data), prec.mat, n)
#' @keywords internal
gaussianloglik <- function(sample.cov, theta, n) {
  if (mean(dim(sample.cov) == dim(theta)) != 1) stop("matrices must have the same dimension")
  if (det(theta) <= 0) stop("precision matrix must be positive definite.")
  if (!isSymmetric(sample.cov)) stop("sample covariance matrix must be symmetric")
  if (n <= 0) stop("number of observations n must be positive")
  p <- nrow(theta)
  return(-p * n * log(2 * pi) / 2 + n * log(det(theta)) / 2 -
           n * sum(diag(sample.cov %*% theta)) / 2)
}

#' Calculates the AIC score of an estimated precision matrix
#'
#' Finds the AIC score of an estimated precision matrix given a sample covariance matrix, in the setting of Gaussian graphical models.
#'
#' @param sample.cov The sample covariance matrix of the observed data.
#' @param theta The estimated precision matrix.
#' @param n  The number of observations.
#'
#' @seealso \code{\link{jointGHS}}
#'
#' @return The AIC score.
#'
#' @export
#'
#' @author Camilla Lingjaerde
#'
#'
#' @examples
#'
#' # example 1: simple example where we check the AIC score of the true precision matrix
#' # generate data
#' set.seed(123)
#' n <- 80
#' p <- 100
#' dat <- huge::huge.generator(n = n, d = p, graph = "scale-free")
#' prec.mat <- dat$omega # true precision matrix
#' gaussianAIC(var(dat$data), prec.mat, n)
gaussianAIC <- function(sample.cov, theta, n) {
  p <- nrow(theta)
  theta2 <- theta
  diag(theta2) <- rep(0, p)
  d <- sum(theta2 != 0) / 2
  return(-2 * gaussianloglik(sample.cov, theta, n) + 2 * d)
}

#' Find the confusion matrix of a graph and an estimate of it
#'
#' Finds the confusion matrix of an inferred graph and the graph it aims to estimate.
#'
#' @param g The adjacency matrix of the true graph.
#'
#' @param g.hat The adjacency matrix of the estimated graph.
#'
#' @seealso \code{\link{jointGHS}}
#'
#'
#' @return The \eqn{2} by \eqn{2} confusion matrix. Element \code{[1,1]} is the number of true positives, element \code{[1,2]} the number of false positives, element  \code{[2,1]} the number of false negatives and  \code{[2,2]} the number of true negatives.
#'
#'
#' @export
#'
#' @author Camilla Lingjaerde
#'
#'
#' @examples
#'
#' # example 1: simple example with two unrelated matrices
#' # generate data and prior matrix
#' set.seed(123)
#' x1 <- matrix(rnorm(40 * 20), ncol = 20)
#' x2 <- matrix(rnorm(40 * 20), ncol = 20)
#' # the precision matrices
#' prec.mat.1 <- solve(var(x1))
#' prec.mat.2 <- solve(var(x2))
#' # make them sparse by threshold values. Do not make diagonal elements zero.
#' prec.mat.1[which(abs(prec.mat.1) < 1.5 & !diag(ncol(prec.mat.1)))] <- 0
#' prec.mat.2[which(abs(prec.mat.2) < 1.5 & !diag(ncol(prec.mat.2)))] <- 0
#' confusion.matrix(prec.mat.1 != 0, prec.mat.2 != 0)
confusion.matrix <- function(g, g.hat) {
  if (mean(dim(g[, ]) == dim(g.hat[, ])) != 1) stop("matrices must have the same dimension")
  if (mean((g[, ] + 0) %in% c(0, 1)) != 1 | mean((g.hat[, ] + 0) %in% c(0, 1)) != 1) stop("g and g.hat must be adjacency matrices with elements in {0,1}")
  p <- nrow(g[, ])
  g <- g[, ]
  g.hat <- g.hat[, ]
  diag(g) <- rep(0, p)
  diag(g.hat[, ]) <- rep(0, p) # Do not include diagonal elements.
  # Divide by 10 to avoid integer overflow.
  tp <- sum(g.hat[, ] == 1 & g[, ] == 1) / 10 # True positives.
  fp <- sum(g.hat[, ] == 1 & g[, ] == 0) / 10 # False positives
  tn <- (sum(g.hat[, ] == 0 & g[, ] == 0) - p) / 10 # True negatives.
  fn <- sum(g.hat[, ] == 0 & g[, ] == 1) / 10 # False negatives.
  return(matrix(10 * c(tp, fp, fn, tn), nrow = 2, byrow = T) / 2)
}

#' Find the recall of an estimated graph
#'
#' Computes the recall of an inferred graph in terms of predicting the edges of the true one. The measure has a value in \eqn{[0,1]} where a value of \eqn{1} means all edges in \code{g} are included in \code{g.hat} while a value of \eqn{0} means no true edges are included.
#'
#' @param g The adjacency matrix of the true graph.
#'
#' @param g.hat The adjacency matrix of the estimated graph.
#'
#' @seealso \code{\link{jointGHS}} and \code{\link{confusion.matrix}}
#'
#'
#' @return Numeric recall.
#'
#' @export
#'
#' @author Camilla Lingjaerde
#'
#'
#' @examples
#'
#' # example 1: simple example with two unrelated matrices
#' # generate data and prior matrix
#' set.seed(123)
#' x1 <- matrix(rnorm(40 * 20), ncol = 20)
#' x2 <- matrix(rnorm(40 * 20), ncol = 20)
#' # the precision matrices
#' prec.mat.1 <- solve(var(x1))
#' prec.mat.2 <- solve(var(x2))
#' # make them sparse by threshold values. Do not make diagonal elements zero.
#' prec.mat.1[which(abs(prec.mat.1) < 1.5 & !diag(ncol(prec.mat.1)))] <- 0
#' prec.mat.2[which(abs(prec.mat.2) < 1.5 & !diag(ncol(prec.mat.2)))] <- 0
#' recall(prec.mat.1 != 0, prec.mat.2 != 0)
recall <- function(g, g.hat) {
  # use confusion.matrix to find FP, TP, TN and FN.
  conf.mat <- confusion.matrix(g, g.hat)
  # check if there are no edges in g and hence no edges to predict. If so, we say the recall is one.
  if (conf.mat[1, 1] == 0 & conf.mat[2, 1] == 0) {
    return(1)
  } # else compute recall
  else {
    return(conf.mat[1, 1] / (conf.mat[1, 1] + conf.mat[2, 1]))
  }
}

#' Find the precision of an estimated graph
#'
#' Computes the precision of an inferred graph in terms of predicting the edges of the true one. The measure has a value in \eqn{[0,1]} where a value of \eqn{1} means perfect prediction while a value of \eqn{0} means completely wrong.
#'
#' @param g The adjacency matrix of the true graph.
#'
#' @param g.hat The adjacency matrix of the estimated graph.
#'
#' @seealso \code{\link{jointGHS}} and \code{\link{confusion.matrix}}
#'
#'
#' @return Numeric precision.
#'
#'
#' @export
#'
#' @author Camilla Lingjaerde
#'
#'
#' @examples
#'
#' # example 1: simple example with two unrelated matrices
#' # generate data and prior matrix
#' set.seed(123)
#' x1 <- matrix(rnorm(40 * 20), ncol = 20)
#' x2 <- matrix(rnorm(40 * 20), ncol = 20)
#' # the precision matrices
#' prec.mat.1 <- solve(var(x1))
#' prec.mat.2 <- solve(var(x2))
#' # make them sparse by threshold values. Do not make diagonal elements zero.
#' prec.mat.1[which(abs(prec.mat.1) < 1.5 & !diag(ncol(prec.mat.1)))] <- 0
#' prec.mat.2[which(abs(prec.mat.2) < 1.5 & !diag(ncol(prec.mat.2)))] <- 0
#' precision(prec.mat.1 != 0, prec.mat.2 != 0)
precision <- function(g, g.hat) {
  # use confusion.matrix to find FP, TP, TN and FN.
  conf.mat <- confusion.matrix(g, g.hat)
  # check if there are no edges in g and hence no edges to predict. If so, we say the precision is one.
  if (conf.mat[1, 1] == 0 & conf.mat[2, 1] == 0) {
    return(1)
  } # check if there are no edges in g.hat and hence no predictions made. If so, we say the precision is one.
  else if (conf.mat[1, 1] == 0 & conf.mat[1, 2] == 0) {
    return(1)
  } # else compute precision
  else {
    return(conf.mat[1, 1] / (conf.mat[1, 1] + conf.mat[1, 2]))
  }
}