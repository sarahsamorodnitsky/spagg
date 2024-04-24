#' Compute aggregation of Ripley's K based on Landau et al. (2004)
#'
#' @param K.vec Vector of Ripley's K estimates across a series of images
#' @param area.vec Vector of sizes of the window for each image
#' @param n.vec Vector of the number of points in each image
#'
#' @return Returns a weighted average of the Ripley's K values based on the sizes of the
#' images and numbers of points
#' @references
#' Landau et al. (2004) Nonparametric One-way Analysis of Variance of Replicated Bivariate Spatial Point Patterns.
#' Biometrical Journal, 46(1), 19-34.
#' @export
#'
#' @examples
#' # Generate simulated Ripley's K values, areas, and number of points
#' K.vec <- rnorm(3, mean = 1000, sd = 100)
#' area.vec <- rnorm(3, mean = 20000, sd = 100)
#' n.vec <- rnorm(3, mean = 300, sd = 50)
#'
#' # Compute a weighted average
#' landau.avg(K.vec, area.vec, n.vec)
landau.avg <- function(K.vec, area.vec, n.vec) {

  # Calculate total number of points
  n <- sum(n.vec)

  # Calculate weighted average
  (sum(area.vec)/n^2) * sum((K.vec * n.vec^2)/area.vec)
}

#' Compute aggregation of Ripley's K based on Diggle et al. (1991)
#'
#' @param K.vec Vector of Ripley's K estimates across a series of images
#' @param n.vec Vector of the number of points in each image
#'
#' @return Returns a weighted average of the Ripley's K values based on the number of points in each image
#' @references Diggle et al. (1991). Analysis of Variance for Replicated Spatial Point Patterns in Clinical Neuroanatomy.
#' Journal of the American Statistical Association, 86(415), 618-625.
#' @export
#'
#' @examples
#' # Generate simulated Ripley's K values, areas, and number of points
#' K.vec <- rnorm(3, mean = 1000, sd = 100)
#' n.vec <- rnorm(3, mean = 300, sd = 50)
#'
#' # Compute a weighted average
#' diggle.avg(K.vec, n.vec)
diggle.avg <- function(K.vec, n.vec) {

  # Calculate the total number of points
  n <- sum(n.vec)

  # Calculate the weighted average
  (1/n) * sum(K.vec * n.vec)
}

#' Compute aggregation of Ripley's K based on Baddeley et al. (1993)
#'
#' @param K.vec Vector of Ripley's K estimates across a series of images
#' @param n.vec Vector of the number of points in each image
#'
#' @return Returns a weighted average of the Ripley's K values based on the number of points in each image
#' @references Baddeley et al. (1993). Analysis of a Three-Dimensional Point Pattern with Replication.
#' Journal of the Royal Statistical Society. Series C (Applied Statistics), 42(4), 641-668.
#' @export
#'
#' @examples
#' # Generate simulated Ripley's K values, areas, and number of points
#' K.vec <- rnorm(3, mean = 1000, sd = 100)
#' n.vec <- rnorm(3, mean = 300, sd = 50)
#'
#' # Compute a weighted average
#' baddeley.avg(K.vec, n.vec)
baddeley.avg <- function(K.vec, n.vec) {

  # Calculate the total number of points squared
  n2 <- sum(n.vec^2)

  # Calculate the weighted average
  (1/n2) * sum(K.vec * n.vec^2)
}
