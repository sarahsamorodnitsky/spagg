% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/means.R
\name{landau.avg}
\alias{landau.avg}
\title{Compute aggregation of Ripley's K based on Landau et al. (2004)}
\usage{
landau.avg(K.vec, area.vec, n.vec)
}
\arguments{
\item{K.vec}{Vector of Ripley's K estimates across a series of images}

\item{area.vec}{Vector of sizes of the window for each image}

\item{n.vec}{Vector of the number of points in each image}
}
\value{
Returns a weighted average of the Ripley's K values based on the sizes of the
images and numbers of points
}
\description{
Compute aggregation of Ripley's K based on Landau et al. (2004)
}
\examples{
# Generate simulated Ripley's K values, areas, and number of points
K.vec <- rnorm(3, mean = 1000, sd = 100)
area.vec <- rnorm(3, mean = 20000, sd = 100)
n.vec <- rnorm(3, mean = 300, sd = 50)

# Compute a weighted average
landau.avg(K.vec, area.vec, n.vec)
}
\references{
Landau et al. (2004) Nonparametric One-way Analysis of Variance of Replicated Bivariate Spatial Point Patterns.
Biometrical Journal, 46(1), 19-34.
}
