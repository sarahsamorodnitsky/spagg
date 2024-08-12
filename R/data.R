#' Simulated point pattern dataset for examples
#'
#' This is a randomly-generated dataset for the purpose
#' of illustrating the usage of the functions contained in
#' this package. See spagg/data-raw/data.R for the code on
#' how to simulate this dataset.
#'
#' @format ## `data`
#' A data frame with 38515 rows and 6 columns
#' \describe{
#'  \item{PID}{Sample ID that groups regions-of-interest (ROI) together}
#'  \item{id}{ROI ID. This corresponds to an image within a larger sample}
#'  \item{cell.id}{Enumerates each cell within each ROI. Each row corresponds to a cell.}
#'  \item{x}{x-coordinate for each cell within each ROI }
#'  \item{y}{y-coordinate for each cell within each ROI}
#'  \item{type}{Cell-type label. Options are 'a' and 'b'. These are randomly generated.}
#'  \item{out}{Sample-level binary outcome. Note that this should be the same within a PID.}
#' }
"data"