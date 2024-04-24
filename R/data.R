#' Subset of non-small cell lung cancer dataset
#'
#' This is the cleaned version of a multiplexed immunohistochemistry (mIHC) dataset
#' obtained from the Johnson et al. (2021) study of non-small cell lung cancer.
#'
#' @format ## `lung_df_tumor`
#' A data frame with 871,171 rows and 9 columns:
#' \describe{
#'   \item{id}{Image ID corresponding to the ROI}
#'   \item{PID}{Patient ID corresponding to the sample (multiple ROIs per sample)}
#'   \item{cell_id}{Label for the cells in each ROI}
#'   \item{x}{The x-coordinate for each cell}
#'   \item{y}{The y-coordinate for each cell}
#'   \item{mhcII_status, mhcII_high}{Binary indicators of major histocompatibility complex II (MHCII).
#'   One column provides the labels as strings and the other as numeric values.}
#'   \item{tissue_category}{Tissue where cells were detected. Cells were subset to just those in the tumor.}
#'   \item{type}{Cell type}
#' }
#' @source <http://juliawrobel.com/MI_tutorial/MI_Data.html>
#' @source <https://www.sciencedirect.com/science/article/pii/S1556086421021754>
"lung_df_tumor"
