#' Ensemble aggregation of spatial summaries (e.g. Ripley's K) using random weighting
#'
#' @param data data.frame or tibble containing data, including patient IDs,
#' image IDs, outcomes, and spatial summaries.
#' @param group Grouping variable for images, e.g. patient IDs.
#' @param outcome Column name for outcome in data.
#' @param cens Column name for event indicator column. If not using a survival outcome, leave as NULL.
#' @param model Model type. Options are "logistic" for logistic regression and
#' "survival" for a Cox proportional hazards model.
#' @param adjustments Column names for columns in data to adjust for in each ensemble replication
#' @param n.ensemble How many ensemble replications should be done? Default is 1000.
#' @param seed A seed for reproducible results. If left NULL, the funciton will generate a seed for you and return it.
#'
#' @return A list of length 2 containing the overall p-value obtained using the Cauchy combination
#' test (`pval`) and the seed used for reproducibility (`seed`). If the user provided a seed, then
#' `seed` will match what was given. Otherwise, the function randomly selects a seed and this is returned.
#' @references Liu et al. (2023) Ensemble methods for testing a global null. Journal of the Royal Statistical Society Series B: Statistical Methodology. 00. 1-26.
#' @export
#' @import spatstat spatstat.explore spatstat.geom dplyr survival ACAT tidyselect
#' @importFrom magrittr %>%
#' @examples
#' # Pick a radius to evaluate Ripley's K
#' r <- 30
#'
#' # Save the image IDs
#' ids <- unique(simdata$id)
#'
#' # Compute Ripley's K for tumor cells at r
#' K.vec <- c()
#' for (i in 1:length(ids)) {
#'  # Save the ith image
#'  image.i <- simdata %>%
#'    dplyr::filter(id == ids[i]) %>%
#'    dplyr::select(x,y,type)
#'
#'  # Convert to a point process object
#'  w <- spatstat.geom::convexhull.xy(image.i$x, image.i$y)
#'  image.i.subset <- image.i %>% dplyr::filter(type == "a")
#'  image.ppp <- spatstat.geom::as.ppp(image.i.subset, W = w)
#'
#'  # Compute Kest
#'  Ki <- spatstat.explore::Kest(image.ppp, r = 0:30)
#'
#'  # Save the result
#'  K.vec[i] <- Ki$iso[31]
#' }
#'
#' data.subset <- simdata %>%
#'  dplyr::select(tidyselect::all_of(c("id", "PID", "out"))) %>%
#'  dplyr::distinct()
#'
#' data.subset$spatial <- K.vec
#'
#' # Remove NaNs
#' data.subset <- data.subset %>% dplyr::filter(!is.na(spatial))
#'
#' # Test
#'ensemble.avg(data = data.subset, group = "PID", outcome = "out", model = "logistic")
ensemble.avg <- function(data, group, outcome, cens = NULL, model = "survival", adjustments = NULL, n.ensemble = 1000, seed = NULL) {

  # Check for a seed
  if (is.null(seed)) {
    seed <- sample(1:.Machine$integer.max, 1)
  }

  set.seed(seed)

  # Create a p-value based on many weighted averages
  p.values <- sapply(1:n.ensemble, function(ii) {

    # Generate weights
    weights <- stats::rnorm(nrow(data))
    scaled.weights <- abs(weights)/sum(abs(weights))

    # Add as a column
    data$weights <- scaled.weights

    # Calculate the weighted average
    data <- data %>%
      dplyr::group_by(get(group)) %>%
      dplyr::mutate(weight.avg = sum(spatial*weights)/sum(weights)) %>%
      dplyr::ungroup() %>%
      dplyr::select(tidyselect::all_of(c(group, outcome, cens, adjustments, "weight.avg"))) %>%
      dplyr::distinct()

    # Fit outcome model
    if (model == "survival") {

      # Fit Cox model
      if (is.null(adjustments)) {
        res <- survival::coxph(survival::Surv(get(outcome), get(cens)) ~ weight.avg, data = data)
      }

      if (!is.null(adjustments)) {
        res <- survival::coxph(survival::Surv(get(outcome), get(cens)) ~ weight.avg + get(adjustments), data = data)
      }

      # Save the p-value
      p.value <- summary(res)$coefficients[1,5]

    }

    if (model == "logistic") {

      # Fit logistic model
      if (is.null(adjustments)) {
        res <- stats::glm(get(outcome) ~ weight.avg, data = data, family = stats::binomial())
      }

      if (!is.null(adjustments)) {
        res <- stats::glm(get(outcome) ~ weight.avg + get(adjustments), data = data, family = stats::binomial())
      }

      # Save the p-value
      p.value <- summary(res)$coefficients[2,4]

    }

    # Return the p-value
    p.value

  })

  # Combine p-values
  pval <- ACAT::ACAT(p.values)

  # Return
  list(pval = pval, seed = seed)

}

#' Ensemble aggregation of spatial summaries (e.g. Ripley's K) using a combination
#' of random weighting and weighting by the number of points
#'
#' @param data data.frame or tibble containing data, including patient IDs,
#' image IDs, outcomes, and spatial summaries.
#' @param group Grouping variable for images, e.g. patient IDs.
#' @param outcome Column name for outcome in data.
#' @param cens Column name for event indicator column. If not using a survival outcome, leave as NULL.
#' @param model Model type. Options are "logistic" for logistic regression and
#' "survival" for a Cox proportional hazards model.
#' @param adjustments Column names for columns in data to adjust for in each ensemble replication
#' @param n.ensemble How many ensemble replications should be done? Default is 1000.
#' @param seed A seed for reproducible results. If left NULL, the funciton will generate a seed for you and return it.
#'
#' @return A list of length 2 containing the overall p-value obtained using the Cauchy combination
#' test (`pval`) and the seed used for reproducibility (`seed`). If the user provided a seed, then
#' `seed` will match what was given. Otherwise, the function randomly selects a seed and this is returned.
#' @export
#' @import spatstat spatstat.explore spatstat.geom dplyr survival ACAT tidyselect
#' @importFrom magrittr %>%
#' @examples
#' # Pick a radius to evaluate Ripley's K
#' r <- 30
#'
#' # Save the image IDs
#' ids <- unique(simdata$id)
#'
#' # Compute Ripley's K for tumor cells at r
#' K.vec <- c()
#' npoints.vec <- c()
#' for (i in 1:length(ids)) {
#'   # Save the ith image
#'   image.i <- simdata %>%
#'     dplyr::filter(id == ids[i]) %>%
#'     dplyr::select(x,y,type)
#'
#'   # Convert to a point process object
#'   w <- spatstat.geom::convexhull.xy(image.i$x, image.i$y)
#'  image.i.subset <- image.i %>% dplyr::filter(type == "a")
#'  image.ppp <- spatstat.geom::as.ppp(image.i.subset, W = w)
#'
#'   # Compute Kest
#'   Ki <- spatstat.explore::Kest(image.ppp, r = 0:30)
#'
#'   # Calculate the number of points
#'   npoints.vec[i] <- spatstat.geom::npoints(image.ppp)
#'
#'   # Save the result
#'   K.vec[i] <- Ki$iso[31]
#' }
#'
#' data.subset <- simdata %>%
#'   dplyr::select(tidyselect::all_of(c("id", "PID", "out"))) %>%
#'   dplyr::distinct()
#'
#' data.subset$spatial <- K.vec
#' data.subset$npoints <- npoints.vec
#'
#' # Remove NaNs
#' data.subset <- data.subset %>% dplyr::filter(!is.na(spatial))
#'
#' # Test
#' combo.weight.avg(data = data.subset, group = "PID", outcome = "out", model = "logistic")
combo.weight.avg <- function(data, group, outcome, cens = NULL, model = "survival", adjustments = NULL, n.ensemble = 1000, seed = NULL) {

  # Check for a seed
  if (is.null(seed)) {
    seed <- sample(1:.Machine$integer.max, 1)
  }

  set.seed(seed)

  # Create a p-value based on many weighted averages
  p.values <- sapply(1:n.ensemble, function(ii) {

    # Generate weights
    weights <- stats::rnorm(nrow(data))
    scaled.weights <- abs(weights)/sum(abs(weights))

    # Add as a column
    data$weights <- scaled.weights

    # Calculate the weighted average
    data <- data %>%
      dplyr::group_by(get(group)) %>%
      dplyr::mutate(weight.avg = sum(spatial*weights*npoints)/sum(weights*npoints)) %>%
      dplyr::ungroup() %>%
      dplyr::select(tidyselect::all_of(c(group, outcome, cens, adjustments, "weight.avg"))) %>%
      dplyr::distinct()

    # Fit outcome model
    if (model == "survival") {

      # Fit Cox model
      if (is.null(adjustments)) {
        res <- survival::coxph(survival::Surv(get(outcome), get(cens)) ~ weight.avg, data = data)
      }

      if (!is.null(adjustments)) {
        res <- survival::coxph(survival::Surv(get(outcome), get(cens)) ~ weight.avg + get(adjustments), data = data)
      }

      # Save the p-value
      p.value <- summary(res)$coefficients[1,5]

    }

    if (model == "logistic") {

      # Fit logistic model
      if (is.null(adjustments)) {
        res <- stats::glm(get(outcome) ~ weight.avg, data = data, family = stats::binomial())
      }

      if (!is.null(adjustments)) {
        res <- stats::glm(get(outcome) ~ weight.avg + get(adjustments), data = data, family = stats::binomial())
      }

      # Save the p-value
      p.value <- summary(res)$coefficients[2,4]

    }

    # Return the p-value
    p.value

  })

  # Combine p-values
  pval <- ACAT::ACAT(p.values)

  # Return
  list(pval = pval, seed = seed)
}
