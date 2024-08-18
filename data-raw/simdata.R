# Generate example data for `spagg` R package

# Packages
library(spatstat)

# Set a seed
set.seed(68490)

# Fix parameters
n <- 20 # Number of samples to generate

# Initialize data.frame to store samples
simdata <- data.frame(
  PID = numeric(),
  id = character(),
  cell.id = numeric(),
  x = numeric(),
  y = numeric(),
  type = character(),
  out = numeric()
)

# Randomly select how many images to generate per sample
n.image.per.PID <- sample(x = 1:5, size = n, replace = TRUE)

# Iterate through the PIDs
for (i in 1:n) {

  # Simulate the outcome
  if (i <= 10) {
    out <- rbinom(1, size = 1, p = 0.75)
  } else {
    out <- rbinom(1, size = 1, p = 0.25)
  }

  # Iterate through the images
  for (j in 1:n.image.per.PID[i]) {

    # Randomly choose a point pattern to simulate
    flip <- sample(0:1, 1)

    # Generate a point pattern --

    # Fix the window
    window <- owin(xrange = c(0, 100), yrange = c(0, 100))

    # Poisson process
    if (flip == 0) {
      xy <- spatstat.random::rpoispp(lambda = 0.1, win = window)
    }

    # Clusters
    if (flip == 1 & i <= 10) {
      xy <- spatstat.random::rMatClust(kappa = 0.001, scale = 10, mu = 25, win = window)
    }

    # Repulsed
    if (flip == 1 & i > 10) {
      xy <- spatstat.random::rStrauss(beta = 0.01, gamma = 0, R = 5, W = window)
    }

    # Save the image
    image.ij <- data.frame(PID = i,
                           id = paste0("PID.", i, ".image.", j),
                           x = xy$x,
                           y = xy$y)

    # Remove duplicates
    image.ij <- image.ij %>% dplyr::distinct()

    # Add cell IDs and labels
    image.ij$cell.id <- 1:xy$n
    image.ij$type<- sample(c("a", "b"), xy$n, replace = TRUE)

    # Save the outcome
    image.ij$out <- out

    # Add in the image information
    simdata <- rbind.data.frame(simdata, image.ij)
  }
}

# Save the resulting dataset
usethis::use_data(simdata, overwrite = TRUE)
