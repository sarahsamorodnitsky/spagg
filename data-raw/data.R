# Generate example data for `spagg` R package

# Packages
library(spatstat)

# Set a seed
set.seed(68490)

# Fix parameters
n <- 20 # Number of samples to generate

# Initialize data.frame to store samples
data <- data.frame(
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
                           cell.id = 1:xy$n,
                           x = xy$x,
                           y = xy$y,
                           type = sample(c("a", "b"), xy$n, replace = TRUE))

    # Save the outcome
    if (i <= 10) {
      image.ij$out <- 1
    } else {
      image.ij$out <- 0
    }

    data <- rbind.data.frame(data, image.ij)
  }
}

# Save the resulting dataset
usethis::use_data(data, overwrite = TRUE)
