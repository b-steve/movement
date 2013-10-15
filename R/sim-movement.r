#' Simulating locations
#'
#' Simulates movement model data.
#'
#' Note that original animal location is (0, 0), and that original bearing is also 0.
#'
#' @param n Number of points to simulate.
#' @param kappa Concentration parameter for the Von-Mises distribution.
#' @param a Mean step length for the gamma distribution.
#' @param b Step length standard deviation for the gamma distribution.
#' @param e Measurement error standard deviation for the bivariate normal distribution.
#' @return A list with two components: \code{X} is a matrix that contains actual animal
#' locations, \code{Xobs} contains observed locations.

sim.movement <- function(n, kappa, a, b, sigma){
  library(CircStats)
  X <- matrix(0, nrow = n + 1, ncol = 2)
  ## Start location is (0, 0); Original bearing is 0
  x <- c(0, 0)
  X[1, ] <- x
  phi <- 0
  for (i in 2:(n + 1)){
    ## Calculate new bearing
    mu <- rvm(1, mean = pi, k = kappa) - pi
    phi <- phi + mu
    ## Adjust for bearings not in (0, 2*pi)
    if (phi > 2*pi){
      phi <- phi - 4*pi
    } else if (phi < -2*pi){
      phi <- phi + 4*pi
    }
    ## Calculate step size
    len <- rgamma(1, shape = a^2/b^2, scale = b^2/a)
    ## Calculate new location
    x.delta <- c(len*sin(phi), len*cos(phi))
    x <- x + x.delta
    X[i, ] <- x
  }
  error <- matrix(rnorm(2*(n + 1), mean = 0, sd = sigma), nrow = n + 1, ncol = 2)
  Xobs <- X + error
  list(X = X, Xobs = Xobs)
}

