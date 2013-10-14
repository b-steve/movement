#' Simulating locations
#'
#' Simulates movement model data.
#'
#' @param n Number of points to simulate.
#' @param kappa Concentration parameter for the Von-Mises distribution.
#' @param a Mean stop length for the gamma distribution.
#' @param b Step length standard deviation for the gamma distribution.
#' @param e Measurement error standard deviation for the bivariate normal distribution.
#' @param phi0 Initial bearing.

sim.movement <- function(n, kappa, a, b, e, phi0){
  library(CircStats)
  X <- matrix(0, nrow = n, ncol = 2)
  x <- c(0, 0)
  phi <- phi0
  for (i in 1:n){
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
  error <- matrix(rnorm(2*n, mean = 0, sd = e), nrow = n, ncol = 2)
  Xobs <- X + error
  list(X = X, Xobs = Xobs)
}

locs <- sim.movement(n = 100, kappa = 3, a = 20, b = 8, e = 15, phi0 = 0)
xlim <- range(sapply(locs, function(x) range(x[, 1])))
ylim <- range(sapply(locs, function(x) range(x[, 2])))
plot(locs$X, type = "n", xlim = xlim, ylim = ylim, asp = 1, col = "blue")
points(locs$X, col = "blue")
lines(locs$Xobs)
points(locs$Xobs)
