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

#' Fitting a movement model
#'
#' Fits a movement model.
#'
#' Parameters to fit are kappa, a, b, and sigma. Note that original animal location
#' (0, 0) and bearing 0 are assumed.
#'
#' @param Xobs A matrix of observed animal locations.
#' @param sv A list of start values. Each component corresponds to a parameter and
#' is named as such.
#' @param dir Directory containing admb movement executable.
fit.movement <- function(Xobs, sv, dir = "."){
  library(R2admb)
  currwd <- getwd()
  setwd(dir)
  n <- nrow(Xobs)
  write_dat("movement", list(n = n, x_obs = Xobs))
  sv$x <- Xobs
  write_pin("movement", sv)
  run_admb("movement", extra.args = "-noinit -l1 1000000000 -l2 1000000000 -l3 1000000000 -nl1 1000000000 -nl2 100000000 -cbs 1000000")
  fit <- read_admb("movement")
  clean_admb("")
  file.remove("movement.dat", "movement.pin")
  fit$Xobs <- Xobs
  setwd(currwd)
  fit
}

#' Plotting a movement model
#'
#' Plots estimated locations from a movement model.
#'
#' @param A fitted movement model.
#' @param Actual animal locations (if known).
plot.movement <- function(fit, X = NULL){
  coords <- coef(fit, type = "random")
  xs <- coords[c(TRUE, FALSE)]
  ys <- coords[c(FALSE, TRUE)]
  if (!is.null(X)){
    Xs <- X[, 1]
    Ys <- X[, 2]
  } else {
    Xs <- NULL
    Ys <- NULL
  }
  Xobs <- fit$Xobs
  obsXs <- Xobs[, 1]
  obsYs <- Xobs[, 2]
  xlim <- range(c(xs, Xs, obsXs))
  ylim <- range(c(ys, Ys, obsYs))
  plot(Xobs, type = "l", xlim = xlim, ylim = ylim, asp = 1,
       col = "lightgrey", xlab = "x", ylab = "y")
  if (!is.null(X)){
    lines(X, col = "blue")
  }
  lines(xs, ys, col = "red")
}
