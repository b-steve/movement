source("../R/movement-funs.r")
## True parameter values
n <- 200
kappa <- 6
a <- 20
b <- 15
sigma <- 15
## Simulate the data
set.seed(5143)
locs <- sim.movement(n = n, kappa = kappa, a = a, b = b, sigma = sigma)
## Pull out observed and actual locations (first row is always (0, 0),
## so removing it)
Xobs <- locs$Xobs[-1, ]
X <- locs$X[-1, ]
## Set up list of starting parameters
sv <- list(kappa = kappa, a = a, b = b, sigma = sigma)
## Fit the model
fit <- fit.movement(Xobs, sv, dir = "../admb")
summary(fit)
## Plot estimated locations
plot.movement(fit, X = X)
