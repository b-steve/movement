set.seed(5439)
library(R2admb)
source("sim-movement.r")
## True parameter values
n <- 20
kappa <- 6
a <- 20
b <- 15
sigma <- 10
## Simulate the data
locs <- sim.movement(n = n, kappa = kappa, a = a, b = b, sigma = sigma)
Xobs <- locs$Xobs[-1, ]
## Generate ADMB .dat and .pin files
write_dat("movement", list(n = n, x_obs = Xobs))
## Set start values here
write_pin("movement", list(kappa = 6, a = 25, b = 15, sigma = 15, x = Xobs))
## No need to compile yourself
##compile_admb("movement", safe = TRUE, re = TRUE, verbose = TRUE)
run_admb("movement", verbose = FALSE, extra.args = "-noinit")
fit <- read_admb("movement")
clean_admb("movement")
## Look at gradient and fit summary
fit$maxgrad
summary(fit)

## Plotting location point estimates (red) with observed locations
## (black) and actual locations (blue)
coords <- coef(fit, type = "random")
xs <- coords[c(TRUE, FALSE)]
ys <- coords[c(FALSE, TRUE)]
xlim <- range(c(range(sapply(locs, function(x) range(x[, 1]))), range(xs)))
ylim <- range(c(sapply(locs, function(x) range(x[, 2])), range(ys)))
plot(locs$X, type = "l", xlim = xlim, ylim = ylim, asp = 1, col = "blue")
points(locs$X, col = "blue")
lines(locs$Xobs)
points(locs$Xobs)
lines(xs, ys, col = "red")
points(xs, ys, col = "red")

## Model summary:
##
## Model file: movement
## Negative log-likelihood:  116.1 	 AIC:  224.2
## Coefficients:
##       Estimate Std. Error z value Pr(>|z|)
## kappa   11.126      6.574   1.693  0.09055 .
## a       15.532      1.859   8.354  < 2e-16 ***
## b        7.772      2.680   2.900  0.00373 **
## sigma   11.780      1.615   7.296 2.96e-13 ***
