set.seed(5434)
library(R2admb)
source("sim-movement.r")
n <- 20
kappa <- 3.5
a <- 20
b <- 15
sigma <- 10
locs <- sim.movement(n = n, kappa = kappa, a = a, b = b, sigma = sigma)

Xobs <- locs$Xobs[-1, ]
write_dat("movement", list(n = n, x_obs = Xobs))
write_pin("movement", list(kappa = kappa, a = a, b = b, sigma = sigma, x = Xobs))
##compile_admb("movement", safe = TRUE, re = TRUE, verbose = TRUE)
run_admb("movement", verbose = FALSE, extra.args = "-noinit")
fit <- read_admb("movement")
clean_admb("movement")

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
## Negative log-likelihood:  119.7 	 AIC:  231.4
## Coefficients:
##       Estimate Std. Error z value Pr(>|z|)
## kappa    6.882      3.979   1.730   0.0837 .
## a       17.698      2.595   6.819 9.16e-12 ***
## b       11.051      2.642   4.183 2.88e-05 ***
## sigma   11.532      1.839   6.272 3.57e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
