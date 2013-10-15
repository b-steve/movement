set.seed(5434)
library(R2admb)
source("sim-movement.r")
n <- 20
kappa <- 60
a <- 20
b <- 15
sigma <- 10
locs <- sim.movement(n = n, kappa = kappa, a = a, b = b, sigma = sigma)

Xobs <- locs$Xobs[-1, ]
write_dat("movement", list(n = n, x_obs = Xobs))
write_pin("movement", list(kappa = kappa, a = a, b = b, sigma = sigma, x = Xobs))
compile_admb("movement", safe = TRUE, re = TRUE, verbose = TRUE)
run_admb("movement", verbose = TRUE, extra.args = "-noinit")
fit <- read_admb("movement")
clean_admb("movement")

fit$maxgrad
summary(fit)

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

## > summary(fit)
## Model file: movement 
## Negative log-likelihood:  105.9 	 AIC:  219.9 
## Coefficients:
##       Estimate Std. Error z value Pr(>|z|)    
## kappa   69.090     55.415   1.247 0.212478    
## a       20.629      2.373   8.691  < 2e-16 ***
## b       10.462      2.730   3.833 0.000127 ***
## sigma   10.163      1.619   6.276 3.46e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
