
## Number of iterations
nMC <- 1000


###############################################################################
                                        #           Define functions          #
###############################################################################

int_fun_exp_joint <- function(theta, mu, lambda, y, sigma) {
  dnorm(y, theta, sigma) * dexp(theta - mu, rate = lambda, log = FALSE)
}


marginal_exp_joint <- function (y, sigma, mu, p, lambda, t) {
  
  density <- (p * integrate(int_fun_exp_joint, lower = mu, upper = Inf,
                           mu = mu, lambda = lambda,
                           y = y, sigma = sigma)$value +
    (1 - p) * dnorm(y, 0, sigma)) * ifelse(abs(y) > t, 1, 0)

  return(density)
  
}

# confidence level
alpha <- 0.10

p <- 0.2

tau <- sqrt(3)

## Number of subjects
Nj <- 2000

## Number of observations per subject
Ni <- 3

## sigma <- 1 * sqrt(Ni - 1)
sigma <- 1 * sqrt(Ni)

## SD for ybar used for confidence intervals
## sigma_mean <- sigma / sqrt(Ni - 1)

sigma_mean <- sigma / sqrt(Ni)

# Acceptance region
tstar <- 2
## t <- tstar * sigma / sqrt(Ni - 1)

t <- tstar * sigma / sqrt(Ni)

## Number of folds

numFolds <- 5

foldid <- rep(1:numFolds, each = Nj / numFolds)

## Clipping
clip <- 0.1


# Confidence level
alpha <- 0.10

p <- 0.2


mu <- 1
lambda <- 1 
tau <- sqrt(3)

