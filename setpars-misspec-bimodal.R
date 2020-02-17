
marginal_nl_joint <- function(y, sigma, t, p, prop, mu, tau) {

  density <- (p * (prop *       dnorm(y, -mu, sqrt(sigma^2 + tau^2)) +
                   (1 - prop) * dnorm(y,  mu, sqrt(sigma^2 + tau^2))) +
                (1 - p) * dnorm(y, 0, sigma) ) *
    ifelse(abs(y) > t, 1, 0)

  return(density)

}

marginal_nll_joint <- function(theta, y, t, p, prop, mu, tau) {
  density <-
    p * prop * pdf_laplace(theta + mu, tau) +
    p * (1 - prop) * pdf_laplace
}




## Number of iterations
nMC <- 1000


p <- 0.2
sigma <- 1 

tau <- sqrt(3)

prop <- 0.50


# Confidence level
alpha <- 0.10


## ------------------------------------------------------------------------
## Generate data



## Number of subjects
Nj <- 2000

## thetaVec <- rss(Nj, p, tau)

## Number of observations per subject
Ni <- 3

sigma <- 1 * sqrt(Ni)

## yMat <- matrix(rnorm(Nj * Ni, thetaVec, sigma), nrow = Nj)
## yVec <- rnorm(Nj, thetaVec, sigma)

# Acceptance region
tstar <- 2

t <- tstar * sigma / sqrt(Ni)

mu <- 3
tau <- 1/2


sigma_mean <- sigma / sqrt(Ni)



numFolds <- 5

foldid <- rep(1:numFolds, each = Nj / numFolds)

## Clipping
clip <- 0.1


