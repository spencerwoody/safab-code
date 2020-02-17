
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

## Number of iterations
nMC <- 1000
