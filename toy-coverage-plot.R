
library(saFAB)
library(here)
library(dplyr)
library(ggplot2)
library(tidyr)
library(latex2exp)
library(grid)
library(gridExtra)
library(cowplot)

theme_set(theme_half_open(font_size = 13))

## ------------------------------------------------------------------------
## Set parameters


p <- 0.1
sigma <- 1

tau <- sqrt(3)

# Acceptance region
t <- 2

# Confidence level
alpha <- 0.10

## Critical value for this confidence level
zstar <- qnorm(1 - alpha / 2)

## ------------------------------------------------------------------------
## Create saFAB intervals

wsaFAB <- make_w_theta(marginal_fun = marginal_ss_joint,
                       sigma = sigma, t = t, p = p, tau = tau,
                       alpha = alpha,  epsilon=1e-10, theta_min = -10,
                       theta_max = 10, num_theta = 1000, verbose = T)

sList <- saFAB(wsaFAB$theta, wsaFAB$w, t, sigma,
               yMin = -8, yMax = 8, yNum = 1e4)

## ------------------------------------------------------------------------
## Grid of theta

numReps <- 1000

thetaVec <- c(rep(0, ceiling((1 - p) / p)),
              seq(0.001, 5, length.out = 100))




intervalDf <- data.frame(
  theta = rep(thetaVec, each = numReps),
  y = NA,
  fLo = NA,
  fHi = NA,
  bLo = NA,
  bHi = NA,
  sIntervals = I(vector("list", length(thetaVec) * numReps)),
  inInterval_f = NA,
  inInterval_b = NA,
  inInterval_s = NA
)

## ------------------------------------------------------------------------
## Simulate data for each theta, and intervals for each one as well

prog <- progress_estimated(nrow(intervalDf))

for (j in 1:nrow(intervalDf)) {

  ## Current theta
  thetaJ <- intervalDf$theta[j]

  ## Generate a y for this theta
  yJ <- qtnorm(runif(1), thetaJ, sigma, t)
 
  ## Unadjuted frequentist interval
  intervalDf$fLo[j] <- yJ - zstar * sigma
  intervalDf$fHi[j] <- yJ + zstar * sigma

  ## Bayesian credible interval
  intervalDf$bLo[j] <- sabayes_find_quantile(yJ, p, sigma, tau, t,
                                             q = alpha / 2)
  intervalDf$bHi[j] <- sabayes_find_quantile(yJ, p, sigma, tau, t,
                                             q = 1 - alpha / 2)

  ## saFAB interval
  yJidx <- which.min(abs(yJ - sList$Cdf$y))
  
  intervalDf$sIntervals[[j]] <- sList$Cdf$intervals[yJidx]

  ## Check for coverage
  intervalDf$inInterval_f[j] <- (thetaJ >= intervalDf$fLo[j] &
                                 thetaJ <= intervalDf$fHi[j]) * 1
  
  intervalDf$inInterval_b[j] <- (thetaJ >= intervalDf$bLo[j] &
                                 thetaJ <= intervalDf$bHi[j]) * 1

  intervalDf$inInterval_s[j] <- inRange(thetaJ, intervalDf$sIntervals[[j]][[1]])

  prog$tick()$print()
  
}

save(list = ls(), file = "cover-sims.Rdata")

load("cover-sims.Rdata")

intervalDf <- intervalDf %>%
  mutate(inInterval_b = (theta >= bLo & theta <= bHi) * 1)

##
##

intervalDf$inInterval_s %>% mean()

intervalDf$inInterval_b %>% mean()

intervalDf$inInterval_f %>% mean()

intervalDf$inInterval_s %>% table()

summaryDf <- intervalDf %>%
  group_by(theta) %>%
  summarize(coverS = mean(inInterval_s),
            coverB = mean(inInterval_b),
            coverF = mean(inInterval_f))

## ------------------------------------------------------------------------
## Plot

## theta = 0

summaryDfZero <- summaryDf %>%
  filter(theta == 0) %>%
  rename("saBayes" = coverB,
         "non-sa UMAU" = coverF,
         "saFAB" = coverS) %>% 
  gather(method, coverage, -theta)

coverplotZero <- summaryDfZero %>%
  filter(theta == 0) %>%
  ggplot() +
  geom_hline(yintercept = 1 - alpha, lty = "dotted", size = 1) +
  geom_col(aes(method, coverage, fill = method), col = "grey90") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  ## scale_fill_manual(name = "method",
  ##                   values = c("grey20", "springgreen4", "dodgerblue2")) +
  scale_fill_manual(name = "method",
                    values = rev(c("black", "grey45", "grey75"))) + 
  labs(title = TeX("Coverage rates for $\\theta = 0$"),
       ## subtitle = TeX(sprintf("$\\alpha = %1.2f$", alpha)),
       x = "method",
       y = "Coverage") +
  theme(legend.position = "none")

coverplotZero

ggsave("comparison_col.pdf", coverplotZero,
       device = "pdf", width = 4, height = 4, units = "in")

## theta != 0
coverplotNonzero <- summaryDf %>%
  filter(theta != 0) %>%
  ggplot() +
  geom_hline(yintercept = 1 - alpha, lty = "dotted", size = 1) + 
  geom_smooth(aes(theta, coverS, col = "\nsaFAB\n",
                  lty = "\nsaFAB\n"),
              se = FALSE, size = 1.4) +
  geom_smooth(aes(theta, coverB, col = "\nsaBayes\n",
                  lty = "\nsaBayes\n"),
              se = FALSE, size = 1.4) +
  geom_smooth(aes(theta, coverF, col = "\nnon-sa\nUMAU\n",
                  lty = "\nnon-sa\nUMAU\n"),
              se = FALSE, size = 1.4) +
  ## geom_point(aes(theta, coverS, col = "\nsaFAB\n")) +
  ## geom_point(aes(theta, coverB, col = "\nsaBayes\n")) +
  ## geom_point(aes(theta, coverF, col = "\nnon-sa\nUMAU\n")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  ## scale_color_manual(name = "method",
  ##                    values = c("grey20", "springgreen4", "dodgerblue2")) +
  scale_color_manual(name = "method",
                    values = rev(c("black", "grey45", "grey75"))) + 
  scale_linetype_manual(name = "method", values = c(1, 4, 2)) + 
  labs(title = TeX("Coverage rates for $\\theta \\neq 0$"),
       ## subtitle = TeX(sprintf("$\\alpha = %1.2f$", alpha)),
       x = TeX("|$\\theta$|"),
       y = "Coverage") +
  guides(color = guide_legend(keywidth = 2)) 

coverplotNonzero

ggsave("comparison_line.pdf", coverplotNonzero,
       device = "pdf", width = 6, height = 4, units = "in")


## ------------------------------------------------------------------------
## Combine the plots

mygrob <- arrangeGrob(coverplotZero, coverplotNonzero, nrow = 1, widths = c(3, 6),
                      heights = 4)

grid.arrange(mygrob)

myscale <- 1.12

ggsave("compare-coverage2.pdf", mygrob,
       width = 9 * myscale, height = 4 * myscale)

ggsave("compare-coverage2.eps", mygrob,
       width = 9 * myscale, height = 4 * myscale)
