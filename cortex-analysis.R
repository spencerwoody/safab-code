
library(saFAB)

library(here)
library(ggplot2)
library(dplyr)
library(readr)
library(latex2exp)
library(grid)
library(gridExtra)
library(cowplot)

theme_set(theme_half_open())

## ------------------------------------------------------------------------
## Read in data

sync <- read_csv("data/synchrony_smithkohn2008.csv")

## Center and scale data
z <- (sync$z - 0.61) / 0.81

t <- 2

## ------------------------------------------------------------------------
## Predictive recurison

mypr <- pr_fitter(z, nulltype = "empirical", sig0 = 1)

## Joint
spendfunj <- mypr$my_fit_joint

spendfunc <- make_w_theta(marginal_fun = mypr$my_fit_joint,
                          theta_min = -20, theta_max = 20, sigma = 1, t = 2, alpha = 0.10,
                          verbose = T, num_theta=1000)



spendfunplot <- spendfunc %>%
  ggplot() +
  ## geom_vline(xintercept = 0) + 
  geom_line(aes(theta, w)) +
  scale_x_continuous(limits = c(-5, 20)) +
  labs(## title = "Spending function for neural synchrony data",
       x = TeX("$\\theta$"),
       y = TeX("$w(\\theta)$"))

spendfunplot

ggsave("cortex-spendfun.pdf", spendfunplot,
       width = 7, height = 3, units = "in")

ggsave("cortex-spendfun.eps", spendfunplot,
       width = 7, height = 3, units = "in")

## ------------------------------------------------------------------------
## Create saFAB intervals



myList <- saFAB(spendfunc$theta,
                spendfunc$w,
                t = t, sigma = 1,
                yMin = -15, yMax = 15, yNum = 5e3, alpha = 0.10)

myListu <- saFAB(spendfunc$theta,
                 rep(0.5, length(spendfunc$theta)),
                 t = t, sigma = 1,
                 yMin = -15, yMax = 15, yNum = 5e3, alpha = 0.10)

## Confidence intervals for observed statistics
myz <- z[abs(z) > t]

## myz <- myz[abs(myz) < quantile(z, 0.99)]

mydf <- data.frame(
  z = myz,
  interval = I(vector("list", length(myz))),
  interval_u = I(vector("list", length(myz))),
  intervalLo = NA,
  intervalHi = NA,
  interval_uLo = NA,
  interval_uHi = NA,
  intervalLength = NA,
  intervalLengthMono = NA,
  intervalLength_u = NA,
  numIntervals = NA,
  error = NA,
  error_u = NA
)

## Get intervals for each z
for (k in 1:nrow(mydf)) {
  zidx <- which.min(abs(myList$Cdf$y - myz[k]))


  zidxu <- which.min(abs(myListu$Cdf$y - myz[k]))

  mydf$interval[[k]] <- myList$Cdf$intervals[[zidx]]
  mydf$intervalLength[k] <- myList$Cdf$intervalLength[zidx]

  mydf$intervalLengthMono[k] <- max(mydf$interval[[k]]) - min(mydf$interval[[k]])
  
  mydf$interval_u[[k]] <- myListu$Cdf$intervals[[zidxu]]
  mydf$intervalLength_u[k] <- myListu$Cdf$intervalLength[zidxu]

  ## Min and max of the interval
  mydf$intervalLo[k] <- min(mydf$interval[[k]])
  mydf$intervalHi[k] <- max(mydf$interval[[k]])

  mydf$interval_uLo[k] <- min(mydf$interval_u[[k]])
  mydf$interval_uHi[k] <- max(mydf$interval_u[[k]])
  
  ## Track number of intervals in each set
  mydf$numIntervals[k] <- myList$Cdf$numIntervals[zidx]

  mydf$error[k] <- myz[k] - myList$Cdf$y[zidx]
  mydf$error_u[k] <- myz[k] - myListu$Cdf$y[zidxu]
}

## Compare the lengths of the intervals
mean(mydf$intervalLength)
mean(mydf$intervalLength_u)

mean(mydf$intervalLength) / mean(mydf$intervalLength_u)

1 - mean(mydf$intervalLength) / mean(mydf$intervalLength_u)

mean(mydf$intervalLength < mydf$intervalLength_u)

ggplot() +
  geom_histogram(aes(mydf$intervalLength_u - mydf$intervalLength),
                 col = "snow", fill = "grey40")

## Plot the given intervals

neuralPlot <- mydf %>%
  arrange(z) %>%
  ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = -t, lty = "dashed") +
  geom_vline(xintercept = t, lty = "dashed") +
  geom_rug(aes(z, z), sides = "b", alpha = 0.25) + 
  ## geom_point(aes((z), z)) + 
  geom_point(aes((z), intervalLo,
                    col = "saFAB"), shape = 20) +
  geom_point(aes((z), intervalHi,
                    col = "saFAB"), shape = 20) +
  geom_point(aes((z), interval_uLo,
                    col = "saUMAU"), shape = 20) +
  geom_point(aes((z), interval_uHi,
                    col = "saUMAU"), shape = 20) +
  scale_color_manual(name = "method", values = c("black", "grey60")) +
  theme_half_open() + 
  theme(legend.position = "top") +
  labs(## title = "Confidence sets for neural synchrony dataset",
       x = "Observed test statistic", y = "Confidence set")

neuralPlot

myscale <- 1.1

ggsave("cortex-intervals.pdf", neuralPlot, width = 6 * myscale,
       height = 4 * myscale, units = "in")

ggsave("cortex-intervals.eps", neuralPlot, width = 6 * myscale,
       height = 4 * myscale, units = "in", device = cairo_ps)
