
# saFAB intervals vs UMAU intervals for joint selection
# Construct confidence intervals along


library(saFAB)

library(here)
library(ggplot2)
library(latex2exp)
library(grid)
library(gridExtra)
library(rootSolve)
library(cowplot)

theme_set(theme_half_open())


# Set parameters ----------------------------------------------------------

# Parameters
p <- 0.1
sigma <- 1

tau <- sqrt(3)

# Acceptance region
t <- 2

# Confidence level
alpha <- 0.10

# Create spending function ------------------------------------------------

wsaFAB <- make_w_theta(marginal_fun = marginal_ss_joint, 
                       sigma = sigma, t = t, p = p, tau = tau, 
                       alpha = alpha,  epsilon=1e-10,
                       theta_min = -7, theta_max = 7, num_theta = 1000, 
                       verbose = T)

spendplot <- wsaFAB %>% ggplot() +
  geom_line(aes(theta, w, col = "saFAB"), size = 0.65) +
  geom_hline(aes(yintercept = 0.5, col = "saUMAU"), size = 0.65) +
  labs(## title = "Spending function",
       x = TeX("$\\theta$"),
       y = TeX("$w(\\theta)$")) +
  ## scale_color_manual(name = "method", values = c("black", "grey80")) +
  scale_color_manual(name = "method", values = c("black", "grey65")) +
  theme(legend.position = "top") 

spendplot

ggsave("toyspendplot.pdf", width=6, height = 3.5, units = "in")

ggsave("toyspendplot.eps", width=6, height = 3.5, units = "in")

## ------------------------------------------------------------------------
## Create saFAB intervals

saFABlist <- saFAB(theta = wsaFAB$theta, w = wsaFAB$w, t = t, sigma =
                   sigma, alpha = alpha, yMin = -5, yMax= 5, yNum =
                   5000, verbose = T)

UMAUlist <- saFAB(theta = wsaFAB$theta, w = rep(0.5, length(wsaFAB$theta)), t = t,
                  sigma = sigma, alpha = alpha, yMin = -5, yMax= 5,
                  yNum = 5000, verbose = T)

tdf <- data.frame(
  t = c(-t, t),
  Name = "selection\nregion"
)

## Plot acceptance intervals
Aplot <- ggplot() + 
  geom_hline(data = tdf, aes(yintercept = t, lty = Name)) +
  geom_line(data = saFABlist$Adf,
            aes(theta, A, group = ul, col = "saFAB")) +
  geom_line(data = UMAUlist$Adf,
            aes(theta, A, group = ul, col = "UMAU")) + 
  scale_linetype_manual(name="", values = "dotted") + 
  scale_color_manual(name = "method", values = c("grey10", "grey65")) +
  ## scale_color_manual(name = "method", values = c("black", "grey75")) +
  ## scale_color_brewer(palette = "Greys") + 
  labs(title = "Acceptance regions", 
       x = TeX("$\\theta$"), 
       y = TeX("$A(\\theta)$")) +
  theme(legend.position = "bottom")

Aplot

## Plot confidence intervals


Cplot <- saFABlist$CdfPlotting %>%
  filter(abs(y) > t) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(data = tdf, aes(xintercept = t, lty = Name)) + 
  geom_linerange(data = UMAUlist$CdfPlotting %>% filter(abs(y) > t),
                 aes(y, ymin = lower, ymax = upper, col = "UMAU"),
                 alpha = 0.1) +
  geom_linerange(aes(y, ymin = lower, ymax = upper, col = "saFAB"),
                 alpha = 0.1) +
  ## scale_color_manual(name = "method", values = c("grey5", "grey60")) +
  scale_color_manual(name = "method", values = c("grey10", "grey75")) +
  scale_linetype_manual(name="", values = "dashed") + 
  labs(title = "Confidence sets for observed y", x = "y", y = "C(y)") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = "bottom")

Cplot

# Concatenate the plots

###############################################################################
                                        #           Rearrange Cplot           #
###############################################################################

## Cdf1 <- saFABlist$Cdf %>%
##   filter(y < -t) %>%
##   mutate(yId = 1:n())

## CdfPlot1 <- saFABlist$CdfPlotting %>%
##   filter(y < -t) %>%
##   left_join(Cdf1 %>% select(y, yId))

## Cdf1 %>% pull(y) %>% range()
## CdfPlot1 %>% pull(y) %>% range()

## firstDivIdx <- NA
## firstDivFound <- FALSE
## lastDivIdx <- NA
## lastDivFound <- FALSE

## for (i in 1:nrow(Cdf1)) {

##   if (Cdf1$numIntervals[i] != 1) {
##     if (!firstDivFound) {
##       firstDivIdx <- i
##       firstDivFound <- TRUE
##     }
##   }

##   if (firstDivFound & Cdf1$numIntervals[i] == 1 & !lastDivFound) {
##     lastDivIdx <- i - 1
##     lastDivFound <- TRUE
##   }

## }

## yVec <- rep(NA, nrow(Cdf1))
## upperVec <- rep(NA, nrow(Cdf1))
## lowerVec <- rep(NA, nrow(Cdf1))

## for (i in 1:nrow(Cdf1)) {
##   if (i <= firstDivIdx) {
    
##   }
## }

## ------------------------------------------------------------------------
## Split up into two

## grob1 <- arrangeGrob(spendplot, Aplot)

## grob2 <- arrangeGrob(Aplot, Cplot)

## grid.arrange(grob1)

## grid.arrange(grob2)

## ggsave("mercury/toy-spend-accept.pdf", grob1)

## ggsave("mercury/toy-accept-conf.pdf", grob2)

## ------------------------------------------------------------------------
## Resume...

## myGrob <- arrangeGrob(spendplot, Aplot, Cplot)

## grid.arrange(spendplot, Aplot, Cplot)

## ggsave(here("tex/fig/triplot.pdf"), plot = myGrob, width = 8, height = 10)

# Simulations: draw from prior, and sampling distribution -----------------

numSamples <- 1e4

thetaVec <- rss(numSamples, p, tau)
  
# IMPORTANT: This is where the joint selection happens...

simsDf <- data.frame(
  theta = thetaVec,
  y = NA,
  intervals_s = I(vector("list", numSamples)),
  intervals_u = I(vector("list", numSamples)),
  intervalLength_s = NA,
  intervalLength_u = NA,
  inInterval_s = NA,
  inInterval_u = NA,
  bayes_l = NA,
  bayes_u = NA,
  sabayes_l = NA,
  sabayes_u = NA
)

prog <- progress_estimated(numSamples)

for (k in 1:numSamples) {

  while (TRUE) {
    simsDf$theta[k] <- rss(1, p, tau)
    my_y <- rnorm(1, simsDf$theta[k], sigma)
    
    if (abs(my_y) > t) break
    
  }
  simsDf$y[k] <- my_y
  
  # Index for generated y
  yIdx <- which.min(abs(simsDf$y[k] - UMAUlist$Cdf$y))

  ## Bayesian credible interval
  simsDf$bayes_l[k] <- bayes_find_quantile(simsDf$y[k], p, sigma, tau, alpha / 2)
  simsDf$bayes_u[k] <- bayes_find_quantile(simsDf$y[k], p, sigma, tau, 1 - alpha / 2)
  
  ## saBayesian credible interval
  simsDf$sabayes_l[k] <- sabayes_find_quantile(simsDf$y[k], p, sigma, tau, t, alpha / 2)
  simsDf$sabayes_u[k] <- sabayes_find_quantile(simsDf$y[k], p, sigma, tau, t, 1 - alpha / 2)


  # Fill in interval and interval length from index
  simsDf$intervals_s[[k]] <- saFABlist$Cdf$intervals[[yIdx]]
  simsDf$intervals_u[[k]] <- UMAUlist$Cdf$intervals[[yIdx]]
  
  simsDf$intervalLength_s[k] <- saFABlist$Cdf$intervalLength[yIdx]
  simsDf$intervalLength_u[k] <- UMAUlist$Cdf$intervalLength[yIdx]
  
  # Check if theta is in the interval
  simsDf$inInterval_s[k] <- inRange(simsDf$theta[k], simsDf$intervals_s[[k]])
  simsDf$inInterval_u[k] <- inRange(simsDf$theta[k], simsDf$intervals_u[[k]])
    
  prog$tick()$print()
}

# Analyze simulations -----------------------------------------------------

simsDf <- simsDf %>%
  mutate(numIntervals_s = sapply(simsDf$intervals_s, length))

simsDf$intervals_s

simsDf$intervals_s %>% sapply(rev) %>% sapply("[[", 1)


## simsDf <- simsDf%>%
##   mutate(
##     s_lo = intervals_s %>%
##       sapply("[[", 1),
##     s_hi = intervals_s %>%
##       sapply(rev) %>%
##       sapply("[[", 1))

## simsDf <- simsDf%>%
##   mutate(
##     s_lo = intervals_s %>%
##       sapply("[[", 1),
##     s_hi = intervals_s %>%
##       sapply(rev) %>%
##       sapply("[[", 1),
##     intervalLength_s_x = s_hi - s_lo,
##     eff_lost = (intervalLength_s_x - intervalLength_s) / intervalLength_s)

## mean(simsDf$eff_lost == 0)

## mean(simsDf$intervalLength_s)

## simsDf %>% filter(eff_lost != 0) %>% pull(eff_lost) %>% mean()

## sum(sapply(simsDf$intervals_s, length) == 4)

## # Check for correct coverage
## simsDf$inInterval_s %>% mean()
## simsDf$inInterval_u %>% mean()

## # Check lengths
## mean(simsDf$intervalLength_s) 
## mean(simsDf$intervalLength_u) 

## (mean(simsDf$intervalLength_u) - mean(simsDf$intervalLength_s)) / 
##   mean(simsDf$intervalLength_u)

## (mean(simsDf$intervalLength_u) - mean(simsDf$intervalLength_s_x)) / 
##   mean(simsDf$intervalLength_u)

## Diff <- abs(saFABlist$Cdf$intervalLength - UMAUlist$Cdf$intervalLength)

## minDiffIdx <- which.min(Diff * ifelse(abs(saFABlist$Cdf$y) > t, 1, 100000))

xint <- abs(saFABlist$Cdf$y[minDiffIdx])



# plot confidence interval length
lengthPlot <- ggplot() +
  geom_line(data = saFABlist$Cdf %>%
              filter(abs(y) > t) %>%
              mutate(y = abs(y)),
            aes(y, intervalLength, col = "saFAB"), size = 1) + 
  geom_line(data = UMAUlist$Cdf %>%
              filter(abs(y) > t) %>%
              mutate(y = abs(y)),
            aes(y, intervalLength, col = "saUMAU"), size = 1) +
  geom_vline(xintercept = xint, lty = "dashed") +
  scale_color_manual(name = "method",
                     values = c("grey10", "grey65")) +
  labs(## title = "Comparison of interval sizes, joint selection setting",
       x = "y", y = "Size of confidence set") + 
  theme(## legend.position = c(0.9, 0.25),
    legend.position = "top", 
    axis.title.x = element_blank())+
  guides(col = guide_legend(keywidth = 2))

lengthPlot


## Plot length vs. density and cumulative density
emp <- ecdf(abs(simsDf$y))

simsDftransform <- simsDf %>%
  mutate(y = abs(y),
         emp = emp(y))

cdf_scale <- 1.5

yHist <- simsDftransform %>% 
  ggplot() + 
  geom_histogram(aes(y, ..density..), fill = "white", col = "grey25", bins = 30) +
  geom_line(aes(y, emp * cdf_scale), size = 0.75) +
  geom_vline(xintercept = xint, lty = "dashed") +
  ## geom_rug(data = simsDftransform %>% filter(numIntervals_s == 4),
  ##          aes(y, 1), sides = "b", alpha = 0.25) + 
  scale_x_continuous(limits = c(t, max(saFABlist$Cdf$y))) + 
  scale_y_continuous(sec.axis = sec_axis(~./cdf_scale, 
                                         name = "cumulative density",
                                         breaks = seq(0, 1, by = 0.5))) +
  labs(x = "|y|")
yHist

pdf(here("length-hist-joint.pdf"), width = 7, height = 5.5)
grid.newpage()
grid.draw(rbind(ggplotGrob(lengthPlot), ggplotGrob(yHist), size = "last"))
dev.off()



## ## ------------------------------------------------------------------------
## ## For Bayesian stuff

## ## Bayesian
## simsDf$bayes_l[k] <- bayes_find_quantile(simsDf$y[k], p, sigma, tau, alpha / 2)
## simsDf$bayes_u[k] <- bayes_find_quantile(simsDf$y[k], p, sigma, tau, 1 - alpha / 2)
  
## ## saBayesian
## simsDf$sabayes_l[k] <- sabayes_find_quantile(simsDf$y[k], p, sigma, tau, t, alpha / 2)
## simsDf$sabayes_u[k] <- sabayes_find_quantile(simsDf$y[k], p, sigma, tau, t, 1 - alpha / 2)

