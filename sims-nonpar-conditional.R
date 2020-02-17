
library(rootSolve)
library(saFAB)
library(dplyr)

source("sims_create_intervals.R")

jobJ <- system("echo $LAUNCHER_JID", intern = TRUE)

jobJpad <- stringr::str_pad(as.numeric(jobJ), 3, "0", side = "left")

Date <- lubridate::today()

filename <- sprintf("output/nonpar-conditional/nonpar-conditional-%s.Rdata", Date)

###############################################################################
                                        #            Set parameters           #
###############################################################################

source("setpars-nonpar-conditional.R")

thetaVec <- rss(Nj, p, tau)

###############################################################################
                                        #         Perform simulations         #
###############################################################################

yArray <- array(dim = c(Nj, Ni, nMC))
thetaMat <- matrix(nrow = Nj, ncol = nMC)

ebList <- vector("list", nMC)

prList <- vector("list", nMC)

wListE <- vector("list", nMC)
wListN <- vector("list", nMC)

simsDfList <- vector("list", nMC)

meancoverVecO <- rep(NA, nMC)
meancoverVecC <- rep(NA, nMC)
meancoverVecE <- rep(NA, nMC)
meancoverVecN <- rep(NA, nMC)
meancoverVecU <- rep(NA, nMC)

meanlengthVecO <- rep(NA, nMC)
meanlengthVecC <- rep(NA, nMC)
meanlengthVecE <- rep(NA, nMC)
meanlengthVecN <- rep(NA, nMC)
meanlengthVecU <- rep(NA, nMC)

## pre-cache the oracle and UMAU intervalsn
wsaFABo <- make_w_theta(marginal_fun = marginal_ss_conditional, 
                        sigma = sigma_mean, t = t, p = p, tau = tau, 
                        alpha = alpha,  epsilon=1e-10,
                        theta_min = -7, theta_max = 7, num_theta = 500, 
                        verbose = T)

## Clipped spending function
wsaFABc <- wsaFABo

wsaFABc$w[wsaFABc$w < clip] <- clip
wsaFABc$w[wsaFABc$w > 1 - clip] <- 1 - clip

## plot(wsaFABo$theta, wsaFABo$w, type = "l")



oList <- saFAB(wsaFABo$theta, wsaFABo$w, t, sigma_mean,
               yMin = -5, yMax = 5, yNum = 1e3)

cList <- saFAB(wsaFABc$theta, wsaFABc$w, t, sigma_mean,
               yMin = -5, yMax = 5, yNum = 1e3)

uList <- saFAB(wsaFABo$theta, rep(0.5, length(wsaFABo$theta)), t, sigma_mean,
               yMin = -5, yMax = 5, yNum = 1e3)


## Loop through simulations
prog <- progress_estimated(nMC)

cat("starting simulations...\n")

k <- 1

while (k <= nMC) {
  
  thetaVecK <- thetaVec
  yMatK <- matrix(rnorm(Nj * Ni, thetaVecK, sigma), nrow = Nj)
  
  thetaMat[, k] <- thetaVecK
  yArray[, , k] <- yMatK

  prList[[k]] <- vector("list", numFolds)
  wListE[[k]] <- vector("list", numFolds)
  wListN[[k]] <- vector("list", numFolds)

  ## ebList[[k]] <- myeb

  simsDf <- vector("list", numFolds)

  foldCount <- 1

  while(foldCount <= numFolds) {

    ## print(foldCount)

    ySpend <- yMatK[foldid != foldCount, ]
    yInt <- yMatK[foldid == foldCount, ]

    thetaInt <- thetaVecK[foldid == foldCount]

    mypr <- try(
    pr_fitter_nosplit(ySpend, nulltype = "empirical", sig0 = sigma_mean)
    )

    if (inherits(mypr, "try-error")) {
      break
    }

    prList[[k]][[foldCount]] <- mypr

    myeb <- eb_fun_knownsigma(rowMeans(ySpend), sigma = sigma_mean)

    ebList[[k]][[foldCount]] <- myeb

    p_hat <- myeb$p_hat
    tau_hat <- myeb$tau_hat

    wsaFABe <- make_w_theta(marginal_fun = marginal_ss_conditional,
                            sigma = sigma_mean, t = t, p = p_hat,
                            tau = tau_hat, alpha = alpha,
                            epsilon=1e-10, theta_min = -7,
                            theta_max = 7, num_theta = 500,
                            verbose = FALSE)

    wsaFABn <- make_w_theta(marginal_fun = mypr$my_fit_cnd,
                            sigma = sigma_mean, t = t, alpha = alpha,
                            epsilon=1e-5, theta_min = -7,
                            theta_max = 7, num_theta = 500,
                            verbose = FALSE, sigma_orig = sigma,
                            n = Ni)

    eList <- saFAB(wsaFABe$theta, wsaFABe$w, t, sigma_mean, verbose = FALSE,
                   yMin = -5, yMax = 5, yNum = 1e3)

    ## BP 5
    nList <- saFAB(wsaFABn$theta, wsaFABn$w, t, sigma_mean, verbose = FALSE,
                   yMin = -5, yMax = 5, yNum = 1e3)

    ## Store the spending functions in lists
    wListE[[k]][[foldCount]] <- wsaFABe
    wListN[[k]][[foldCount]] <- wsaFABn

    ## Fit samples to intervals
    yMeanK <- rowMeans(yInt)
    selectIdx <- which(abs(yMeanK) > t)
    ySelect <- yMeanK[selectIdx]
    thetaSelect <- thetaInt[selectIdx]

    simsDf[[foldCount]] <- sims_create_intervals(k, thetaSelect,
                                                 ySelect, selectIdx,
                                                 t, oList, cList, eList,
                                                 nList) %>%
      mutate(fold = foldCount)



    foldCount <- foldCount + 1

  }

  if (inherits(mypr, "try-error")) {
    next
  }

  simsDf <- plyr::rbind.fill(simsDf)

  simsDfList[[k]] <- simsDf
  
    ## Coverage check
  meancoverVecO[k] <- (simsDf$inInterval_o) %>% mean()
  meancoverVecC[k] <- (simsDf$inInterval_c) %>% mean()
  meancoverVecE[k] <- simsDf$inInterval_e %>% mean()
  meancoverVecN[k] <- simsDf$inInterval_n %>% mean()
  meancoverVecU[k] <- simsDf$inInterval_u %>% mean()
  
  ## Length
  meanlengthVecO[k] <- simsDf$intervalLength_o %>% mean()
  meanlengthVecC[k] <- simsDf$intervalLength_c %>% mean()
  meanlengthVecE[k] <- simsDf$intervalLength_e %>% mean()
  meanlengthVecN[k] <- simsDf$intervalLength_n %>% mean()
  meanlengthVecU[k] <- simsDf$intervalLength_u %>% mean()
  
  ## Print progress report
  prog$tick()$print()
 
  ## if (k %% 10 == 0) cat(sprintf("\n%i out of %i simulations...\n", k, nMC))
  
  cat(sprintf("\n%i out of %i simulations at %s...\n", k, nMC, lubridate::now()))

  k <- k + 1

}

simsDfCat <- plyr::rbind.fill(simsDfList) %>% mutate(job = jobJ)

save(list=ls(), file = filename)

load("output/nonpar-conditional/nonpar-conditional-2019-10-15.Rdata")

nMC
Ni
Nj

meancoverVecO %>% mean()
meancoverVecC %>% mean()
meancoverVecE %>% mean()
meancoverVecN %>% mean()
meancoverVecU %>% mean()

## sqrt(meancoverVecO * (1 - meancoverVecO) / nMC)
## sqrt(meancoverVecC * (1 - meancoverVecC) / nMC)
## sqrt(meancoverVecE * (1 - meancoverVecE) / nMC)
## sqrt(meancoverVecN * (1 - meancoverVecN) / nMC)
## sqrt(meancoverVecU * (1 - meancoverVecU) / nMC)

mean(meanlengthVecO)
mean(meanlengthVecC) 
mean(meanlengthVecE) 
mean(meanlengthVecN) 
mean(meanlengthVecU) 

mean(meanlengthVecO) / mean(meanlengthVecU)
mean(meanlengthVecC) / mean(meanlengthVecU)
mean(meanlengthVecE) / mean(meanlengthVecU)
mean(meanlengthVecN) / mean(meanlengthVecU)
mean(meanlengthVecU) / mean(meanlengthVecU)

###############################################################################
                                        #                xtable                #
###############################################################################

meancoverO <- meancoverVecO %>% mean()
meancoverE <- meancoverVecE %>% mean()
meancoverN <- meancoverVecN %>% mean()
meancoverU <- meancoverVecU %>% mean()

## SE for average coverage
meancoverSE_O <- sqrt(mean(meancoverVecO) * (1 - mean(meancoverVecO)) / nMC)
meancoverSE_E <- sqrt(mean(meancoverVecE) * (1 - mean(meancoverVecE)) / nMC)
meancoverSE_N <- sqrt(mean(meancoverVecN) * (1 - mean(meancoverVecN)) / nMC)
meancoverSE_U <- sqrt(mean(meancoverVecU) * (1 - mean(meancoverVecU)) / nMC)

## Average mean length
meanlengthO <- mean(meanlengthVecO) 
meanlengthE <- mean(meanlengthVecE) 
meanlengthN <- mean(meanlengthVecN) 
meanlengthU <- mean(meanlengthVecU) 

## SE of Average mean length
meanlengthSE_O <- sd(meanlengthVecO) / sqrt(nMC) 
meanlengthSE_E <- sd(meanlengthVecE) / sqrt(nMC) 
meanlengthSE_N <- sd(meanlengthVecN) / sqrt(nMC) 
meanlengthSE_U <- sd(meanlengthVecU) / sqrt(nMC) 

## Relative length
meanrellengthO <- mean(meanlengthVecO) / mean(meanlengthVecU)
meanrellengthE <- mean(meanlengthVecE) / mean(meanlengthVecU)
meanrellengthN <- mean(meanlengthVecN) / mean(meanlengthVecU)
meanrellengthU <- mean(meanlengthVecU) / mean(meanlengthVecU)

## SE for relative length
meanrellengthSE_O <- 1 / mean(meanlengthVecU) * sd(meanlengthVecO) / sqrt(nMC) 
meanrellengthSE_E <- 1 / mean(meanlengthVecU) * sd(meanlengthVecE) / sqrt(nMC) 
meanrellengthSE_N <- 1 / mean(meanlengthVecU) * sd(meanlengthVecN) / sqrt(nMC) 
meanrellengthSE_U <- 1 / mean(meanlengthVecU) * sd(meanlengthVecU) / sqrt(nMC) 

## Combine into summary, print to console and LaTeX

summaryDf <- data.frame(
  Method = c("Oracle", "PEB", "NPEB", "UMAU"),
  Coverage = c(
    sprintf("%1.4f (%1.4f)", meancoverO, meancoverSE_O),
    sprintf("%1.4f (%1.4f)", meancoverE, meancoverSE_E),
    sprintf("%1.4f (%1.4f)", meancoverN, meancoverSE_N),
    sprintf("%1.4f (%1.4f)", meancoverU, meancoverSE_U)
  ),
  `Average size` = c(
    sprintf("%1.4f (%1.4f)", meanlengthO, meanlengthSE_O),
    sprintf("%1.4f (%1.4f)", meanlengthE, meanlengthSE_E),
    sprintf("%1.4f (%1.4f)", meanlengthN, meanlengthSE_N),
    sprintf("%1.4f (%1.4f)", meanlengthU, meanlengthSE_U)
  ),
  Rel_average_size = c(
    sprintf("%1.4f (%1.4f)", meanrellengthO, meanrellengthSE_O),
    sprintf("%1.4f (%1.4f)", meanrellengthE, meanrellengthSE_E),
    sprintf("%1.4f (%1.4f)", meanrellengthN, meanrellengthSE_N),
    sprintf("%1.4f (%1.4f)", meanrellengthU, meanrellengthSE_U)
  )
)

summaryDfXtable <- xtable(summaryDf, floating = FALSE,
                          booktabs = TRUE,
                          caption = "Well specified prior case, conditional selection")

print(summaryDfXtable, include.rownames=FALSE)
