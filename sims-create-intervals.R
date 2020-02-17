
sims_create_intervals <- function(k, thetaSelect, ySelect, selectIdx, t,
                                  oList, cList, eList, nList) {

  simsDf <- data.frame(
    trial = k,
    theta = thetaSelect,
    y = ySelect,
    intervals_o = I(vector("list", length(selectIdx))),
    intervals_c = I(vector("list", length(selectIdx))),
    intervals_e = I(vector("list", length(selectIdx))),
    intervals_n = I(vector("list", length(selectIdx))),
    intervals_u = I(vector("list", length(selectIdx))),
    inInterval_o = NA,
    inInterval_e = NA,
    inInterval_n = NA,
    inInterval_u = NA,
    intervalLength_o = NA,
    intervalLength_e = NA,
    intervalLength_n = NA,
    intervalLength_u = NA
  )

  for (j in 1:nrow(simsDf)) {
    
    yIdx <- which.min(abs(simsDf$y[j] - oList$Cdf$y))
    
    simsDf$intervals_o[[j]] <- oList$Cdf$intervals[yIdx]
    simsDf$intervals_c[[j]] <- cList$Cdf$intervals[yIdx]
    simsDf$intervals_e[[j]] <- eList$Cdf$intervals[yIdx]
    simsDf$intervals_n[[j]] <- nList$Cdf$intervals[yIdx]
    ## simsDf$intervals_u[[j]] <- uList$Cdf$intervals[yIdx]
    simsDf$intervals_u[[j]] <- sa_umau_intervals(simsDf$y[j], thetaRange = c(-20, 20),
                                                 sd = sigma_mean, t = t)
    
    simsDf$inInterval_o[j] <- inRange(simsDf$theta[j], simsDf$intervals_o[[j]][[1]])
    simsDf$inInterval_c[j] <- inRange(simsDf$theta[j], simsDf$intervals_c[[j]][[1]])
    simsDf$inInterval_e[j] <- inRange(simsDf$theta[j], simsDf$intervals_e[[j]][[1]])
    simsDf$inInterval_n[j] <- inRange(simsDf$theta[j], simsDf$intervals_n[[j]][[1]])
    ## simsDf$inInterval_u[j] <- inRange(simsDf$theta[j], simsDf$intervals_u[[j]][[1]])
    simsDf$inInterval_u[j] <- inRange(simsDf$theta[j], simsDf$intervals_u[[j]])
    
    simsDf$intervalLength_o[j] <- oList$Cdf$intervalLength[yIdx]
    simsDf$intervalLength_c[j] <- cList$Cdf$intervalLength[yIdx]
    simsDf$intervalLength_e[j] <- eList$Cdf$intervalLength[yIdx]
    simsDf$intervalLength_n[j] <- nList$Cdf$intervalLength[yIdx]
    simsDf$intervalLength_u[j] <- uList$Cdf$intervalLength[yIdx]  

    ## print(j)

  }

  simsDf

}
