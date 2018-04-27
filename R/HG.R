HG <- function(datInput){
  library(tidyverse)
  library(parallel)
  datCnt <-
    datInput %>%
    mutate(`NormalSpec` = `Peptide_cnt`/`Length`) %>%
    group_by(`Run_id`) %>%
    mutate(`SumNS` = sum(`NormalSpec`)) %>%
    mutate(`NSAF` = `NormalSpec`/`SumNS`) %>%
    group_by(`Run_id`) %>%
    mutate(`NormalNSAF` = `NSAF`/min(`NSAF`)) %>%
    mutate(`Tn` = as.integer(sqrt(`NormalNSAF`)))
  d <-
    spread(datCnt[, c("Run_id", "Prey", "Tn")], `Run_id`, `Tn`)
  g <-
    as.matrix(d[, -1])
  rownames(g) <-
    d$Prey
  g.list <-
    setNames(split(g, seq(nrow(g))), rownames(g))
  names(g.list) <-
    rownames(g)
  preys <-
    rep(0, length(names(g.list)))
  names(preys) <-
    names(g.list)
  pps <-
    combn(names(preys), 2)
  numCores <- detectCores()
  cl <- makeCluster(numCores)
  # clusterExport(cl, "g.list")
  system.time(
    ppsTN <-
      parApply(cl, pps, 2, function(x){
        sum(apply(rbind(g.list[[x[1]]], g.list[[x[2]]]), 2, min), na.rm = TRUE)
      })
  )
  stopCluster(cl)
  ppiTN <-
    ppsTN[ppsTN != 0]
  ppi <-
    pps[, ppsTN != 0]
  datPPI <-
    data.frame(cbind(t(ppi), ppiTN), stringsAsFactors = FALSE) %>%
    mutate(`ppiTN` = as.numeric(`ppiTN`))
  colnames(datPPI) <-
    c("InteractorA", "InteractorB", "ppiTN")
  tnInteractorA <-
    datPPI[, c("InteractorA", "ppiTN")]
  colnames(tnInteractorA) <-
    c("UniprotID", "ppiTN")
  tnInteractorB <-
    datPPI[, c("InteractorB", "ppiTN")]
  colnames(tnInteractorB) <-
    c("UniprotID", "ppiTN")
  tnProtein <-
    bind_rows(tnInteractorA, tnInteractorB) %>%
    group_by(`UniprotID`) %>%
    summarise(`minTn` = sum(`ppiTN`))
  sumMinTnInteractorA <-
    tnProtein
  colnames(sumMinTnInteractorA) <-
    c("InteractorA", "tnA")
  sumMinTnInteractorB <-
    tnProtein
  colnames(sumMinTnInteractorB) <-
    c("InteractorB", "tnB")
  scorePPI <-
    datPPI %>%
    left_join(., sumMinTnInteractorA, by = "InteractorA") %>%
    left_join(., sumMinTnInteractorB, by = "InteractorB") %>%
    mutate(`NMinTn` = sum(tnProtein$minTn)/2) %>%
    mutate(`HG` = -phyper(`ppiTN`, `tnA`, `NMinTn`-`tnB`,`tnB`, lower.tail=FALSE, log.p=TRUE))
  return(scorePPI)
}
