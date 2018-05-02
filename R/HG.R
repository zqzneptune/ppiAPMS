#' HG: Calculates HG scoresa given input table of AP-MS runs
#'
#'
#'
#' @title HG
#' @param datInput A table of Run_id, Prey, Spectral counts, Protein Length.
#' @return A data frame containing the following columns:
#'         InteractorA, InteractorB, HG scores
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom dplyr n
#' @importFrom tidyr spread
#' @importFrom magrittr %>%

#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel parApply
#' @importFrom parallel stopCluster
#' @importFrom utils combn
#' @importFrom dplyr bind_rows
#' @importFrom stats setNames
#' @importFrom stats phyper
#'
#' @export
#' @author Qingzhou Zhang

HG <- function(datInput){
  . <- NULL
  Run_id <- NULL
  Peptide_cnt <- NULL
  Length <- NULL
  NormalSpec <- NULL
  SumNS <- NULL
  NSAF <- NULL
  NormalNSAF <- NULL
  Tn <- NULL
  UniprotID <- NULL
  tnA <- NULL
  tnB <- NULL
  NMinTn <- NULL
  HG <- NULL
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
  ppsTN <-
    parApply(cl, pps, 2, function(x){
      # sum(apply(rbind(g.list[[x[1]]], g.list[[x[2]]]), 2, min), na.rm = TRUE)
      gX <-
        g.list[[x[1]]]
      gY <-
        g.list[[x[2]]]
      sum(c(gX[gX < gY], gY[gY <= gX]), na.rm = TRUE)
    })
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
