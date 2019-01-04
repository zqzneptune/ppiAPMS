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
#' @importFrom utils combn
#' @importFrom dplyr bind_rows
#' @importFrom stats setNames
#' @importFrom stats phyper
#' @importFrom ppiAPMS GetPPN
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
  g[is.na(g)] <- 0
  rownames(g) <-
    d$Prey
  pps <-
    combn(d$Prey, 2)
  PPN <-
    GetPPN(t(g))
  CppPPN <-
    PPN[lower.tri(PPN,diag = F)]
  datPPI <-
    data.frame(cbind(t(pps[, CppPPN != 0]),
                     CppPPN[CppPPN != 0]),
               stringsAsFactors = FALSE)
  colnames(datPPI) <-
    c("InteractorA", "InteractorB", "ppiTN")
  datPPI$ppiTN <-
    as.numeric(datPPI$ppiTN)
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
      mutate(`HG` = -phyper(`ppiTN`, `tnA`, `NMinTn`-`tnB`,`tnB`,
                            lower.tail = FALSE, log.p = TRUE))
  return(scorePPI)
}
