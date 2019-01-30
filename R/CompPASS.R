#' CompPASS
#' Comparative Proteomic Analysis Software Suite (CompPASS) is based on spoke
#' model. This algorithm was developed by Dr. Mathew Sowa for defining the human
#' deubiquitinating enzyme interaction landscape (Sowa, Mathew E., et al., 2009)
#' and later BioPlex 1.0 (Huttlin, Edward L., et al., 2015). The implementation
#' of this algorithm was inspired by Dr. Sowa's online tutorial
#'  (\url{http://besra.hms.harvard.edu/ipmsmsdbs/cgi-bin/tutorial.cgi}).
#' The output includes Z-score, S-score, D-score and WD-score.
#'
#' @title CompPASS
#' @param datInput A dataframe with column names: idRun, idBait, idPrey, countPrey.
#' Each row represent one unique protein captured in one pull-down experiment
#' @return A data frame consists of unique bait-prey pairs with Z-score, S-score,
#' D-score and WD-score indicating interacting probabilities.
#'
#' @author Qingzhou Zhang, \email{zqzneptune@hotmail.com}
#' @references Huttlin, Edward L., et al. "The BioPlex network: a systematic exploration of the human interactome." Cell 162.2 (2015): 425-440. \url{https://doi.org/10.1016/j.cell.2015.06.043}
#' @references Sowa, Mathew E., et al. "Defining the human deubiquitinating enzyme interaction landscape." Cell 138.2 (2009): 389-403. \url{https://doi.org/10.1016/j.cell.2009.04.042}
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom dplyr n
#' @importFrom tidyr spread
#' @importFrom magrittr %>%
#' @importFrom stats quantile
#' @importFrom stats sd
#' @export

CompPASS <- function(datInput){
  idBait <- NULL
  idPrey <- NULL
  countPrey <- NULL
  AvePSM <- NULL
  Stsc <- NULL
  Mtsc <- NULL
  f_sum <- NULL
  . <- NULL

  k <-
    length(unique(datInput$idBait))
  statsTbl <-
    datInput %>%
    group_by(`idBait`, `idPrey`) %>%
    summarise(`AvePSM` = mean(`countPrey`)) %>%
    mutate(`BP` = paste(`idBait`, `idPrey`, sep = "~"))
  stats <-
    spread(statsTbl[, c("idBait", "idPrey", "AvePSM")], `idBait`, `AvePSM`)
  statsTable <-
    as.matrix(stats[, -1])
  rownames(statsTable) <-
    stats$idPrey
  statsTable[is.na(statsTable)] <- 0
  z <-
    t(apply(statsTable, 1, function(x){(x-mean(x))/sd(x)}))
  zTable <-
    data.frame(`idBait` = unlist(lapply(colnames(z), function(x){rep(x, nrow(z))})),
               `idPrey` = rep(rownames(z), ncol(z)),
               `Z_score` = c(z),
               stringsAsFactors = FALSE) %>%
    mutate(`BP` = paste(`idBait`, `idPrey`, sep = "~"))
  w <-
    data.frame(`idPrey` = rownames(statsTable),
               `Mtsc` = apply(statsTable, 1, mean),
               `Stsc` = apply(statsTable, 1, sd),
               stringsAsFactors = FALSE) %>%
    mutate(`w` = ifelse(`Stsc`/`Mtsc` <=1, 1, `Stsc`/`Mtsc`))
  f <-
    unique(datInput[, c("idBait", "idPrey")]) %>%
    group_by(`idPrey`) %>%
    summarise(`f_sum` = n())
  p <-
    datInput[, c("idBait", "idPrey")] %>%
    group_by(`idBait`, `idPrey`) %>%
    summarise(`p` = n()) %>%
    mutate(`BP` = paste(`idBait`, `idPrey`, sep = "~"))
  scoreTbl <-
    statsTbl %>%
    left_join(., f, by = "idPrey") %>%
    left_join(., p[, c("BP", "p")], by = "BP") %>%
    left_join(., w, by = "idPrey") %>%
    mutate(`k` = k) %>%
    mutate(`S_score` = sqrt((`AvePSM`)*(`k`)/(`f_sum`))) %>%
    mutate(`D_score` = sqrt((`AvePSM`)*(((`k`)/(`f_sum`))^`p`))) %>%
    mutate(`WD_inner` = (`k` / `f_sum`) * (`Stsc` / `Mtsc`)) %>%
    mutate(`WD_score` =  sqrt((`AvePSM`)*(((`k`)/(`f_sum`)*`w`)^`p`))) %>%
    left_join(., zTable[, c("BP", "Z_score")], by = "BP")
  return(scoreTbl)
}
