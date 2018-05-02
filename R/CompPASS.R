#' CompPASS: Calculates Z-Scores, WD scores, and entropy for a given input
#' table of AP-MS runs
#'
#'
#' @title CompPASS
#' @param raw_dat A table of Run_id, Bait, Prey, and Spectral counts.
#' @return A data frame containing the following columns:
#'         Bait, Prey, AvePSM, Entropy, Z, S, D, WD scores
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
#' @author Qingzhou Zhang

CompPASS <- function(raw_dat){
  Bait <- NULL
  Prey <- NULL
  Peptide_cnt <- NULL
  AvePSM <- NULL
  Stsc <- NULL
  Mtsc <- NULL
  f_sum <- NULL
  . <- NULL

  k <-
    length(unique(raw_dat$Bait))
  statsTbl <-
    raw_dat %>%
    group_by(`Bait`, `Prey`) %>%
    summarise(`AvePSM` = mean(`Peptide_cnt`)) %>%
    mutate(`BP` = paste(`Bait`, `Prey`, sep = "~"))
  stats <-
    spread(statsTbl[, c("Bait", "Prey", "AvePSM")], `Bait`, `AvePSM`)
  statsTable <-
    as.matrix(stats[, -1])
  rownames(statsTable) <-
    stats$Prey
  statsTable[is.na(statsTable)] <- 0
  z <-
    t(apply(statsTable, 1, function(x){(x-mean(x))/sd(x)}))
  zTable <-
    data.frame(`Bait` = unlist(lapply(colnames(z), function(x){rep(x, nrow(z))})),
               `Prey` = rep(rownames(z), ncol(z)),
               `Z_score` = c(z),
               stringsAsFactors = FALSE) %>%
    mutate(`BP` = paste(`Bait`, `Prey`, sep = "~"))
  w <-
    data.frame(`Prey` = rownames(statsTable),
               `Mtsc` = apply(statsTable, 1, mean),
               `Stsc` = apply(statsTable, 1, sd),
               stringsAsFactors = FALSE) %>%
    mutate(`w` = ifelse(`Stsc`/`Mtsc` <=1, 1, `Stsc`/`Mtsc`))
  f <-
    unique(raw_dat[, c("Bait", "Prey")]) %>%
    group_by(`Prey`) %>%
    summarise(`f_sum` = n())
  p <-
    raw_dat[, c("Bait", "Prey")] %>%
    group_by(`Bait`, `Prey`) %>%
    summarise(`p` = n()) %>%
    mutate(`BP` = paste(`Bait`, `Prey`, sep = "~"))
  scoreTbl <-
    statsTbl %>%
    left_join(., f, by = "Prey") %>%
    left_join(., p[, c("BP", "p")], by = "BP") %>%
    left_join(., w, by = "Prey") %>%
    mutate(`k` = k) %>%
    mutate(`S_score` = sqrt((`AvePSM`)*(`k`)/(`f_sum`))) %>%
    mutate(`D_score` = sqrt((`AvePSM`)*(((`k`)/(`f_sum`))^`p`))) %>%
    mutate(`WD_inner` = (`k` / `f_sum`) * (`Stsc` / `Mtsc`)) %>%
    mutate(`WD_score` =  sqrt((`AvePSM`)*(((`k`)/(`f_sum`)*`w`)^`p`))) %>%
    left_join(., zTable[, c("BP", "Z_score")], by = "BP")
  return(scoreTbl)
}
