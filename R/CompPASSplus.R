#' CompPASSplus: Calculates Z-Scores, WD scores, and entropy for a given input
#' table of AP-MS runs
#'
#'
#' @title CompPASSplus
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
#' @export
#' @author Qingzhou Zhang

CompPASSplus <- function(raw_dat){
  . <- NULL
  Prey <- NULL
  Bait <- NULL
  Run_id <- NULL
  Peptide_cnt <- NULL
  MaxTSC <- NULL
  nBait <- NULL
  f_sum <- NULL
  AvePSM <- NULL
  MeanDiff <- NULL
  Mean <- NULL
  SD <- NULL
  WD_inner <- NULL
  normFactor <- 0.98

  # Use total number of baits instead of actual AP-MS runs as
  # total number of experiments.
  # Multiple runs for the same bait were considered as replicated
  k <-
    length(unique(raw_dat$Bait))
  # f_sum: number of runs capturing the same prey
  f <-
    unique(raw_dat[, c("Bait", "Prey")]) %>%
    group_by(`Prey`) %>%
    summarise(`f_sum` = n())
  # p: number of replicates the bait-prey pair captured
  p <-
    raw_dat[, c("Bait", "Prey")] %>%
      group_by(`Bait`, `Prey`) %>%
      summarise(`p` = n()) %>%
      mutate(`BP` = paste(`Bait`, `Prey`, sep = "~"))
  # Using MAXIMUM spectral count to compute Entropy
  e <-
    raw_dat %>%
      group_by(`Bait`, `Prey`, `Run_id`) %>%
      summarise(`MaxTSC` = max(`Peptide_cnt`)) %>%
      mutate(`Entropy` = entropy(`MaxTSC`)) %>%
      mutate(`BP` = paste(`Bait`, `Prey`, sep = "~"))
  # AvePSM = average speactral counts across replicates.
  stats <-
    raw_dat %>%
      group_by(`Bait`, `Prey`) %>%
      summarise(`AvePSM` = mean(`Peptide_cnt`)) %>%
      mutate(`BP` = paste(`Bait`, `Prey`, sep = "~"))
  dat <-
    stats %>%
      left_join(., unique(e[, c("BP", "Entropy")]), by = "BP") %>%
      left_join(., f, by = "Prey") %>%
      left_join(., p[, c("BP", "p")], by = "BP") %>%
      mutate(`k` = k)

  preyWithBaitOnly <-
    stats %>%
      group_by(`Prey`) %>%
      summarise(`nBait` = n()) %>%
      filter(`nBait` == 1) %>%
      .$Prey
  dat <-
    dat %>%
    mutate(`f_num_no_prey` =
             ifelse((`Bait` == `Prey`)&(`Prey` %in% preyWithBaitOnly),
                    k, k-`f_sum`))
  df_num_no_prey <-
    unique(dat[, c("Prey", "f_num_no_prey")])
  f_num_no_prey <-
    df_num_no_prey$f_num_no_prey
  names(f_num_no_prey) <-
    df_num_no_prey$Prey
  statsM <-
    stats %>%
      mutate(`AvePSM` = ifelse(`Bait` != `Prey`,
                                `AvePSM`,
                                ifelse(`Prey` %in% preyWithBaitOnly,
                                       `AvePSM`,
                                       NA)))
  statsMatrix <-
    spread(`statsM`[, c("Bait", "Prey", "AvePSM")], `Bait`, `AvePSM`)
  m <-
    as.matrix(statsMatrix[, -1])
  rownames(m) <-
    statsMatrix$Prey

  prey.mean <-
    apply(m, 1, function(x){
      prey.mean <-
        sum(x, na.rm = TRUE)/k
    })

  prey.stats <-
    data.frame(`Prey` = names(prey.mean),
               `Mean` = prey.mean,
               `MeanDiff` = rowSums((m - prey.mean)^2, na.rm = TRUE),
               `f_num_no_prey` = f_num_no_prey[names(prey.mean)],
               stringsAsFactors = FALSE) %>%
    mutate(`SD` =  sqrt((`MeanDiff` + ((`Mean`^2) * (`f_num_no_prey`)))/(k - 1)))
  avePSM <-
    left_join(`dat`, prey.stats[, c("Prey", "Mean", "SD")], by = "Prey")
  output <-
    avePSM %>%
    mutate(`Z_score` = (`AvePSM` - `Mean`) / (`SD`)) %>%
    mutate(`S_score` = sqrt((`AvePSM`) * (`k`) / (`f_sum`))) %>%
    mutate(`D_score` = sqrt((`AvePSM`) * (((`k`) / (`f_sum`))^`p`))) %>%
    mutate(`WD_inner` = (`k` / `f_sum`) * (`SD` / `Mean`)) %>%
    mutate(`WD_raw` = sqrt(`AvePSM` * (`WD_inner`^`p`)))
  output[, "WD_score"] <-
    normalize.wd(output$WD_raw, normFactor)
  output <-
    as.data.frame(output)
  return(output)
}

## Calculates Shannon entropy for a list of values. Because the log of zero is
## undefined, a fractional pseudocount is added to each value. This pseudocount
## is set to 1/ # of values.
## Orginated from  devtools::install_github("dnusinow/cRomppass")

## @param xs A vector of values to calculate the entropy for
## @return The calculated entropy value

entropy <- function(xs) {
  p <- (xs + 1/length(xs)) / (sum(xs) + 1)
  ent <-
    sum(sapply(p, function(x){
      -1*x*log(x, 2)
    }))
  return(ent)
}
## Normalizes a vector of WD scores for a given normalization factor
## Orginated from devtools::install_github("dnusinow/cRomppass")
## @param xs A vector of unnormalized WD scores
## @param norm.factor A number between 0 and 1 corresponding to the
##    quantile of the xs to normalize to. Defaults in the comppass function
##    to the 98th percentile.
## @return A vector of WD scores divided by whatever score is the percentile
##    specified by norm.factor

normalize.wd <- function(xs, norm.factor) {
  xs.f <-
    Filter(function(x){ !(is.nan(x) || is.na(x)) }, xs)
  return(xs / quantile(xs.f, norm.factor)[1])
}
