#' CompPASSplus
#' Implemented a naive Bayes classifier that learns to distinguish true
#'  interacting proteins from non-specific background and false positive
#'  identifications on the basis of CompPASS scores.
#'  The source code for this function was based on the source code. \url{https://github.com/dnusinow/cRomppass}
#'
#'
#' @title CompPASSplus
#' @param datInput A dataframe with column names: idRun, idBait, idPrey, countPrey.
#' Each row represent one unique protein captured in one pull-down experiment
#' @return A data frame consists of unique bait-prey pairs with Z-score, S-score,
#' D-score and WD-score indicating interacting probabilities.

#' @author Qingzhou Zhang, \email{zqzneptune@hotmail.com}
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
#' @export


CompPASSplus <- function(datInput){
  . <- NULL
  idPrey <- NULL
  idBait <- NULL
  idRun <- NULL
  countPrey <- NULL
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
    length(unique(datInput$idBait))
  # f_sum: number of runs capturing the same prey
  f <-
    unique(datInput[, c("idBait", "idPrey")]) %>%
    group_by(`idPrey`) %>%
    summarise(`f_sum` = n())
  # p: number of replicates the bait-prey pair captured
  p <-
    datInput[, c("idBait", "idPrey")] %>%
      group_by(`idBait`, `idPrey`) %>%
      summarise(`p` = n()) %>%
      mutate(`BP` = paste(`idBait`, `idPrey`, sep = "~"))
  # Using MAXIMUM spectral count to compute Entropy
  e <-
    datInput %>%
      group_by(`idBait`, `idPrey`, `idRun`) %>%
      summarise(`MaxTSC` = max(`countPrey`)) %>%
      mutate(`Entropy` = entropy(`MaxTSC`)) %>%
      mutate(`BP` = paste(`idBait`, `idPrey`, sep = "~"))
  # AvePSM = average speactral counts across replicates.
  stats <-
    datInput %>%
      group_by(`idBait`, `idPrey`) %>%
      summarise(`AvePSM` = mean(`countPrey`)) %>%
      mutate(`BP` = paste(`idBait`, `idPrey`, sep = "~"))
  dat <-
    stats %>%
      left_join(., unique(e[, c("BP", "Entropy")]), by = "BP") %>%
      left_join(., f, by = "idPrey") %>%
      left_join(., p[, c("BP", "p")], by = "BP") %>%
      mutate(`k` = k)

  preyWithBaitOnly <-
    stats %>%
      group_by(`idPrey`) %>%
      summarise(`nBait` = n()) %>%
      filter(`nBait` == 1) %>%
      .$idPrey
  dat <-
    dat %>%
    mutate(`f_num_no_prey` =
             ifelse((`idBait` == `idPrey`)&(`idPrey` %in% preyWithBaitOnly),
                    k, k-`f_sum`))
  df_num_no_prey <-
    unique(dat[, c("idPrey", "f_num_no_prey")])
  f_num_no_prey <-
    df_num_no_prey$f_num_no_prey
  names(f_num_no_prey) <-
    df_num_no_prey$idPrey
  statsM <-
    stats %>%
      mutate(`AvePSM` = ifelse(`idBait` != `idPrey`,
                                `AvePSM`,
                                ifelse(`idPrey` %in% preyWithBaitOnly,
                                       `AvePSM`,
                                       NA)))
  statsMatrix <-
    spread(`statsM`[, c("idBait", "idPrey", "AvePSM")], `idBait`, `AvePSM`)
  m <-
    as.matrix(statsMatrix[, -1])
  rownames(m) <-
    statsMatrix$idPrey

  prey.mean <-
    apply(m, 1, function(x){
      prey.mean <-
        sum(x, na.rm = TRUE)/k
    })

  prey.stats <-
    data.frame(`idPrey` = names(prey.mean),
               `Mean` = prey.mean,
               `MeanDiff` = rowSums((m - prey.mean)^2, na.rm = TRUE),
               `f_num_no_prey` = f_num_no_prey[names(prey.mean)],
               stringsAsFactors = FALSE) %>%
    mutate(`SD` =  sqrt((`MeanDiff` + ((`Mean`^2) * (`f_num_no_prey`)))/(k - 1)))
  avePSM <-
    left_join(`dat`, prey.stats[, c("idPrey", "Mean", "SD")], by = "idPrey")
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
## Orginated from https://github.com/dnusinow/cRomppass

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
## Orginated from https://github.com/dnusinow/cRomppass
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
