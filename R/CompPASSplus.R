entropy <- function(xs) {
  p <- (xs + 1/length(xs)) / (sum(xs) + 1)
  ent <- sum(sapply(p, function(x) { -1*x*log(x, 2) }))
  return(ent)
}
normalize.wd <- function(xs, norm.factor) {
  xs.f <- Filter(function(x) { ! (is.nan(x) || is.na(x)) }, xs)
  return(xs / quantile(xs.f, norm.factor)[1])
}
CompPASSplus <- function(raw_dat){
  library(tidyverse)
  normFactor <- 0.98
  k <-
    length(unique(raw_dat$Bait))
  f <-
    unique(raw_dat[, c("Bait", "Prey")]) %>%
    group_by(Prey) %>%
    summarise(f_sum = n())
  p <-
    raw_dat[, c("Bait", "Prey")] %>%
    group_by(Bait, Prey) %>%
    summarise(p = n()) %>%
    mutate(BP = paste(Bait, Prey, sep = "~"))
  e <-
    raw_dat %>%
    group_by(Bait, Prey, Run_id) %>%
    summarize(MaxTSC = max(Peptide_cnt)) %>%
    mutate(Entropy = entropy(MaxTSC)) %>%
    mutate(BP = paste(Bait, Prey, sep = "~"))
  stats <-
    raw_dat %>%
    group_by(Bait, Prey) %>%
    summarize(AvePSM = mean(Peptide_cnt)) %>%
    mutate(BP = paste(Bait, Prey, sep = "~"))
  dat <-
    stats %>%
    left_join(., unique(e[, c("BP", "Entropy")]), by = "BP") %>%
    left_join(., f, by = "Prey") %>%
    left_join(., p[, c("BP", "p")], by = "BP") %>%
    mutate(k = k)

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
    mutate(`AvePSM` = ifelse(`Bait` != `Prey`, `AvePSM`,
                             ifelse(`Prey` %in% preyWithBaitOnly, `AvePSM`, NA)))
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
    left_join(dat, prey.stats[, c("Prey", "Mean", "SD")], by = "Prey")
  output <-
    avePSM %>%
    mutate(Z_score = (AvePSM - Mean)/(SD)) %>%
    mutate(S_score = sqrt((AvePSM) * (k)/(f_sum))) %>%
    mutate(D_score = sqrt((AvePSM) * (((k)/(f_sum))^p))) %>%
    mutate(WD_inner = (k/f_sum) * (SD/Mean)) %>%
    mutate(WD_raw = sqrt(AvePSM * (WD_inner^p)))
  output[, "WD_score"] <- normalize.wd(output$WD_raw, normFactor)
  return(output)
}
