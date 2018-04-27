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
  normFactor <-
    0.98
  k <-
    length(unique(raw_dat$Bait))
  f <-
    unique(raw_dat[, c("Bait", "Prey")]) %>%
    group_by(`Prey`) %>%
    summarise(`f_sum` = n())
  p <-
    raw_dat[, c("Bait", "Prey")] %>%
    group_by(`Bait`, `Prey`) %>%
    summarise(`p` = n()) %>%
    mutate(`BP` = paste(`Bait`, `Prey`, sep = "~"))

  e <-
    raw_dat %>%
    group_by(`Bait`, `Prey`, `Run_id`) %>%
    summarize(`MaxTSC` = max(`Peptide_cnt`)) %>%
    mutate(`Entropy` = entropy(`MaxTSC`)) %>%
    mutate(`BP` = paste(`Bait`, `Prey`, sep = "~"))
  stats <-
    raw_dat %>%
    group_by(`Bait`, `Prey`) %>%
    summarize(`AvePSM` = mean(`Peptide_cnt`)) %>%
    mutate(`BP` = paste(`Bait`, `Prey`, sep = "~"))
  dat <-
    stats %>%
    left_join(., unique(e[, c("BP", "Entropy")]), by = "BP") %>%
    left_join(., f, by = "Prey") %>%
    left_join(., p[, c("BP", "p")], by = "BP") %>%
    mutate(`k` = k)
  library(parallel)
  numCores <- detectCores()
  cl <- makeCluster(numCores)
  clusterEvalQ(cl, {
    library(tidyverse)
  })
  # clusterExport(cl, "dat")
  # clusterExport(cl, "k")
  pp <-
    parApply(cl, dat, 1, function(x){
      d <-
        dat %>%
        filter((`Bait` != x[2])&(`Prey` == x[2]))
      f_num_no_prey <-
        k-unique(d$f_sum)
      if(nrow(d) == 0){
        d <-
          dat %>%
          filter((`Bait` == x[2])&(`Prey` == x[2]))
        f_num_no_prey <-
          k
      }
      prey.sum <-
        sum(d$AvePSM)
      prey.mean <-
        prey.sum/k
      prey.sum.squared.err <-
        sum((d$AvePSM - prey.mean)^2) + ((prey.mean^2) * (f_num_no_prey))
      prey.sd <- sqrt(prey.sum.squared.err / (k - 1))
      return(c("SumAPSM" = prey.sum,
               "Mean" = prey.mean,
               "SD" = prey.sd))
    })
  stopCluster(cl)
  avePSM <-
    bind_cols(dat, data.frame(t(pp), stringsAsFactors = FALSE))
  output <-
    avePSM %>%
    mutate(`Z_score` = (`AvePSM` - `Mean`)/(`SD`)) %>%
    mutate(`S_score` = sqrt((`AvePSM`)*(`k`)/(`f_sum`))) %>%
    mutate(`D_score` = sqrt((`AvePSM`)*(((`k`)/(`f_sum`))^`p`))) %>%
    mutate(`WD_inner` = (`k` / `f_sum`) * (`SD` / `Mean`)) %>%
    mutate(`WD_raw` = sqrt(`AvePSM` * (`WD_inner`^`p`)))
  output[, "WD_score"] <- normalize.wd(output$`WD_raw`, normFactor)
  return(output)
}
