# ppiAPMS
This is an R package that implements statistical modeling approaches that assign confidence scores to protein-protein interaction data generated using affinity purification–mass spectrometry (AP-MS) data.

## Download

The development version can be installed through github:

    devtools::install_github(repo="zqzneptune/ppiAPMS")
    library(ppiAPMS)

## Quick start

1. CompPASS and CompPASS-Plus

Summarize your AP-MS data from proteome database search into the following format:

|idRun|idBait|idPrey|countPrey|
|-----|:----:|:----:|:-------:|
|Unique ID of affinity run|Bait ID|Prey ID|Peptide count|
