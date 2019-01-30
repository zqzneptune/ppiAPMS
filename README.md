# ppiAPMS
This is an R package that implements statistical modeling approaches that assign confidence scores to protein-protein interaction data generated using affinity purification–mass spectrometry (AP-MS) data.

## Installation

The development version can be installed through github:

    devtools::install_github(repo="zqzneptune/ppiAPMS")
    library(ppiAPMS)

## Quick start

1. CompPASS and CompPASS-Plus

Summarize your AP-MS data from proteome database search into the dataframe *datInput* with the following format:

|idRun|idBait|idPrey|countPrey|
|-----|:----:|:----:|:-------:|
|Unique ID of one affinity run|Bait ID|Prey ID|Peptide count|

Then run:

CompPASS(datInput)

CompPASSplus(datInput)

2. HGScore

For *datInput*, we need more column 'lenPrey', while 'idBait' is not necessary:

|idRun|idPrey|countPrey|lenPrey|
|-----|:----:|:----:|:----:|:----:|
|Unique ID of one affinity run|Prey ID|Peptide count|Prey protein length|

Then run:

HG(datInput)
