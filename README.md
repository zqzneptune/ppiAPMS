# ppiAPMS[Deprecated]
This is an R package that implements statistical modeling approaches that assign confidence scores to protein-protein interaction data generated using affinity purification–mass spectrometry (AP-MS) data.Please use SMAD(http://bioconductor.org/packages/release/bioc/html/SMAD.html)

## Installation

The development version can be installed through github:
```{r}
 devtools::install_github(repo="zqzneptune/ppiAPMS")
 library(ppiAPMS)
```
## Quick start

1. CompPASS and CompPASS-Plus

Summarize your AP-MS data from proteome database search into the dataframe *datInput* with the following format:

|idRun|idBait|idPrey|countPrey|
|-----|:----:|:----:|:-------:|
|Unique ID of one AP-MS run|Bait ID|Prey ID|Prey peptide count|

Then run:

```{r}
CompPASS(datInput)
```

```{r}
CompPASSplus(datInput)
```

2. HGScore

For *datInput*, we need more column 'lenPrey', while 'idBait' is not necessary:

|idRun|idPrey|countPrey|lenPrey|
|-----|:----:|:----:|:----:|
|Unique ID of one AP-MS run|Prey ID|Prey peptide count|Prey protein length|

Then run:

```{r}
HG(datInput)
```
## License

MIT @ Qingzhou Zhang
