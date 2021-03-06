% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CompPASSplus.R
\name{CompPASSplus}
\alias{CompPASSplus}
\title{CompPASSplus}
\usage{
CompPASSplus(datInput)
}
\arguments{
\item{datInput}{A dataframe with column names: idRun, idBait, idPrey, countPrey.
Each row represent one unique protein captured in one pull-down experiment}
}
\value{
A data frame consists of unique bait-prey pairs with Z-score, S-score,
D-score and WD-score indicating interacting probabilities.
}
\description{
CompPASSplus
Implemented a naive Bayes classifier that learns to distinguish true
 interacting proteins from non-specific background and false positive
 identifications on the basis of CompPASS scores.
 The source code for this function was based on the source code. \url{https://github.com/dnusinow/cRomppass}
}
\author{
Qingzhou Zhang, \email{zqzneptune@hotmail.com}
}
