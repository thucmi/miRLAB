% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/miRLAB.R
\name{ValidateAll}
\alias{ValidateAll}
\title{Validate the targets of all miRNA using both experimentally confirmed and transfection data}
\usage{
ValidateAll(CEmatrix, topk, groundtruth, LFC, downreg = TRUE)
}
\arguments{
\item{CEmatrix}{the matrix of correlation/causal effects/scores with columns are miRNAs and rows are mRNAs}

\item{topk}{the number of targets of each miRNA that are being validated.}

\item{groundtruth}{the csv file containing the ground truth.}

\item{LFC}{the log fold change threshold for the transfection data. The targets that have the absolute value of log fold change greater than the LFC
will be regarded as the confirmed targets.}

\item{downreg}{if TRUE the negative correlation/causal effect/score values will be ranked on the top of the ranking. This is to favour the down regulations.}
}
\value{
a list of matrices that contains the confirmed interactions by both provided ground truth and built-in transfection data.
}
\description{
Given the predicted target of all miRNA, the function returns a list of targets of each miRNA that are  confirmed
based on the experimentally validated interactions or curated transfection data. Users need to download the file logFC.imputed.rda from nugget.unisa.edu.au/Thuc/miRLAB/ and place it in the working directory (this file is obtained from the TargetScoreData package)
}
\examples{
print("ps=Pearson(dataset, cause=1:3, effect=4:18)")
print("results=ValidateAll(ps, 10, groundtruth, LFC=0.5, downreg=TRUE)")
}

