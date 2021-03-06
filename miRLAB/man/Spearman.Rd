% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/miRLAB.R
\name{Spearman}
\alias{Spearman}
\title{miRNA target prediction with the Spearman correlation coefficient method}
\usage{
Spearman(datacsv, cause, effect, targetbinding = NA)
}
\arguments{
\item{datacsv}{the input dataset in csv format}

\item{cause}{the column range that specifies the causes (miRNAs), e.g. 1:35}

\item{effect}{the column range that specifies the effects (mRNAs), e.g. 36:2000}

\item{targetbinding}{the putative target, e.g. "TargetScan.csv". If targetbinding is not specified, only expression data is used.
If targetbinding is specified, the prediction results using expression data with be intersected with the interactions in the target binding file.}
}
\value{
A  matrix that includes the Spearman correlation coefficients. Columns are miRNAs, rows are mRNAs.
}
\description{
Calculate the Spearman correlation coefficient of each pair of miRNA-mRNA,and return a matrix of correlation coefficients with columns are miRNAs and rows are mRNAs.
}
\examples{
dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
results=Spearman(dataset, 1:3, 4:18)
}
\references{
Spearman, C. (1904) General intelligence, objectively determined and measured. Am. J. Psychol., 15, 201 - 92.
}

