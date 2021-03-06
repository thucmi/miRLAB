% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/miRLAB.R
\name{ReadExtResult}
\alias{ReadExtResult}
\title{Read results from other methods}
\usage{
ReadExtResult(datacsv, cause, effect, ExtCEcsv)
}
\arguments{
\item{datacsv}{the input dataset in csv format}

\item{cause}{the column range that specifies the causes (miRNAs), e.g. 1:35}

\item{effect}{the column range that specifies the effects (mRNAs), e.g. 36:2000}

\item{ExtCEcsv}{score matrix predicted by an external matrix with columns are miRNAs and rows are mRNAs.}
}
\value{
a matrix of scores predicted by an external matrix and ready for further validation and comparison tasks.
}
\description{
Read the results predicted by external methods (methods that are not in this package and may not be implemented in R). Consequently, we can compare the results
predicted by the external methods and results predicted by the methods in the miRLAB package.
}
\examples{
print("GenemiR=ReadExtResult(dataset, cause=1:3, effect=4:18, 'genemirresults.csv')")
}

