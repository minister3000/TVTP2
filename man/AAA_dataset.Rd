\name{AAA_dataset}
\alias{AAA_dataset}
\docType{data}
\title{
Final set of appropriately modified observations for the dependent and explanatory variables for estimating the \code{AAA} spread system
}
\description{
Final dataset of observations after the application of lag and log operators and first differencing where appropriate. This dataset can be used directly in the log-likelihood function for estimating \code{AAA} credit spread dynamics. The dataset contains explanatory variables impacting credit spreads directly, and explanatory variables governing the regime transition process. The final column in the dataset is BBB spread changes month-over-month.
}
\usage{data("AAA_dataset")}
\format{
  The format is:
 data.frame "AAA_dataset"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(AAA_dataset)
## str(AAA_dataset) ; plot(AAA_dataset) ...
}
\keyword{datasets}
