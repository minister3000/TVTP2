\name{RawEconomicData}
\alias{RawEconomicData}
\docType{data}
\title{
(Raw) Economic data time series used in dissertation effort
}
\description{
Monthly economic data series from January 1920 to December 2018 as provided by the ultimate data sources described in the 'Data' section of the dissertation. The only exception to this is the calculation of the monthly default data and inflation rates. Annual observations provided by the ultimate data source have been divided by twelve to obtain average monthly default estimates. Annual inflation rates are calculated as the as the 12th difference in the natural logs of the original 'Consumer Price Index' index time series. For additional information on data sources and data assembly refer to the 'Data' section in the dissertation.
}
\usage{data("RawEconomicData")}
\format{
  The format is:
 data.frame "RawEconomicData"
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
data(RawEconomicData)
## str(RawEconomicData) ; plot(RawEconomicData) ...
}
\keyword{datasets}
