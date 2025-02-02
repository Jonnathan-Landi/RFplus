\name{BD_Insitu}
\alias{BD_Insitu}
\docType{data}
\title{
Precipitation Station Measurement Dataset
}
\description{
This dataset contains daily measurements from several precipitation stations. The first column represents the measurement date, and the following columns correspond to the measurements from each station on that date. The station columns are labeled with unique identifiers for each station, and the number of stations may vary depending on the dataset configuration.
}
\usage{data("BD_Insitu")}
\format{
  A `data.table` object with station measurements. The dataset includes the following columns:
  \describe{
    \item{\code{Date}}{The measurement date (type \code{Date}).}
    \item{\code{Station_ID_1, Station_ID_2, ...}}{Measurements from the stations (numeric values). Each column after \code{Date} represents the measurements of a precipitation station for the corresponding date. The columns are labeled with unique identifiers (e.g., \code{Station_ID_1}, \code{Station_ID_2}, etc.) for each station, and the number of stations (columns) may vary.}
  }
}
\details{
The data represents daily measurements taken from several precipitation stations. The first column contains the measurement dates, and the following columns represent the measurements of each station on those dates. The number of stations may vary depending on the dataset, and each station is uniquely identified by its column name (e.g., \code{Station_ID_1}, \code{Station_ID_2}, etc.).
}
\source{
The data was generated for use in the bias correction model for satellite products, RFplus.
}
\examples{
data(BD_Insitu)
## You can use str(BD_Insitu) to get a description of the structure
## or view some of the first rows using head(BD_Insitu)
}
\keyword{datasets}
