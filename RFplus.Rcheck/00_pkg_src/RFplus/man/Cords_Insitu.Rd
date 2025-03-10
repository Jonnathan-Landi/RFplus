\name{Cords_Insitu}
\alias{Cords_Insitu}
\docType{data}
\title{
Precipitation Station Coordinates Dataset
}
\description{
This dataset contains the coordinates (in UTM format) of several precipitation stations. Each station is uniquely identified by the \code{Cod} column, which corresponds to the station identifiers used in the \code{BD_Insitu} dataset. The coordinates of each station are provided in two columns: \code{X} for the Easting (longitude) and \code{Y} for the Northing (latitude).
}
\usage{data("Cords_Insitu")}
\format{
  A `data.table` object with station coordinates. The dataset includes the following columns:
  \describe{
    \item{\code{Cod}}{The unique identifier for each station. This should correspond to the station columns in the \code{BD_Insitu} dataset.}
    \item{\code{X}}{The Easting (X-coordinate) of the station in UTM format (numeric).}
    \item{\code{Y}}{The Northing (Y-coordinate) of the station in UTM format (numeric).}
  }
}
\details{
The data represents the geographic coordinates of precipitation stations used in the analysis. The first column, \code{Cod}, contains the unique identifiers of the stations, which should match the column names in the \code{BD_Insitu} dataset. The subsequent columns, \code{X} and \code{Y}, contain the UTM coordinates for each station, representing the station's location on the Earth's surface.
}
\source{
The data was generated for use in the bias correction model for satellite products, RFplus.
}
\examples{
data(Cords_Insitu)
## You can use str(Cords_Insitu) to get a description of the structure
## or view some of the first rows using head(Cords_Insitu)
}
\keyword{datasets}
