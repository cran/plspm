\name{ropes}
\alias{ropes}
\docType{data}
\title{Climbing Ropes dataset}
\description{
  This dataset gives the measurements of 101 climbing ropes available in the market by spring 2011. The data was collected from the brands websites.
}
\usage{data(ropes)}
\format{
  A data frame with 101 observations on the following 7 variables.
  \describe{
    \item{\code{brand}}{a factor with the brand names}
    \item{\code{diameter}}{the diameter of the ropes measured in mm}
    \item{\code{weight}}{the weight measured in grams per meter}
    \item{\code{falls}}{the number of UIAA falls}
    \item{\code{imp.force}}{the impact force value}
    \item{\code{stat.elong}}{the static elongation value}
    \item{\code{dyn.elong}}{the dynamic elongation value}
  }
}
\references{
  Sanchez, G. (2012) PLS Modeling (In Press)
}
\examples{
   data(ropes)
   summary(ropes)
}
\keyword{datasets}
