\name{russett}
\alias{russett}
\docType{data}
\title{Spanish football dataset}
\description{
  Data set from Russett (1964) about agricultural inequality, industrial development and political instability.
}
\usage{russett}
\format{
  A data frame with 47 observations on the following 11 variables. The variables may be used to construct three latent concepts: 1) \code{AGRIN}=Agricultural Inequality, 2) \code{INDEV}=Industrial Development, 3) \code{POLINS}=Political Instability. 

  \code{AGRIN}: Recommended variables gini, farm, rent. \cr
  \code{INDEV}: Recommended variables gnpr, labo. \cr
  \code{POLINS}: Recommended variables inst, ecks, death, demostab, demoinst, dictator. \cr
  \describe{
    \item{\code{gini}}{Inequality of land distribution}
    \item{\code{farm}}{Percentage of farmers that own half of the land}
    \item{\code{rent}}{Percentage of farmers that rent all their land}
    \item{\code{gnpr}}{Gross national product per capita}
    \item{\code{labo}}{Percentage of labor force employed in agriculture}
    \item{\code{inst}}{Instability of executive (45-61)}
    \item{\code{ecks}}{Number of violent internal war incidents (41-61)}
    \item{\code{death}}{Number of people killed as a result of civic group violence (50-62)}
    \item{\code{demostab}}{Political regime: stable democracy}
    \item{\code{demoinst}}{Political regiem: unstable democracy}
    \item{\code{dictator}}{Political regime: dictatorship}
  }
}
\references{
  Russett B.M. (1964) Inequality and Instability: The Relation of Land Tenure to Politics. \emph{World Politics} \bold{16:3}, pp. 442-454.
}
\examples{
  data(russett)
  russett
}
\keyword{datasets}
