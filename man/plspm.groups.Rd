\name{plspm.groups}
\alias{plspm.groups}
\title{Group Comparison in PLS-PM}
\description{
  This function performs a group comparison test for comparing path coefficients between two groups 
}
\usage{
  plspm.groups(x, pls, g, method = "bootstrap", reps = NULL)
}
\arguments{
  \item{x}{A numeric matrix or data frame containing the manifest variables.}
  \item{pls}{An object of class \code{"plspm"} returned by the function \code{\link{plspm}}}
  \item{g}{A factor with 2 levels indicating the groups to be compared}
  \item{method}{The method to be used in the test. Possible values are \code{"bootstrap"} or \code{"permutation"} }
  \item{reps}{An integer indicating the number of either bootstrap resamples or number of permutations. If \code{NULL} then \code{reps}=100}
}
\details{
  The function \code{plspm.groups} performs a two groups comparison test in PLS-PM for comparing path coefficients between two groups. For the moment, only two methods are available: 1) bootstrap, and 2) permutation. \cr
  The bootstrap test is an adapted t-test based on bootstrap standard errors. \cr
  The permutation test is a randomization test which provides a non-parametric option. \cr

  The null and alternative hypotheses to be tested are: \cr
  H0: path coefficients are not significantly different \cr
  H1: path coefficients are significantly different \cr
  
  In order to use the function \code{plspm.groups} the argument \code{x} must be the same argument used in the function \code{\link{plspm}} to calculate the object \code{pls}. \cr  
}
\value{
  An object of class \code{"plspm.groups"}, basically a list with the following elements:
  \item{test}{Table with the results of the applied test. Includes: path coefficients of the global model, path coeffs of group1, path coeffs of group2, (absolute) difference of path coeffs between groups, and the test results with the p-value.}
  \item{global}{List with results of the inner model for the global model}
  \item{group1}{List with results of the inner model for group1}
  \item{group2}{List with results of the inner model for group2}
}
\references{ 
  Chin, W.W. (2003) A permutation procedure for multi-group comparison of PLS models. In: Vilares, M., Tenenhaus, M., Coelho, P., Esposito Vinzi, V., Morineau, A. (Eds.) \emph{PLS and Related Methods - Proceedings of the International Symposium PLS03.} Decisia, pp. 33-43. 

  Chin, W.W. (2000) Frequently Asked Questions, Partial Least Squares PLS-Graph. Available from: \cr
  \url{http://disc-nt.cba.uh.edu/chin/plsfaq/multigroup.htm}
}
\author{Gaston Sanchez}
\seealso{\code{\link{print.plspm.groups}},\code{\link{plot.plspm.groups}}}
\examples{
  \dontrun{
  ## example with customer satisfaction analysis
  ## group comparison based on the segmentation variable "gender"
  data(satisfaction)
  IMAG <- c(0,0,0,0,0,0)
  EXPE <- c(1,0,0,0,0,0)
  QUAL <- c(1,1,0,0,0,0)
  VAL  <- c(0,1,1,0,0,0)
  SAT  <- c(1,1,1,1,0,0) 
  LOY  <- c(1,0,0,0,1,0)
  sat.mat <- rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
  sat.sets <- list(1:5,6:10,11:15,16:19,20:23,24:27)
  sat.mod <- rep("A",6)   ## reflective indicators
  pls <- plspm(satisfaction, sat.mat, sat.sets, sat.mod, scheme="factor", 
                scaled=FALSE)
  ## permutation test with 100 permutations
  res.group <- plspm.groups(satisfaction, pls, satisfaction$gender, 
               method="permutation", reps=100)
  res.group
  plot(res.group)
  }
}
