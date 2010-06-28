\name{local.models}
\alias{local.models}
\title{Calculates PLS-PM for global and local models}
\description{
  Calculates PLS-PM for global and local models from a given partition.
}
\usage{
  local.models(pls, y, scheme=NULL, scaled=NULL, boot.val=FALSE, br=NULL)
}
\arguments{
  \item{pls}{An object of class \code{"plspm"}}
  \item{y}{One object of the following classes: \code{"rebus"}, \code{"integer"}, or \code{"factor"}, that provides the class partitions. }
  \item{scheme}{Possible values are \code{"centroid"}, \code{"factor"} and \code{"path"}.}
  \item{scaled}{A logical value indicating whether scaling data is performed.}
  \item{boot.val}{A logical value indicating whether bootstrap validation is performed (\code{FALSE} by default).}
  \item{br}{An integer indicating the number bootstrap resamples. Used only when \code{boot.val=TRUE}.}

}
\details{
  The function \code{local.models} calculates PLS-PM for the global model (i.e. over all observations) as well as PLS-PM for local models (i.e. observations of different partitions). \cr
  When \code{y} is an object of class \code{"rebus"}, the function \code{local.models} is applied to the classes obtained from the REBUS algorithm. \cr
  When \code{y} is an \code{integer} vector or a \code{factor}, the values or levels are assumed to represent the group to which each observation belongs. In this case, the function \code{local.models} calculates PLS-PM for the global model, as well as PLS-PM for each group (local models). \cr
  If \code{scheme=NULL}, then the original scheme from the object \code{pls} will be taken. \cr
  If \code{scaled=NULL}, then the original scaled from the object \code{pls} will be taken. \cr
  When bootstrap validation is performed, the default number of re-samples is 100. However, \code{br} can be specified in a range from 100 to 1000.
}
\value{
  An object of class \code{"local.models"}, basically a list of length \code{k+1}, where \code{k} is the number of classes. The list contains the following elements:
  \item{glob.model}{PLS-PM of the global model}
  \item{loc.model.1}{PLS-PM of segment (class) 1}
  \item{loc.model.2}{PLS-PM of segment (class) 2}
  \item{loc.model.k}{PLS-PM of segment (class) k}
  Each element of the list is an object of class \code{"plspm"}. Thus, in order to examine the
  results for each local model, it is necessary to use the \code{summary} function. See \code{examples} below.

  Note that if \code{scheme} and/or \code{scaled} differ from the original arguments of \code{pls}, the results obtained in \code{local.models} will be different from those obtained in \code{pls} as well as in \code{y} (if it is of class \code{"rebus"}).
}
\author{Laura Trinchera, Gaston Sanchez}
\seealso{\code{\link{rebus.pls}}}
\examples{
  \dontrun{
  ## example of rebus analysis
  data(sim.data)
  ## First compute GLOBAL model
  sim.mat <- matrix(c(0,0,0,0,0,0,1,1,0),3,3,byrow=TRUE)
  dimnames(sim.mat) <- list(c("Price","Quality","Satisfaction"),
                            c("Price","Quality","Satisfaction"))
  sim.sets <- list(c(1,2,3,4,5),c(6,7,8,9,10),c(11,12,13)) 
  sim.mod <- c("A","A","A")  ## reflective indicators
  sim.global <- plspm(sim.data, inner.mat=sim.mat, 
                      sets=sim.sets, modes=sim.mod)
  sim.global
  ## Then compute cluster on residual from global model
  sim.res.clus <- res.clus(sim.global)
  ## To conclude run iteration algorithm
  rebus.sim <- it.reb(sim.global, sim.res.clus, nk=2, 
                     stop.crit = 0.005, iter.max = 100 )
  ## Computation of local models 
  local.rebus <- local.models(sim.global, rebus.sim)
  ## Display plspm summary for first local model 
  summary(local.rebus$loc.model.1)
  }
}
