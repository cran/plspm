\name{innerplot}
\alias{innerplot}
\title{Plot inner model}
\usage{
  innerplot(x, colpos = "#6890c4BB", colneg = "#f9675dBB",
    box.prop = 0.55, box.size = 0.08, box.cex = 1,
    box.col = "gray95", lcol = "gray95", box.lwd = 2,
    txt.col = "gray50", shadow.size = 0, curve = 0,
    lwd = 3, arr.pos = 0.5, arr.width = 0.2, arr.lwd = 3,
    cex.txt = 0.9, show.values = FALSE, ...)
}
\arguments{
  \item{x}{Either a matrix defining an inner model or an
  object of class \code{"plspm"}.}

  \item{colpos}{Color of arrows for positive path
  coefficients.}

  \item{colneg}{Color of arrows for negative path
  coefficients.}

  \item{box.prop}{Length/width ratio of ellipses.}

  \item{box.size}{Size of ellipses.}

  \item{box.cex}{Relative size of text in ellipses.}

  \item{box.col}{fill color of ellipses,}

  \item{lcol}{border color of ellipses.}

  \item{box.lwd}{line width of the box.}

  \item{txt.col}{color of text in ellipses.}

  \item{shadow.size}{Relative size of shadow of label box.}

  \item{curve}{arrow curvature.}

  \item{lwd}{line width of arrow.}

  \item{arr.pos}{Relative position of arrowheads on
  arrows.}

  \item{arr.width}{arrow width.}

  \item{arr.lwd}{line width of arrow, connecting two
  different points, (one value, or a matrix with same
  dimensions as \code{x}).}

  \item{cex.txt}{Relative size of text on arrows.}

  \item{show.values}{should values be shown when \code{x}
  is a matrix.}

  \item{\dots}{Further arguments passed on to
  \code{\link{plotmat}}.}
}
\description{
  Plot the inner (structural) model for objects of class
  \code{"plspm"}, as well as path matrices
}
\note{
  \code{innerplot} uses the function \code{\link{plotmat}}
  in package \code{diagram}. \cr
  \url{https://cran.r-project.org/package=diagram/vignettes/diagram.pdf}
}
\seealso{
  \code{\link{plot.plspm}}, \code{\link{outerplot}}
}

