\name{outerplot}
\alias{outerplot}
\title{Plot outer model}
\usage{
  outerplot(x, what = "loadings", colpos = "#6890c4BB",
    colneg = "#f9675dBB", box.prop = 0.55, box.size = 0.08,
    box.cex = 1, box.col = "gray95", lcol = "gray95",
    box.lwd = 2, txt.col = "gray40", shadow.size = 0,
    curve = 0, lwd = 2, arr.pos = 0.5, arr.width = 0.15,
    cex.txt = 0.9, ...)
}
\arguments{
  \item{x}{An object of class \code{"plspm"}.}

  \item{what}{What to plot: \code{"loadings"} or
  \code{"weights"}.}

  \item{colpos}{Color of arrows for positive path
  coefficients.}

  \item{colneg}{Color of arrows for negative path
  coefficients.}

  \item{box.prop}{Length/width ratio of ellipses and
  rectangles.}

  \item{box.size}{Size of ellipses and rectangles.}

  \item{box.cex}{Relative size of text in ellipses and
  rectangles.}

  \item{box.col}{fill color of ellipses and rectangles.}

  \item{lcol}{border color of ellipses and rectangles.}

  \item{box.lwd}{line width of the box.}

  \item{txt.col}{color of text in ellipses and rectangles.}

  \item{shadow.size}{Relative size of shadow of label box.}

  \item{curve}{arrow curvature.}

  \item{lwd}{line width of arrow.}

  \item{arr.pos}{Relative position of arrowheads on
  arrows.}

  \item{arr.width}{arrow width.}

  \item{cex.txt}{Relative size of text on arrows.}

  \item{\dots}{Further arguments passed on to
  \code{\link{plotmat}}.}
}
\description{
  Plot either outer weights or loadings in the outer model
  for objects of class \code{"plspm"}
}
\note{
  \code{outerplot} uses the function \code{\link{plotmat}}
  of package \code{diagram}. \cr
  \url{https://cran.r-project.org/package=diagram/vignettes/diagram.pdf}
}
\seealso{
  \code{\link{innerplot}}, \code{\link{plot.plspm}},
  \code{\link{plspm}}
}

