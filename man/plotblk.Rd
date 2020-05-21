\name{plot_blk}
\alias{plot_blk}
\title{Plot a blockmodel}
\description{
Displays a plot of a blockmodel. Based on \code{\link[sna]{plot.blockmodel}} (\pkg{sna}
must be installed), but reformats the plot to be square and removes the mandatory title.
}
\usage{
plot_blk(x, labels = FALSE, ...)
}
\arguments{
  \item{x}{An object of class \code{blockmodel}}
  \item{labels}{If \code{TRUE}, vertex ids are displayed as row/column/diagonal labels.
  If \code{FALSE}, no node labels are shown.}
  \item{...}{Further arguments passed to or from other methods}
}
\details{
This is a modification of the \code{\link[sna]{plot.blockmodel}} function. The original
displays vertex ids as row, column, and diagonal labels, which can be unreadable for
larger networks. \code{\link[sna]{plot.blockmodel}} also adds a title of the form
"Relation - 1", which this version omits.
}
\value{
Returns \code{NULL}, invisibly.
}
\references{
Carter T. Butts (2019). sna: Tools for Social Network Analysis. R package version 2.5.
\url{https://CRAN.R-project.org/package=sna}
}
\author{
Carter T. Butts (buttsc@uci.edu)

Modified by Tyme Suda
}

\seealso{
\code{\link[sna]{plot.blockmodel}}
}
\examples{
g1 <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(g1) <- c("a", "b", "c", "d")
colnames(g1) <- c("a", "b", "c", "d")

gl <- list(g1)
bm <- make_blk(gl, 1)[[1]]
plot_blk(bm)
}

