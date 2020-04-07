\name{plot_blk}
\alias{plot_blk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot.blk
%%  ~~function to do ... ~~
}
\description{
Displays a plot of a blockmodel. Is \code{\link[sna]{plot.blockmodel}} but reformats the shape to be square and removed the "relation - " form the titling of the output.
}
\usage{
plot.blk(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{blockmodel}}
  \item{...}{Further arguments passed to or from other methods}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Plots a blockmodel that is square.
}
\references{
Carter T. Butts (2019). sna: Tools for Social Network Analysis. R package version 2.5.
https://CRAN.R-project.org/package=sna
}
\author{
Carter T. Butts buttsc@uci.edu

edited by Tyme Suda
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[sna]{plot.blockmodel}}
}
\examples{
g1 <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(g1) <- c("a", "b", "c", "d")
colnames(g1) <- c("a", "b", "c", "d")

gl <- list(g1)
bm <- make_blk(gl, 1)[[1]]
plot_plk(bm)
}

