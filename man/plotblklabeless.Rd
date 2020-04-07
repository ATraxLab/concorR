\name{plot_blk_labeless}
\alias{plot_blk_labeless}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot_blk_labeless
%%  ~~function to do ... ~~
}
\description{
Runs \code{plot_blk} but gets rid of pesky labels.
}
\usage{
plot_blk_labeless(bm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bm}{An object of class \code{blockmodel}
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Plots a square blockmodel with no labels.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[sna]{plot.blockmodel}}, \code{\link{plot_blk}}
}
\examples{
g1 <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(g1) <- c("a", "b", "c", "d")
colnames(g1) <- c("a", "b", "c", "d")

gl <- list(g1)
bm <- make_blk(gl, 1)[[1]]
plot_plk_labeless(bm)
}

