\name{concor_plot}
\alias{concor_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{concor.plot
%%  ~~function to do ... ~~
}
\description{
Plots a network of interest with a concor output as the vertex color. The \code{igraph} object used should be created using \code{\link{make_igraph}} or \code{\link{concor_igraph_apply}}.
}
\usage{
concor.plot(iobject, p=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iobjet}{The igraph object with \code{concor} splitting as a vertex attribute. Should be created using \code{\link{make_igraph}} or \code{\link{concor_igraph_apply}}}
  \item{p}{The \code{concor} splitting of interest, \code{p} must correspond to a split number added to the igraph object}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Plots the network with the concor blocking as color.
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
\code{\link{concor}}, \code{\link{make_igraph}}, \code{\link{concor_igraph_apply}}
}
\examples{
a <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(a) <- c("a", "b", "c", "d")
colnames(a) <- c("a", "b", "c", "d")

i_out <- make_igraph(list(a))
concor_plot(i_out[[1]], 1)
}

