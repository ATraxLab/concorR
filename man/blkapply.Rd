\name{blk_apply}
\alias{blk_apply}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{blk.apply
%%  ~~function to do ... ~~
}
\description{
Adds the split taken from the output of \code{concor} to its \code{igraph} object as a vertex attribute.
}
\usage{
blk_apply(iobject, split, v = "cat")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iobject}{The igraph object representing same data used in \code{concor}}
  \item{split}{The output form \code{concor} for a specific split}
  \item{v}{What the vertex attribute will be called, defaults to "cat"}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the igraph object with the new vertex attribute v.
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
\code{\link{concor}}
}
\examples{
library(igraph)
a <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(a) <- c("a", "b", "c", "d")
colnames(a) <- c("a", "b", "c", "d")
a_igraph <- igraph::graph_from_adjacency_matrix(a)

con_out <- concor(list(a))

blk_apply(a_igraph, con_out, "split_1")
}
