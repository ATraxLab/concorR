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
blk.apply(iobject, split, v="cat")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iobject}{The igraph object representing same data used in \cose{concor}}
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
a=matrix(c(0,0,0,1,1,0,1,0,1,1,0,1,1,0,1,0), ncol=4)
rownames(a)=c("a","b","c","d")
colnames(a)=c("a","b","c","d")
a.igraph=graph_from_adjacency_matrix(a)

con.out=concor(list(a))

blk.apply(a.igraph, con.out, "split_1")
}
