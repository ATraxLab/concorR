\name{concor_igraph_apply}
\alias{concor_igraph_apply}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{concor_igraph_apply
%%  ~~function to do ... ~~
}
\description{
Takes a list of graph objects and runs \code{concor} on the desired split \code{nsplit} and adds the outputted split to each \code{igraph} object as an vertex attribute.
}
\usage{
concor.igraph.apply(igraph.list, nsplit=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{igraph.list}{The list of igraph objects for use in \code{concor}}
  \item{nsplit}{The split number desired, defaults to 1}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list of igraph objects with a vertex attribute \code{csplit_n} where \code{n} is \code{nsplit} that contain the block information obtained from \code{concor}.
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
\code{\link{concor}}, \code{\link{blk_apply}}, \code{\link{concor_make_igraph}}
}
\examples{
a <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(a) <- c("a", "b", "c", "d")
colnames(a) <- c("a", "b", "c", "d")
a.igraph <- igraph::graph_from_adjacency_matrix(a)

concor_igraph_apply(list(a.igraph))
}
