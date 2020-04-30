\name{concor_igraph_apply}
\alias{concor_igraph_apply}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find CONCOR partition and add to a list of igraph objects
%%  ~~function to do ... ~~
}
\description{
Takes a list of graph objects, runs \code{concor} until the desired number of splits \code{nsplit}, and adds the final split to each \code{igraph} object as a vertex attribute.
}
\usage{
concor_igraph_apply(igraph_list, nsplit = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{igraph.list}{The list of igraph objects for use in \code{concor}}
  \item{nsplit}{Number of times to split, defaults to 1}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function is a wrapper for a common task: Running \code{\link{concor}} on one or more igraph objects and adding the resulting partition to each object as a vertex attribute. If multiple igraph objects are included in the input list, they should be multiple relations for the same nodes.
}
\value{
Returns a list of igraph objects, each with a vertex attribute \code{csplitn}, where \code{n} is \code{nsplit} that contains the block assignment obtained from \code{concor}.
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
