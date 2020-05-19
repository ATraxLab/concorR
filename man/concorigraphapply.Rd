\name{concor_igraph_apply}
\alias{concor_igraph_apply}
\title{Find CONCOR partition and add to a list of \pkg{igraph} objects}
\description{
Takes a list of \pkg{igraph} objects, runs \code{concor} until the desired number of splits \code{nsplit}, and adds the final split to each object as a vertex attribute.
}
\usage{
concor_igraph_apply(igraph_list, nsplit = 1)
}
\arguments{
  \item{igraph_list}{The list of \pkg{igraph} objects for use in \code{concor}}
  \item{nsplit}{Number of times to split, defaults to 1}
}
\details{
This function is a wrapper for a common task: Running \code{\link{concor}} on one or more \pkg{igraph} objects and adding the resulting partition to each object as a vertex attribute. If multiple \pkg{igraph} objects are included in the input list, they should be multiple relations for the same nodes.
}
\value{
Returns a list of \pkg{igraph} objects, each with a vertex attribute \code{csplitn}, where \code{n} is \code{nsplit} that contains the block assignment obtained from \code{concor}.
}

\seealso{
\code{\link{concor}}, \code{\link{concor_make_igraph}}
}
\examples{
a <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(a) <- c("a", "b", "c", "d")
colnames(a) <- c("a", "b", "c", "d")
a.igraph <- igraph::graph_from_adjacency_matrix(a)

concor_igraph_apply(list(a.igraph))
}
