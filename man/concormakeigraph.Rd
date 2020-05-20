\name{concor_make_igraph}
\alias{concor_make_igraph}
\title{Find CONCOR partition and create a list of \code{igraph} objects}
\description{
Run \code{concor} on a list of adjacency matrices and create a list of \code{igraph}
objects with a specific CONCOR split added as a vertex attribute.
}
\usage{
concor_make_igraph(adj_list, nsplit = 1)
}
\arguments{
  \item{adj_list}{A list of adjacency matrices, each representing a different relation
  for the same nodes.}

  \item{nsplit}{The total number of splits CONCOR will (attempt to) perform.}
}
\details{
This function is a wrapper for a common task: Running \code{\link{concor}} on one or more
adjacency matrices, then saving the networks as \code{igraph} objects with the CONCOR
partition as a vertex attribute. If multiple adjacency matrices are included in the input
list, they should be multiple relations for the same nodes.
}
\value{
Returns a list of \code{igraph} objects with a vertex attribute \code{csplitn}, where
\code{n} is \code{nsplit} that contains the block assignment obtained from \code{concor}.
}

\seealso{
\code{\link{concor}}, \code{\link{concor_igraph_apply}}
}
\examples{
a <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(a) <- c("a", "b", "c", "d")
colnames(a) <- c("a", "b", "c", "d")

concor_make_igraph(list(a))
}
