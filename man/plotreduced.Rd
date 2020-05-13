\name{plot_reduced}
\alias{plot_reduced}
\title{Plot a reduced network}
\description{
Plot a reduced network using \pkg{igraph}, with nodes colored by block number.
}
\usage{
plot_reduced(blk)
}
\arguments{
  \item{blk}{The igraph object representing the reduced network}
}
\details{
Plots the reduced network, where each node represents all the nodes assigned to that position by CONCOR. Node colors on the reduced matrix plot are assigned by position, so if CONCOR output is also used for vertex color on a sociogram, the node colors will align between the plots.
}
\value{
Returns \code{NULL}, invisibly.
}

\seealso{
\code{\link{make_reduced}}, \code{\link{make_reduced_igraph}}
}
\examples{
library(igraph)
g1 <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(g1) <- c("a", "b", "c", "d")
colnames(g1) <- c("a", "b", "c", "d")

#make the reduced matrix
r_mat <- make_reduced(list(g1), nsplit = 1)
r_igraph <- make_reduced_igraph(r_mat$reduced_mat[[1]])

plot_reduced(r_igraph)
}
