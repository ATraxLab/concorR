\name{make_reduced_from_partition}
\alias{make_reduced_from_partition}
\title{Output reduced adjacency matrices given a partitioning of the nodes into appropriate groups}
\description{
Takes an adjacency matrix and a partition of the nodes, and returns a
reduced adjacency matrix and the corresponding cutoff density/degree. The 
\pkg{sna} package must be installed, but does not need to be attached.
}
\usage{
make_reduced_from_partition(adj_mat, partition, stat = 'density')
}
\arguments{
  \item{adj_mat}{A adjacency matrix}
  \item{partition}{A vector that gives the desired partitioning of the nodes. }
  \item{stat}{The statistic that will be used to determine if edges should be included in a reduced network.  Default value is \code{density}}
}
\details{
A reduced network represents each identified position in the network as a single node.
Links (or self-links) exist if the density (or normalized degree) of ties from that block to a target block is greater than a threshold density (or degree). In the default implementation, the density of the whole network is used as the threshold for each block.  In the degree implementation, the normalized degree of the network is used as the threshhold.
}
\value{
\item{reduced_mat}{A reduced matrix}
\item{dens}{The cut-off density used (equal to the edge density of \code{adj_mat}).  Only for \code{stat="density"}}
\item{deg}{The cut-off normalized outdegree used (equal to the mean outdegree of \code{adj_mat}).  Only for \code{stat="degree"}.}
}
\references{
S. Wasserman and K. Faust, \emph{Social Network Analysis: Methods and Applications} (Cambridge University Press, 1994).
}

\seealso{
  \code{\link{concor}}, \code{\link{make_blk}}
}
\examples{
g1 <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(g1) <- c("a", "b", "c", "d")
colnames(g1) <- c("a", "b", "c", "d")

g2 <- matrix(c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0), ncol = 4)
rownames(g2) <- c("a", "b", "c", "d")
colnames(g2) <- c("a", "b", "c", "d")

make_reduced(list(g1, g2), nsplit = 1)
make_reduced(list(g1, g2), nsplit = 1, stat="degree")
}
