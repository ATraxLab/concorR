\name{make_reduced}
\alias{make_reduced}
\title{Run CONCOR and output reduced adjacency matrices}
\description{
Takes a list of adjacency matrices, partitions using \code{\link{concor}}, and returns a
list of reduced adjacency matrices and their corresponding cutoff densities. Each reduced
matrix corresponds to one input relation. The \pkg{sna} package must be installed, but
does not need to be attached.
}
\usage{
make_reduced(adj_list, nsplit = 1, stat = 'density')
}
\arguments{
  \item{adj_list}{A list of adjacency matrices, each representing a different relation
  for the same nodes.}

  \item{nsplit}{The total number of splits CONCOR will (attempt to) perform.}
  \item{stat}{The statistic that will be used to determine if edges should be included in a reduced network.  Default value is \code{density}}
}
\details{
A reduced network represents each identified position in the network as a single node.
Links (or self-links) exist if the density (or normalized degree) of ties from that block to a target block is greater than a threshold density (or degree). In the default implementation, the density of the whole network is used as the threshold for each block.  In the degree implementation, the normalized degree of the network is used as the threshhold.

In the list of input matrices \code{adj_list}, each should correspond to a different
relation for the same nodes. Each adjacency matrix is partitioned with the CONCOR
algorithm, continuing for \code{nsplit} divisions. After the threshold density is
applied, each entry in the reduced matrix has values of 0 or 1.
}
\value{
\item{reduced_mat}{A list of reduced matrices, one for each input matrix.}
\item{dens}{A vector of the cut-off densities used (equal to the edge density of each entry in \code{adj_list}).  Only for \code{stat="density"}}
\item{deg}{A vector of the cut-off normalized degrees used (equal to the mean normalized degree of each entry in \code{adj_list}).  Only for \code{stat="degree"}.}
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
