\name{make_reduced}
\alias{make_reduced}
\title{Run CONCOR and output reduced adjacency matrices}
\description{
Takes a list of adjacency matrices, partitions using \code{concor}, and returns a list of reduced adjacency matrices and their corresponding cutoff densities. Each reduced matrix corresponds to one input relation. Needs the \pkg{sna} package installed to work (does not need to be attached).
}
\usage{
make_reduced(adj_list, nsplit = 1)
}
\arguments{
  \item{adj_list}{List of adjacency matrices for use in \code{concor}}
  \item{nsplit}{The number of times \code{concor} splits each matrix}
}
\details{
A reduced network shows each identified position in the network as a single node. Links (or self-links) exist if the density of ties in that block is greater than a threshold density. In this implementation, the density of the whole network is used as the threshold for each block.

In the list of input matrices \code{adj_list}, each should correspond to a different relation for the same nodes. Each matrix will be partitioned with the CONCOR algorithm, continuing through split \code{nsplit}. After the threshold density is applied, each entry in a reduced matrix has values of 0 or 1.
}
\value{
\item{reduced_mat}{A list of reduced matrices, one for each input matrix.}
\item{dens}{A vector of the cut-off densities used (equal to the edge density of each entry in \code{adj_list}).}
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
}
