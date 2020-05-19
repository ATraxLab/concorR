\name{concor}
\alias{concor}
\title{Find CONCOR partition for a graph}
\description{
Use the CONCOR (CONvergence of iterated CORrelations) algorithm to identify roles within social network data.
}
\usage{
concor(m_list, nsplit = 1, self_ties = FALSE, cutoff = .9999999, max_iter = 50)
}
\arguments{
  \item{m_list}{A list of adjacency matrices. Matrices must be square, of the same dimensions, and have row/column names (node names). Each matrix represents a different relationship of the same nodes. If there is only one relationship of interest, \code{m_list} is a list of that single matrix.}
  \item{nsplit}{The number of times the input matrices will be partitioned. Each split divides a partition in two, resulting in 2^\code{nsplit} roles identified.}
  \item{self_ties}{A Boolean representing whether self-ties can be present in the data, \code{TRUE} allows for self-ties, \code{FALSE} does not.}
  \item{cutoff}{The cutoff point of the iterated correlations. Both this and \code{max_iter} can be lowered slightly to improve speed.}
  \item{max_iter}{The maximum number of times the correlation will be run for a split. Both this and \code{cutoff} can be lowered slightly to improve speed.}
}
\details{
This version works for multiple relationships, assuming they all are for the same data (same size of input matrices), and can be used with isolates present. It requires further testing on weighted networks (but appears to successfully split such networks). It will attempt to split the network \code{nsplit} times, causing there to be 2^nsplit partitions, plus one for isolated nodes (if they exist), unless a structurally equivalent node group or singular node group is present. If the algorithm attempts to split such a node group the function will ignore that group and continue to split all other blocks until specified. If a higher number of splits, \code{nsplit}, are requested than are possible to apply to the specified data (due to structurally equivalent node groups being present or all blocks being singular nodes) the code will warn the user that split \code{nsplit} was the same as split \code{i}, the final possible split, and stop.
}
\value{
A data frame with two columns: \code{block} is the block or role identified by concor, and \code{vertex} is the node names.
}
\references{
R. L. Breiger, S. A. Boorman, P. Arabie, An algorithm for clustering relational data with applications to social network analysis and comparison with multidimensional scaling. \emph{J. of Mathematical Psychology}. \strong{12}, 328 (1975). \doi{10.1016/0022-2496(75)90028-0}

S. Wasserman and K. Faust, \emph{Social Network Analysis: Methods and Applications} (Cambridge University Press, 1994).
}

\examples{
a <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(a) <- c("a", "b", "c", "d")
colnames(a) <- c("a", "b", "c", "d")

concor(list(a))
}

