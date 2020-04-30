\name{concor}
\alias{concor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find CONCOR partition for a graph
%%  ~~function to do ... ~~
}
\description{
Use the CONCOR (CONvergence of iterated CORrelations) algorithm to identify roles within social network data.
}
\usage{
concor(m0, p = 1, self_ties = FALSE, cutoff = .9999999, max_iter = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m0}{A list of adjacency matrices. Matrices must be square, of the same dimensions, and have row/column names (node names). Each matrix represents a different relationship of the same nodes. If there is only one relationship of interest, \code{m0} is a list of that single matrix.}
  \item{p}{The number of times the input matrices will be partitioned. Each split divides a partition in two, resulting in 2^p roles identified.}
  \item{self_ties}{A Boolean representing whether self-ties can be present in the data, \code{TRUE} allows for self-ties, \code{FALSE} does not.}
  \item{cutoff}{The cutoff point of the iterated correlations. Both this and \code{max_iter} can be lowered slightly to improve speed.}
  \item{max_iter}{The maximum number of times the correlation will be run for a split. Both this and \code{cutoff} can be lowered slightly to improve speed.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This version works for multiple relationships, assuming they all are for the same data (same size of input matrices), and can be used with isolates present. It requires further testing on weighted networks (but there is no evidence of it failing on them yet). It will attempt to split the network \code{p} times, causing there to be 2^p partitions, plus one for isolated nodes (if they exist). If the algorithm attempts to split a structurally equivalent node group, or a single-node group, the function will fail. The current implementation cannot ignore such a group and continue to split subsequent groups, so in these cases, the number of desired splits (\code{p}) must be lowered.
}
\value{
A data frame with two columns: \code{block} is the block or role identified by concor, and \code{vertex} is the node names.
}
\references{
%% ~put references to the literature/web site here ~
R. L. Breiger, S. A. Boorman, P. Arabie, An algorithm for clustering relational data with applications to social network analysis and comparison with multidimensional scaling. \emph{J. of Mathematical Psychology}. \strong{12}, 328 (1975). \doi{10.1016/0022-2496(75)90028-0}

S. Wasserman and K. Faust, \emph{Social Network Analysis: Methods and Applications} (Cambridge University Press, 1994).
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
a <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(a) <- c("a", "b", "c", "d")
colnames(a) <- c("a", "b", "c", "d")

concor(list(a))
}

