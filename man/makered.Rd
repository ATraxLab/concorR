\name{make_reduced}
\alias{make_reduced}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Run CONCOR and output reduced adjacency matrices
%%  ~~function to do ... ~~
}
\description{
Takes a list of adjacency matrices, partitions using \code{concor}, and returns a list of reduced adjacency matrices and their corresponding cutoff densities. Each reduced matrix corresponds to one input relation and can be weighted or unweighted. Needs the \pkg{sna} package installed to work (does not need to be attached).
}
\usage{
make_reduced(adj_list, splitn = 1, weighted = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adj_list}{List of adjacency matrices for use in \code{concor}}
  \item{splitn}{The number of times \code{concor} splits each matrix}
  \item{weighted}{Whether edges in reduced matrices should be weighted (\code{TRUE}) or unweighted (\code{FALSE})}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
A reduced network shows each identified position in the network as a single node. Links (or self-links) exist if the density of ties in that block is greater than the network density.

In the list of input matrices, each should correspond to a different relation for the same nodes. Each matrix will be partitioned with the CONCOR algorithm, continuing through split \code{splitn}. The default \code{weighted} value of \code{FALSE} means that each reduced matrix has values of 0 or 1. If \code{weighted} is \code{TRUE}, the reduced matrices instead stores the density of ties within each block. Comparing density values is only useful if the blocks are of similar sizes, and otherwise can be misleading.
}
\value{
\item{reduced_mat}{A list of reduced matrices, one for each input matrix.}
\item{dens}{A vector of the cut-off densities used (equal to the edge density of each entry in \code{adj_list}).}
}
\references{
%% ~put references to the literature/web site here ~
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
  \code{\link{concor}}, \code{\link{make_blk}}
}
\examples{
g1 <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(g1) <- c("a", "b", "c", "d")
colnames(g1) <- c("a", "b", "c", "d")

g2 <- matrix(c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0), ncol = 4)
rownames(g2) <- c("a", "b", "c", "d")
colnames(g2) <- c("a", "b", "c", "d")

make_reduced(list(g1, g2), splitn = 2)
make_reduced(list(g1, g2), splitn = 2, weighted = TRUE)
}
