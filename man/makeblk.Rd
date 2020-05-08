\name{make_blk}
\alias{make_blk}
\title{Make \pkg{blockmodel} objects using CONCOR partition}
\description{
Run \code{concor} for a specified split and output \pkg{blockmodel} object(s).
}
\usage{
make_blk(adj_list, splitn = 1)
}
\arguments{
  \item{adj_list}{A list of adjacency matrices, each representing a different relation for the same nodes}
  \item{splitn}{The total number of splits CONCOR will preform}
}
\details{
Run \code{concor} for a specified split, \code{splitn}, and create a \pkg{blockmodel} type object from each input matrix in \code{adj_list} using the membership vector created through \code{concor}.
}
\value{
A list of \pkg{blockmodel} type objects, one for each input adjacency matrix, using the split returned by \code{concor}. Each blockmodel in the list corresponds to one relation. See the help file for \code{\link[sna]{blockmodel}} in the \pkg{sna} package for more details.
}

\seealso{
  \code{\link[sna]{blockmodel}},
  \code{\link{concor}}
}
\examples{
g1 <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(g1) <- c("a", "b", "c", "d")
colnames(g1) <- c("a", "b", "c", "d")

g2 <- matrix(c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0), ncol = 4)
rownames(g2) <- c("a", "b", "c", "d")
colnames(g2) <- c("a", "b", "c", "d")

gl <- list(g1, g2)
make_blk(gl)
}

