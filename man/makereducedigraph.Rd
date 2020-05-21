\name{make_reduced_igraph}
\alias{make_reduced_igraph}
\title{Build an \code{igraph} object for a reduced network}
\description{
Turns a reduced adjacency matrix (usually output from \code{\link{make_reduced}}) into an
\code{igraph} object. This function requires \pkg{igraph} to work.
}
\usage{
make_reduced_igraph(reduced_mat)
}
\arguments{
  \item{reduced_mat}{A reduced network adjacency matrix (typically outputted from
  \code{\link{make_reduced}} in the \code{reduced_mat} list)}
}
\value{
A directed and unweighted \code{igraph} object for the reduced matrix.
}

\seealso{
\code{\link{make_reduced}}
}
\examples{
a <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(a) <- c("a", "b", "c", "d")
colnames(a) <- c("a", "b", "c", "d")

r_mat <- make_reduced(list(a), nsplit = 1)
make_reduced_igraph(r_mat$reduced_mat[[1]])
}

