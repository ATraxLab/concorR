\name{make_reduced_igraph}
\alias{make_reduced_igraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{make_reduced_igraph
%%  ~~function to do ... ~~
}
\description{
Turns a reduced adjacency matrix (usually output from \code{\link{make_reduced}}) into an \pkg{igraph} object. This function requires \pkg{igraph} to work.
}
\usage{
make_reduced_igraph(reduced_mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reduced_mat}{A reduced network adjacency matrix, outputted from \code{\link{make_reduced}} as \code{$reduced_mat[[i]]}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Either a weighted or unweighted \pkg{igraph} object depending on if the input matrix has any values greater than 1. If such values exist, the output will be a weighted igraph object.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{make_reduced}}
}
\examples{
a <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(g1) <- c("a", "b", "c", "d")
colnames(g1) <- c("a", "b", "c", "d")

#make the reduced matrix
r_mat <- make_reduced(list(a), splitn = 1, TRUE)

make_reduced_igraph(r_mat$reduced_mat[[1]])
}

