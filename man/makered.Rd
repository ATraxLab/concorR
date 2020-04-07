\name{make_reduced}
\alias{make_reduced}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{make_reduced
%%  ~~function to do ... ~~
}
\description{
Takes a list of adjacency matrices, each corresponding to a different relationship for the same data-set and uses the \code{concor} splitting for the desired split \code{splitn} and returns a list of reduced adjacency matrices and their corresponding cutoff densities. Each reduced matrix corresponds to one input relation and can be weighted or unweighted. Needs the \code{sna} package installed to work (is not needed to be attached).
}
\usage{
make_reduced(adj_list, splitn = 1, weighted = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adj_list}{List of adjacency matrices for use in \code{concor}}
  \item{splitn}{The number of times \code{concor} splits the data}
  \item{weighted}{Whether the output should be weighted, \code{TRUE}, or unweighted, \code{FALSE}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list of reduced matrices saved in \code{$reduced_mat} followed by a vector of the cut-off densities calculated (edge densities) in \code{$dens}.
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
make_reduced(list(g1, g2), splitn = 2, TRUE)
}
