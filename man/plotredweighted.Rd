\name{plot_red_weighted}
\alias{plot_red_weighted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot_red_weighted
%%  ~~function to do ... ~~
}
\description{
Plot a reduced graph with weighted ties using \code{igraph}.
}
\usage{
plot_red_weighted(blk)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{blk}{The igraph object representing the weighted reduced graph}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Plots the reduced graph with node color being different for each role (and the same for the overall network) and tie thickness representing the weight of that tie.
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
\code{\link{make_reduced}}, \code{\link{make_reduced_igraph}}, \code{\link{plot_red_unweighted}}
}
\examples{
library(igraph)
g1 <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(g1) <- c("a", "b", "c", "d")
colnames(g1) <- c("a", "b", "c", "d")

#make the reduced matrix
r_mat <- make_reduced(list(g1), nsplit = 1, TRUE)
r_igraph <- make_reduced_igraph(r_mat$reduced_mat[[1]])

plot_red_weighted(r_igraph)
}
