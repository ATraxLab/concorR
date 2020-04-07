\name{plot_red_unweighted}
\alias{plot_red_unweighted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot_red_unweighted
%%  ~~function to do ... ~~
}
\description{
Plot a reduced graph with unweighted ties using \code{igraph}.
}
\usage{
plot_red_unweighted(blk)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{blk}{The igraph object representing the unweighted reduced graph}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Plots the reduced graph with node color being different for each role (and the same for the overall network)
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
\code{\link{make_reduced}}, \code{\link{make_reduced_igraph}}, \code{\link{plot_red_weighted}}
}
\examples{
library(igraph)
g1 <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(g1) <- c("a", "b", "c", "d")
colnames(g1) <- c("a", "b", "c", "d")

#make the reduced matrix
r_mat <- make_reduced(list(g1), splitn = 1, FALSE)
r_igraph <- make_reduced_igraph(r_mat$reduced_mat[[1]])

plot_red_unweighted(r_igraph)
}
