\name{plot.red.weighted}
\alias{plot.red.weighted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot.red.weighted
%%  ~~function to do ... ~~
}
\description{
Plot a reduced graph with weighted ties using \code{igraph}.
}
\usage{
plot.red.weighted(x)
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
\code{\link{make.red}}, \code{\link{make.red.igraph}}, \code{\link{plot.red.unweighted}}
}
\examples{
library(igraph)
a=matrix(c(0,0,0,1,1,0,1,0,1,1,0,1,1,0,1,0), ncol=4)
rownames(a)=c("a","b","c","d")
colnames(a)=c("a","b","c","d")

#make the reduced matrix
r.mat=make.red(list(a),splitn=1, TRUE)
r.igraph=make.red.igraph(r.mat$red.mat[[1]])

plot.red.weighted(r.igraph)
}
