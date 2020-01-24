\name{plot.red.unweighted}
\alias{plot.red.unweighted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot.red.unweighted
%%  ~~function to do ... ~~
}
\description{
Plot a reduced graph with unweighted ties using \code{igraph}.
}
\usage{
plot.red.unweighted(blk)
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
\code{\link{make.red}}, \code{\link{make.red.igraph}}, \code{\link{plot.red.weighted}}
}
\examples{
library(igraph)
a=matrix(c(0,0,0,1,1,0,1,0,1,1,0,1,1,0,1,0), ncol=4)
rownames(a)=c("a","b","c","d")
colnames(a)=c("a","b","c","d")

#make the reduced matrix
r.mat=make.red(list(a),splitn=1, FALSE)
r.igraph=make.red.igraph(r.mat$red.mat[[1]])

plot.red.unweighted(r.igraph)
}
