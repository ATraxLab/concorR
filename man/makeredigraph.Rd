\name{make.red.igraph}
\alias{make.red.igraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{make.red.igraph
%%  ~~function to do ... ~~
}
\description{
A function to turn a adjacency matrix outputted from \code{\link{make.red}} into an \code{igraph} object. This function requires \code{igraph} to function.
}
\usage{
make.red.igraph(red.mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{red.mat}{The outputted reduced network adjacency matrix outputted from \code{\link{make.red}} as \code{$red.mat[[i]]}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns either a weighted or weighted igraph object depending on if the input matrix has any values greater than 1. If such values exist the output will be a weighted igraph object.
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
\code{\link{make.red}}
}
\examples{
library(igraph)
a=matrix(c(0,0,0,1,1,0,1,0,1,1,0,1,1,0,1,0), ncol=4)
rownames(a)=c("a","b","c","d")
colnames(a)=c("a","b","c","d")

#make the reduced matrix
r.mat=make.red(list(a),splitn=1, TRUE)

make.red.igraph(r.mat$red.mat[[1]])
}

