\name{plot.blk}
\alias{plot.blk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot.blk
%%  ~~function to do ... ~~
}
\description{
Displays a plot of a blockmodel. Is \code{\link[sna]{plot.blockmodel}} but reformats the shape to be square and removed the "relation - " form the titling of the output. Needs the \code{sna} package installed to work (is not needed to be attached).
}
\usage{
plot.blk(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{blockmodel}}
  \item{...}{Further arguments passed to or from other methods}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Plots a blockmodel that is square.
}
\references{
https://CRAN.R-project.org/package=sna
}
\author{
Carter T. Butts buttsc@uci.edu

edited by Tyme Suda
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[sna]{plot.blockmodel}}
}
\examples{
#Create a random graph with some edge structure
g.p1<-sapply(runif(20,0,1),rep,20)  #Create a matrix of edge probabilities
g1<-rgraph(20,tprob=g.p1)           #Draw from a Bernoulli graph distribution4
colnames(g1) <- paste( LETTERS, "col", sep ="")[1:20]
rownames(g1) <- paste( LETTERS, "col", sep ="")[1:20]

gl=list(g1)
bm=make.blk(gl,2)[[1]]
plot.blk(bm)
}

