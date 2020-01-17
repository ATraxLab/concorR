\name{make.red}
\alias{make.red}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{make.red
%%  ~~function to do ... ~~
}
\description{
Takes a list of adjacency matrices, each corresponding to a different relationship for the same data-set and uses the \code{concor} splitting for the desired split \code{splitn} and returns a list of reduced adjacency matrices and their corresponding cutoff densities. Each reduced matrix corresponds to one input relation and can be weighted or unweighted. Needs the \code{sna} package installed to work (is not needed to be attached).
}
\usage{
make.red(adj.list, splitn=1, weighted=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adj.list}{List of adjacency matrices for use in \code{concor}}
  \item{splitn}{The number of times \code{concor} splits the data}
  \item{weighted}{Whether the output should be weighted, \code{TRUE}, or unweighted, \code{FALSE}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list of reduced matrices saved in \code{$red.mat} followed by a vector of the cut-off densities calculated using \code{edge.dens} in \code{$dens}.
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
  \code{\link{concor}}, \code{\link{make.blk}}, \code{\link{edge.dens}}
}
\examples{
#Create a random graph with some edge structure
g.p1<-sapply(runif(20,0,1),rep,20)  #Create a matrix of edge probabilities
g1<-rgraph(20,tprob=g.p1)           #Draw from a Bernoulli graph distribution
colnames(g1) <- paste( LETTERS, "col", sep ="")[1:20]
rownames(g1) <- paste( LETTERS, "col", sep ="")[1:20]

g.p2<-sapply(runif(20,0,1),rep,20)
g2<-rgraph(20,tprob=g.p2)
colnames(g2) <- paste( LETTERS, "col", sep ="")[1:20]
rownames(g2) <- paste( LETTERS, "col", sep ="")[1:20]

make.red(list(g1,g2),splitn=2)
make.red(list(g1,g2),splitn=2, TRUE)

}
