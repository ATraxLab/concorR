\name{make.blk}
\alias{make.blk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{make.blk
%%  ~~function to do ... ~~
}
\description{
Run concor for a specified split, \code{splitn}, and create a \code{blockmodel} type object from each input matrix in \code{adj.list} using the membership vector created through \code{concor}. Needs the \code{sna} package installed to work (is not needed to be attached).
}
\usage{
make.blk(adj.list, splitn=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adj.list}{The list of input matrices to CONCOR, each representing a specific relationship}
  \item{splitn}{The total number of splits CONCOR will preform}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list of blockmodel type objects based off the input matrices and the concor split. Each blockmodel in the list corresponds to one relationship. See the help file for \code{\link[sna]{blockmodel}} in the sna package for more details.
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
  \code{\link[sna]{blockmodel}},
  \code{\link{concor}}
}
\examples{
#Create a random graph with some edge structure
g.p1<-sapply(runif(20,0,1),rep,20)  #Create a matrix of edge probabilities
g1<-rgraph(20,tprob=g.p1)           #Draw from a Bernoulli graph distribution4
colnames(g1) <- paste( LETTERS, "col", sep ="")[1:20]
rownames(g1) <- paste( LETTERS, "col", sep ="")[1:20]

g.p2<-sapply(runif(20,0,1),rep,20)
g2<-rgraph(20,tprob=g.p2)
colnames(g2) <- paste( LETTERS, "col", sep ="")[1:20]
rownames(g2) <- paste( LETTERS, "col", sep ="")[1:20]

gl=list(g1,g2)
make.blk(gl)
}

