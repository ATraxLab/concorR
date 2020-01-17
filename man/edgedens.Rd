\name{edge.dens}
\alias{edge.dens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{edge.dens
%%  ~~function to do ... ~~
}
\description{
Calculate the edge density of an adjacency matrix where self ties are not possible.
}
\usage{
edge.dens(adj.mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adj.mat}{The adjacency matrix of interest
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the numeric value for edge density, for unweighted matrices this is between 0 and 1.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Create a random graph with some edge structure
g.p1<-sapply(runif(20,0,1),rep,20)  #Create a matrix of edge probabilities
g1<-rgraph(20,tprob=g.p1)           #Draw from a Bernoulli graph distribution

edge.dens(g1)
}

