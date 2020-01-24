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
g1=matrix(c(0,0,0,1,1,0,1,0,1,1,0,1,1,0,1,0), ncol=4)

edge.dens(g1)
}

