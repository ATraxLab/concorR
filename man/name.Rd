\name{name}
\alias{name}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{name
%%  ~~function to do ... ~~
}
\description{
gives names \code{v001}, \code{v002}... to each column/row of an input matrix. These are treated as the node names of the network represented by that matrix.
}
\usage{
name(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{The matrix with no names
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns the input matrix with new column/row names, \code{v001}, \code{v002}...
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
a=matrix(c(0,0,0,1,1,0,1,0,1,1,0,1,1,0,1,0), ncol=4)
name(a)
}

