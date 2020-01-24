\name{val.diag}
\alias{val.diag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{val.diag
%%  ~~function to do ... ~~
}
\description{
A subfunction of \code{\link{concor}}.

Sets the value of the diagonal of a matrix, for use in apply functions.
}
\usage{
valdiag(m, value=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{The matrix}
  \item{value}{The value for the diagonals, defaults to \code{NA}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the matrix with the new diagonal values.
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
\code{\link{concor}}
}
\examples{
a=matrix(c(1,2,3,4),2)
b=matrix(c(4,5,6,7),2)
c=list(a,b)

lapply(c, function(x) val.diag(x,0))
}
