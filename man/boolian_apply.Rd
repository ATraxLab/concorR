\name{boolian_apply}
\alias{boolian_apply}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{boolian_apply
%%  ~~function to do ... ~~
}
\description{
A sub-function of concor.

Take a matrix with n columns and a boolean vector with n entries and remove the columns that correspond to the false entries to of the vector.
}
\usage{
boolian_apply(boolian, mat_stack)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{boolian}{The Boolean vector}
  \item{mat_stack}{The matrix}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the shrunken matrix
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
b=c(TRUE, TRUE, FALSE, TRUE)

boolian_apply(b,a)
}

