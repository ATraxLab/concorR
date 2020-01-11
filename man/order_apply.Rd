\name{order_apply}
\alias{order_apply}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{order_apply
%%  ~~function to do ... ~~
}
\description{
A sub-function of concor.

A simple function for rearranging the row and columns of a matrix.
}
\usage{
order_apply(order,mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{order}{The new order for the rows/columns of the matrix}
  \item{mat}{The matrix}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the rearranged matrix
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
o=c(1,3,2,4)
a=matrix(c(0,0,0,1,1,0,1,0,1,1,0,1,1,0,1,0), ncol=4)
colnames(a)=c("a","b","C","D")
rownames(a)=c("a","b","C","D")

order_apply(o,a)
}
