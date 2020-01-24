\name{make_sub_boolian}
\alias{make_sub_boolian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{make_sub_boolian
%%  ~~function to do ... ~~
}
\description{
A sub-function of \code{\link{concor}}.

Makes a Boolean vector from the first column of a matrix where TRUE a positive value and FALSE is negative.
}
\usage{
make_sub_boolian(cor_mat_orderd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cor_mat_orderd}{The  input matrix
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the Boolean vector.
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
a=matrix(c(0,0,0,1,1,0,1,0,1,1,0,1,1,0,1,0), ncol=4)
b=rbind(a, t(a))
c=concor1(b)
c

make_sub_boolian(c)
}

