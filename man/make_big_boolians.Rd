\name{make_big_boolians}
\alias{make_big_boolians}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{make_big_boolians
%%  ~~function to do ... ~~
}
\description{
A sub-function of \code{\link{concor}}, you will probably never use this function outside of the \code{\link{concor}} function.

A function to take a list of booleans in a specified order and to make large booleans that are the length of all the boolean vectors combined. These large booleans are all \code{FALSE} except for the location of one of the booleans.
}
\usage{
make_big_boolians(bool_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bool_list}{The list of input booleans, they should be in the order desired.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list of booleans of length double the number of input booleans each the length of all input booleans combined. Each element of the output will be the combination of one of the booleans (or its notted form) placed within a vector of FALSEs at the location of the boolean of interest in the input list.
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
b1=c(FALSE,FALSE,TRUE)
b2=c(TRUE,FALSE,FALSE,TRUE)

make_big_boolians(list(b1,b2))
}

