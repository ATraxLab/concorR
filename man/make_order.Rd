\name{make_order}
\alias{make_order}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{make_order
%%  ~~function to do ... ~~
}
\description{
A sub-function of \code{\link{concor}}.

Combines a list of suborders into an overall order.
}
\usage{
make_order(order_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{order_list}{list of suborders, they are in the order they will be combined in
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the overall order that is combined form the smaller sub-orders
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
#find the order of the top row of two matrixies then conbine for some reason
#is actually usfull in the concor funciton
a=matrix(c(3,2,3,19,1,8,1,1,1),3)
b=matrix(c(4,5,6,7),2)
o1=order(a[1,])
o2=order(b[1,])

make_order(list(o1,o2))

}
