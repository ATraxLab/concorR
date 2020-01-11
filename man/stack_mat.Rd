\name{stack_mat}
\alias{stack_mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{stack_mat
%%  ~~function to do ... ~~
}
\description{
A subfunction of concor.

Stacks all matrices in a list and their transposes atop each other.
}
\usage{
stack_mat(m_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m_list}{The list of input matrices
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a very tall matrix formed from each input matrix stacked atop its transpose
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
a=matrix(c(1,2,3,4),2)
b=matrix(c(4,5,6,7),2)
c=list(a,b)

d=stack_mat(c)
d
}
