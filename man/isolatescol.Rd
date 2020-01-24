\name{isolates.col}
\alias{isolates.col}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{isolates.col
%%  ~~function to do ... ~~
}
\description{
A subfunction of\code{\link{concor}}.

Find all empty columns (all values zero), for a stack of an adjacency matrix and its transpose these will correspond to the isolates in the network.
}
\usage{
isolates.col(m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{The input matrix, must have column names
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the column names of the empty columns.
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
m=matrix(c(0,0,0,1,0,1,1,0,0),3)
colnames(m)=c("a","b", "C")
rownames(m)=c("a","b", "C")
m

isolates.col(m)
}

