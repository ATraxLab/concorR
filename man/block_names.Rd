\name{block_names}
\alias{block_names}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{block_names
%%  ~~function to do ... ~~
}
\description{
A sub-function of concor.

Create a data frame from a list of matrices with column names. The data frame will include each input matrix as the block, and the column names as vertex.
}
\usage{
block_names(mat_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat_list}{The list of input matrices. They must all include column names.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a data frame with the first column and arbitrary index number, the second column being the block membership, and the third being the node name (column name).
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
m1=matrix(c(0,0,0,1,0,1,1,0,0),3)
colnames(m1)=c("a","b", "C")
m2=matrix(c(0,0,0,1),2)
colnames(m2)=c("d","e")

block_names(list(m1,m2))
}

