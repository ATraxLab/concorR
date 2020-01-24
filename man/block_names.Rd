\name{block_names}
\alias{block_names}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{block_names
%%  ~~function to do ... ~~
}
\description{
A sub-function of \code{\link{concor}}.

Create list a data frames from a list of matrices with column names. The data frames will each include the input matrix as the block, and the column names as vertex.
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
Returns list of data frames each with the first column and arbitrary index number, the second column being the block membership (1 for the first 2 for the second, etc), and the third being the node name (column name from the input matrix).
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
#make the block names for two matrixies and conbine
m1=matrix(c(0,0,0,1,0,1,1,0,0),3)
colnames(m1)=c("a","b", "C")
m2=matrix(c(0,0,0,1),2)
colnames(m2)=c("d","e")

do.call(rbind, block_names(list(m1,m2)))
}

