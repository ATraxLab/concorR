\name{make_blk}
\alias{make_blk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{make_blk
%%  ~~function to do ... ~~
}
\description{
Run concor for a specified split, \code{splitn}, and create a \code{blockmodel} type object from each input matrix in \code{adj.list} using the membership vector created through \code{concor}.
}
\usage{
make.blk(adj.list, splitn=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adj.list}{The list of input matrices to CONCOR, each representing a specific relationship}
  \item{splitn}{The total number of splits CONCOR will preform}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list of blockmodel type objects based off the input matrices and the concor split. Each blockmodel in the list corresponds to one relationship. See the help file for \code{\link[sna]{blockmodel}} in the sna package for more details.
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
  \code{\link[sna]{blockmodel}},
  \code{\link{concor}}
}
\examples{
g1=matrix(c(0,0,0,1,1,0,1,0,1,1,0,1,1,0,1,0), ncol=4)
rownames(g1)=c("a","b","c","d")
colnames(g1)=c("a","b","c","d")

g2=matrix(c(0,1,0,0,0,0,1,0,0,0,1,1,1,1,0,0), ncol=4)
rownames(g2)=c("a","b","c","d")
colnames(g2)=c("a","b","c","d")

gl=list(g1,g2)
make.blk(gl)
}

