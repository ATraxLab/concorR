\name{concor}
\alias{concor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{concor
%%  ~~function to do ... ~~
}
\description{
Utilize the CONCOR (CONvergence of iterated CORrelations) algorithm to identify roles within social network data. This version works for an arbitrary number of relationships, assuming they all are for the same data (same size of input matrices), and is able to be used with isolates present. It requires further testing on weighted networks (but there is no evidence of it failing on them yet). It will attempt split the network p times, causing there to be 2^p partitions with one extra that is just the isolated nodes (if they exist). If the algorithm attempts to split a structurally equivalent node group, or a singular node group the function will fail. Due to the current function being unable to ignore such a group and continue to split subsequent groups one must lower the number of desired splits, p, to allow the function to succeed.
}
\usage{
concor(m0, cutoff=.9999999, max_iter=50, p=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m0}{A list of adjacency matrices. All matrices MUST be of the square and of the same dimensions, they must also have row/column names (node names). These should correspond to different relationships of the same nodes. If there is only one relationship of interest m0 must be a list of that singular matrix.}
  \item{cutoff}{The cutoff point of the iterated correlations, this defaults to .9999999. If concor seems to be running slowly both this can max_iter can be lowered slightly.}
  \item{max_iter}{The maximum number of times the correlation can be run for each split. If concor seems to be running slowly both this can max_iter can be lowered slightly.}
  \item{p}{The number of times the input matrices will be partitioned. Defaults to p=1 and you will get 2^p roles identified.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The output will be a single large data-frame with two important columns the first being $block which is the block or role identified by concor and the second, $vertex which is vector of the nodes that correspond to each block.
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
rownames(a)=c("a","b","c","d")
colnames(a)=c("a","b","c","d")

concor(list(a))
}

