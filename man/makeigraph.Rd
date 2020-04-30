\name{concor_make_igraph}
\alias{concor_make_igraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{concor_make_igraph
%%  ~~function to do ... ~~
}
\description{
Takes the list of adjacency matrices used in \code{concor} and creates a list of \code{igraph} objects with a specific concor split added as an vertex attribute.
}
\usage{
concor_make_igraph(adj_list, nsplit = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adj_list}{The list of input matrices used for \code{concor}}
  \item{nsplit}{The split number desired, defaults to 1}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list of igraph objects with a vertex attribute \code{csplit_n} where \code{n} is \code{nsplit} that contain the block information obtained from \code{concor}.
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
\code{\link{concor}}, \code{\link{concor_igraph_apply}}
}
\examples{
a <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(a) <- c("a", "b", "c", "d")
colnames(a) <- c("a", "b", "c", "d")

concor_make_igraph(list(a))
}
