\name{concor_plot}
\alias{concor_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot network colored by CONCOR partition
%%  ~~function to do ... ~~
}
\description{
Plots a network of interest with CONCOR partition as vertex color. The \pkg{igraph} object used should be created using \code{\link{concor_make_igraph}} or \code{\link{concor_igraph_apply}}.
}
\usage{
concor_plot(iobject, p = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iobject}{An \pkg{igraph} object with \code{concor} splitting as vertex attribute. Should be created using \code{\link{concor_make_igraph}} or \code{\link{concor_igraph_apply}}}
  \item{p}{Split number of interest; \code{csplitp} should be a vertex attribute of \code{iobject}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This is a wrapper to plot an \pkg{igraph} object with usually-readable settings. It looks for the input \code{iobject} to have a vertex attribute \code{csplitp} that holds the CONCOR partition assignment for split \code{p}.
}
\value{
Returns NULL, invisibly.
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
\code{\link{concor}}, \code{\link{concor_make_igraph}}, \code{\link{concor_igraph_apply}}
}
\examples{
a <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
rownames(a) <- c("a", "b", "c", "d")
colnames(a) <- c("a", "b", "c", "d")

i_out <- concor_make_igraph(list(a))
concor_plot(i_out[[1]], 1)
}

