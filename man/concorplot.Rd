\name{concor_plot}
\alias{concor_plot}
\title{Plot network colored by CONCOR partition}
\description{
Plots a network of interest using CONCOR partition as vertex color. Uses an \pkg{igraph}
object, normally created by \code{\link{concor_make_igraph}} or
\code{\link{concor_igraph_apply}}.
}
\usage{
concor_plot(iobject, nsplit = NULL)
}
\arguments{
  \item{iobject}{An \pkg{igraph} object with \code{concor} split as appropriately-named
  vertex attribute.}

  \item{nsplit}{Split number to use as vertex color.}
}
\details{
This is a wrapper to plot an \pkg{igraph} object with usually-readable settings. It looks
for the input \code{iobject} to have a vertex attribute called \code{csplit(nsplit)} that
holds the CONCOR partition assignment (e.g., if \code{nsplit = 2}, then
\code{concor_plot} expects a vertex attribute named \code{csplit2}).
}
\value{
Returns NULL, invisibly.
}

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

