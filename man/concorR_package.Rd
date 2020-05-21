\name{concorR}
\alias{concorR}

\title{The concorR Package}

\description{
\pkg{concorR} implements the CONCOR (CONvergence of iterated CORrelations) algorithm as
introduced by Breiger, Boorman, and Arabie (1975) and detailed by Wasserman and Faust
(1994). This package includes the \code{\link{concor}} algorithm itself, and some useful
functions for plotting and interpreting the outputs. The plotting methods included are:
plotting the overall network with the \code{\link{concor}} partitions as colors (using
the \pkg{igraph} package), plotting a blockmodel using the \code{\link{concor}}
partitions (using the \pkg{sna} package), and plotting reduced network graphs showing the
connections between \code{\link{concor}} partitions (again using \pkg{igraph}).
}

\author{
Tyme Suda \email{suda.4@wright.edu}

Adrienne Traxler \email{adrienne.traxler@wright.edu}
}

\references{
R. L. Breiger, S. A. Boorman, P. Arabie, An algorithm for clustering relational data with
applications to social network analysis and comparison with multidimensional scaling.
\emph{J. of Mathematical Psychology}. \strong{12}, 328 (1975).
\doi{10.1016/0022-2496(75)90028-0}

S. Wasserman and K. Faust, \emph{Social Network Analysis: Methods and Applications}
(Cambridge University Press, 1994).
}
