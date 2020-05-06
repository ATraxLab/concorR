\name{krack_report}
\alias{krack_report}
\docType{data}
\title{
Krackhardt High-Tech Managers: Reporting
}
\description{
A network of connections between 21 managers at a high-tech company taken by David Krackhardt in 1987. This network (\code{krack_report}) uses who each manager reports to as to define connections between managers.
}
\usage{krack_advice}
\format{
An directed and unweighted \code{igraph} object, including each of the 21 managers.

Vertex Attributes:

  \code{Name}: Vertex names (v1, v2, ..., v21)

  \code{x}: A stet of x coordinates, the same for \code{krack_friend}, \code{krack_advice}, and \code{krack_report}

  \code{y}: A stet of y coordinates, the same for \code{krack_friend}, \code{krack_advice}, and \code{krack_report}

  \code{Age}: The age of each manager

  \code{Tenure}: The tenure of each manager

  \code{Level}: The level in the corporate hierarchy (1 = CEO, 2 = Vice President, 3 = manager)

  \code{Departmen}t: What department each manager is in (1, 2, 3, 4, or 0 for the CEO)
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
The data was found at \href{http://vlado.fmf.uni-lj.si/pub/networks/data/WaFa/default.htm}{http://vlado.fmf.uni-lj.si/pub/networks/data/WaFa/default.htm} and reformat for use in R.
}
\references{
Krackhardt D. (1987). Cognitive social structures. Social Networks, 9, 104-134.

Stanley Wasserman, Katherine Faust: Social Network Analysis: Methods and Applications. CUP, 1994.
}
\examples{
}
\keyword{datasets}
