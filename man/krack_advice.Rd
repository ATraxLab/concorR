\name{krack_advice}
\alias{krack_advice}
\alias{krack_friend}
\alias{krack_report}
\docType{data}
\title{
Krackhardt High-Tech Managers data
}
\description{
A network of connections between 21 managers at a manufacturing company taken by David
Krackhardt in 1987. The three networks available are advice, friendship, and report
structure.
}
\usage{
krack_advice

krack_friend

krack_report
}
\format{
A directed and unweighted \pkg{igraph} object.

Vertex attributes:
\itemize{
  \item \code{Name}: Vertex names (character; "v1", "v2", ..., "v21").

  \item \code{x,y}: A set of plotting coordinates, used if no other layout is supplied
  (numeric).

  \item \code{Age}: The age of each manager (integer).

  \item \code{Tenure}: The tenure of each manager (numeric).

  \item \code{Level}: The level in the corporate hierarchy (integer; 1 = CEO, 2 = Vice
  President, 3 = manager).

  \item \code{Department}: What department each manager is in (integer; 1, 2, 3, 4, or 0
  for the CEO).
}
}
\details{
Edges in the \code{krack_advice} network come from the question 'Who would X go to for
help or advice at work?' Edges in the \code{krack_friend} network use the question 'Who
is a friend of X?' The \code{krack_report} network uses the formal organization chart to
define connections between managers.
}
\source{
The data was found at
\href{http://vlado.fmf.uni-lj.si/pub/networks/data/WaFa/default.htm}{http://vlado.fmf.uni-lj.si/pub/networks/data/WaFa/default.htm}
and reformatted for use in R.
}
\references{
D. Krackhardt, Cognitive social structures. \emph{Social Networks}. \strong{9}, 104
(1987). \doi{10.1016/0378-8733(87)90009-8}

S. Wasserman and K. Faust, \emph{Social Network Analysis: Methods and Applications}
(Cambridge University Press, 1994).
}
\keyword{datasets}
