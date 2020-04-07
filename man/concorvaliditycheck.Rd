\name{concor_validitycheck}
\alias{concor_validitycheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{concor.validitycheck
%%  ~~function to do ... ~~
}
\description{
Checks whether a list of input matrices are valid for \code{concor}. Specifically makes sure matrix sizes/vertex names match and adds vertex names if they don't exist.
}
\usage{
concor_validitycheck(m_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m_list}{The list of input matrices  for use in \code{concor}
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the same list of matrices if everything checks out, fails if size mismatch or node names mismatch, or if node names are simply missing adds them and warns you.
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
\code{\link{concor}}, \code{\link{name}}
}
\examples{
a <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
b <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
c <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1), ncol = 3)
rownames(a) <- c("a", "b", "c", "d")
colnames(a) <- c("a", "b", "c", "d")

#fails name mismatch
concor_validitycheck(list(a, b))

#fails size mismach
concor_validitycheck(list(b, c))

#passes but adds node names
concor_validitycheck(list(b, b))

#passes
concor_validitycheck(list(a, a))
}

