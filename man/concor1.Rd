\name{concor1}
\alias{concor1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find +1/-1 convergence points of an input matrix
%%  ~~function to do ... ~~
}
\description{
A subfunction of \code{\link{concor}}.

Runs the Pearson correlation on the matrix stack repeatedly until a square matrix consisting of only +1 and -1 is formed.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
concor1(m_stack, cutoff = .9999999, max_iter = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m_stack}{The stack of input matrices, needs to include all transposes and relations}
  \item{cutoff}{The absolute value of the convergence threshold for each matrix element}
  \item{max_iter}{Maximum number of times to run the correlation while seeking convergence}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
For network data with \code{R} relations, the input matrix "stack" is a \code{2NR x N} matrix, which consists of each relation's adjacency matrix, then the transpose of that matrix, appended to each other vertically. The correlation is run until either all matrix entries have absolute values greater than \code{cutoff} or the maximum number of iterations \code{max_iter} is reached. On its own, this function does not execute the whole CONCOR method, but is listed separately from \code{concor} to demonstrate the core step.
}
\value{
Returns a square matrix with number of rows/columns equal to the number of columns of \code{m_stack} where the values are either +1 or -1. The input matrix can later be split into two "roles" based off the locations of the positive/negative values.
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
a <- matrix(c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), ncol = 4)
b <- rbind(a, t(a))

concor1(b)
}
