\name{concor1}
\alias{concor1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{concor1
%%  ~~function to do ... ~~
}
\description{
A subfunction of \code{\link{concor}}.

Runs the Pearson correlation on the matrix stack repeatedly until a square matrix consisting of only +1 and -1 is formed. This will only run the correlation until either all values are absolutely greater than that of the cutoff or the maximum number of iterations max_iter is reached.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
concor1(m_stack, cutoff = .9999999, max_iter = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m_stack}{The stack of input matrices, needs to include all transposes and relations}
  \item{cutoff}{The absolute cutoff value for each matrix element}
  \item{max_iter}{The maximum number of times that the correlation can be ran}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a sure matrix with number of rows/columns equal to the number of columns of \code{m_stack} where the values are either +1 or -1, the input matrix can later be split into two "roles" based off the locations of the positive/negative values.
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
