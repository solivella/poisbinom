\name{Poisson-Binomial}
\alias{dpoisbinom}
\alias{qpoisbinom}
\alias{ppoisbinom}
\alias{rpoisbinom}
\title{The Poisson Binomial Distribution}
\description{
Probability mass, distribution, quantile and function, and random number generator for the Poisson-Binomial distribution with parameter vector \code{pp} (the probability parameter of the component Binomial random variables).  
}
\usage{
dpoisbinom(x, pp, log_d = FALSE)
ppoisbinom(q, pp, lower_tail = TRUE, log_p = FALSE)
qpoisbinom(p, pp, lower_tail = TRUE, log_p = FALSE)
rpoisbinom(n,pp)
}
\arguments{
  \item{x, q}{
  vector of quantiles.
}
\item{p, pp}{
vector of probabilities.
}
\item{n}{
number of random deviates.
}
\item{log_d, log_p}{
logical; if TRUE, probabilities are given in the log scale.
}
\item{lower_tail}{
logical; if TRUE (default), probabilities are \eqn{Pr(X \leq x)},
otherwise, \eqn{Pr(X > x)}.
}
}
\details{
The Poisson-Binomial distribution is the distribution of a sum of \eqn{n} independent and \emph{not} identically distributed Binomial random variables. It is parameterized by the vector of \eqn{n} possibly distinct probability parameters of these Binomial distributions, and is computed using a discrete Fourier transform. See Hong (2013) for details.
}
\value{
\code{dpoisbinom} gives the mass, \code{ppoisbinom} gives the distribution function, \code{qpoisbinom} gives the quantile function and \code{rpoisbinom} generates random deviates. 

If \code{pp} contains values outside of [\eqn{0,1}], an error is returned. 

The length of the result is determined by \code{n} in \code{rpoisbinom}, and is the length of the first argument for all other functions.  
}
\references{
Hong, Y. (2013) ``On computing the distribution function for the Poisson
binomial distribution''. \emph{Computational Statistics and Data Analysis}, 59, 41--51.
}
\author{
Shiraito, Y. and Olivella, S. (2017).}

\examples{
## Binomial probabilities
pp <- runif(500)

## PMF
dpoisbinom(36, pp)

## CDF
ppoisbinom(36, pp)

## Quantile function
qpoisbinom(0.3, pp)

## Random deviates
rpoisbinom(5, pp)
}
\keyword{ distribution }% use one of  RShowDoc("KEYWORDS")
