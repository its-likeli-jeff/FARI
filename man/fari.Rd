\name{fari}
\alias{fari}
\alias{FARI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fari: Computes the Frobenius Rand index and Frobenius adjusted Rand index
}
\description{
Efficient computation of fuzzy generalizations of the Rand and adjusted Rand indices.
}
\usage{
fari(a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
an n (sample) by G1 (number of groups) matrix of hard or fuzzy cluster memberships
}
  \item{b}{
an n (sample) by G2 (number of groups) matrix of hard or fuzzy cluster memberships
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{fri}{Measure of Frobenius Rand index between a and b}
\item{fari}{Measure of Frobenius adjusted Rand index between a and b}
%% ...
}
\references{
Andrews, Jeffrey L., Ryan Browne, and Chelsey Hvingelby. "On assessments of agreement between fuzzy partitions". Submitted.
}
\author{
Jeffrey L. Andrews and Ryan Browne
}

\examples{
set.seed(351237342)
a <- matrix(runif(600), nrow=200, ncol=3)
a <- a/rowSums(a)
b <- matrix(runif(600), nrow=200, ncol=3)
b <- b/rowSums(b)
fari(a, b)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fri}% use one of  RShowDoc("KEYWORDS")
\keyword{ fari }% __ONLY ONE__ keyword per line
\keyword{ ari }