\name{aci}
\alias{aci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
aci: computes NDC and ACI measures
}
\description{
Efficient computation of fuzzy generalizations of the Rand and adjusted Rand indices.
}
\usage{
aci(a, b)
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
\item{ndc}{Normalized degree of concordance between a and b}
\item{aci}{Adjusted concordance index between a and b}
%% ...
}
\references{
Hullermeier, Eyke, Maria Rifqi, Sascha Henzgen, and Robin Senge. "Comparing fuzzy partitions: A generalization of the rand index and related measures." IEEE Transactions on Fuzzy Systems 20, no. 3 (2011): 546-556.

D'Ambrosio, Antonio, Sonia Amodio, Carmela Iorio, Giuseppe Pandolfo, and Roberta Siciliano. "Adjusted Concordance Index: an Extension of the Adjusted Rand Index to Fuzzy Partitions." Journal of Classification (2020): 1-17.

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
aci(a, b)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ndc}% use one of  RShowDoc("KEYWORDS")
\keyword{aci}% __ONLY ONE__ keyword per line
\keyword{ari}