\name{FARI-package}
\alias{FARI-package}
\docType{package}
\title{
\packageTitle{FARI}
}
\description{
\packageDescription{FARI}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{FARI}
\packageIndices{FARI}
See fari and aci functions
}
\author{
\packageAuthor{FARI}

Maintainer: \packageMaintainer{FARI}
}
\references{
Hullermeier, Eyke, Maria Rifqi, Sascha Henzgen, and Robin Senge. "Comparing fuzzy partitions: A generalization of the rand index and related measures." IEEE Transactions on Fuzzy Systems 20, no. 3 (2011): 546-556.

D'Ambrosio, Antonio, Sonia Amodio, Carmela Iorio, Giuseppe Pandolfo, and Roberta Siciliano. "Adjusted Concordance Index: an Extension of the Adjusted Rand Index to Fuzzy Partitions." Journal of Classification (2020): 1-17.

Andrews, Jeffrey L., Ryan Browne, and Chelsey Hvingelby. "On assessments of agreement between fuzzy partitions". Submitted.
}
\keyword{ package }
\examples{
set.seed(351237342)
a <- matrix(runif(600), nrow=200, ncol=3)
a <- a/rowSums(a)
b <- matrix(runif(600), nrow=200, ncol=3)
b <- b/rowSums(b)
fari(a, b)
aci(a, b)
}
