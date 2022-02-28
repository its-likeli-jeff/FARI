# FARI
An R package for computing efficient versions of the Adjusted Concordance Index (D’Ambrosio  et  al., 2020) and the Frobenius Adjusted Rand Index (Andrews, Browne, and Hvingelby, to appear). Both compute fuzzy generalizations of the Adjusted Rand Index (Hubert and Arabie, 1985).


# Installation
library(devtools)

install_github("its-likeli-jeff/FARI")

library(FARI)

#two functions, aci(a,b) and fari(a,b)


# References
Andrews, Jeffrey L., Ryan Browne, and Chelsey Hvingelby. "On Assessments of Agreement Between Fuzzy Partitions". Journal of Classification (to appear, available online). https://link.springer.com/article/10.1007/s00357-021-09407-3

D'Ambrosio, Antonio, Sonia Amodio, Carmela Iorio, Giuseppe Pandolfo, and Roberta Siciliano. "Adjusted Concordance Index: an Extension of the Adjusted Rand Index to Fuzzy Partitions." Journal of Classification (2020): 1-17.

Hubert, Lawrence, and Phipps Arabie. "Comparing partitions." Journal of classification 2, no. 1 (1985): 193-218.
