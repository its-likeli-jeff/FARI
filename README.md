# FARI
An R package for computing efficient versions of the Adjusted Concordance Index (D’Ambrosio  et  al., 2020) and the Frobenius Adjusted Rand Index (Andrews, Browne, and Hvingelby, Submitted). Both compute fuzzy generalizations of the Adjusted Rand Index (Hubert and Arabie, 1985).


# Installation
library(devtools)

install_github("its-likeli-jeff/FARI")

library(FARI)

#two functions, aci(a,b) and fari(a,b)


# References
Andrews, Jeffrey L., Ryan Browne, and Chelsey Hvingelby. Submitted.

D'Ambrosio, Antonio, Sonia Amodio, Carmela Iorio, Giuseppe Pandolfo, and Roberta Siciliano. "Adjusted Concordance Index: an Extension of the Adjusted Rand Index to Fuzzy Partitions." Journal of Classification (2020): 1-17.

Hubert, Lawrence, and Phipps Arabie. "Comparing partitions." Journal of classification 2, no. 1 (1985): 193-218.
