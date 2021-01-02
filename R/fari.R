fari <-
function(a, b){
  #A, B bonding matrices
  n <- nrow(a)
  A <- a %*% t(a)
  B <- b %*% t(b)
  j <- matrix(1, n, n)
  Na <- (sum(A*j)/sum(A*A))*A
  Nb <- (sum(B*j)/sum(B*B))*B
  ri <- (sum(Na*Nb) + sum((j-Na)*(j-Nb)) - n)/(2*choose(n,2))
  M <- j/n
  R <- diag(n) - M
  Eri <- ((2*sum(A*j)*sum(B*j)/(sum(A*A)*sum(B*B))) * (sum(M*A)*sum(M*B)+((1/(n-1))*sum(R*A)*sum(R*B))) - (sum(A*j)^2)/sum(A*A) - (sum(B*j)^2)/sum(B*B) + n^2 - n)/(2*choose(n,2))
  ari <- (ri - Eri)/(1 - Eri)
  store <- list(fri=ri, fari=ari)
  store
}
