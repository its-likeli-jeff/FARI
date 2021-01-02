aci <-
function(a, b) {
  p = 1/2*as.numeric(dist(a,diag=FALSE, upper=FALSE, method="manhattan")) 
  q= 1/2*as.numeric(dist(b,diag=FALSE, upper=FALSE, method="manhattan"))
  
  ndc = 1 - mean( abs( p - q ) )
  
  p = sort(p)
  q = sort(q)
  
  sum.ndc.1 = .Call("sumabsallsortedpairs", p, q )  
  avg.ndc    = 1 - sum(sum.ndc.1)/length(p)^2
  
  ratio = (ndc - avg.ndc)/(1-avg.ndc)    
  list(ndc=ndc, aci=ratio )
}
