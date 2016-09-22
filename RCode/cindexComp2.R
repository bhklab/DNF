


cindexComp2 <-
  function(cindex1, cindex2, x1, x2) {
    
    if(cindex1['n'] != cindex2['n']) { stop("the concordance indices are computed from different number of samples!") }
    cindex1SE <- cindex1['S.D.']/2
    cindex2SE <- cindex2['S.D.']/2
    if(is.na(cindex1SE) || is.na(cindex2SE)){stop("the concordance indices must be computed using method noether!")}
    eps <- 1E-15
    
    n <- cindex1['n']
    r <- cor(x1, x2, use="complete.obs", method="spearman")
    if((1 - abs(r)) > eps) {
      t.stat <- (cindex1['C Index'] - cindex2['C Index']) / sqrt(cindex1SE^2 + cindex2SE^2 - 2 * r * cindex1SE * cindex2SE)
      diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
    } else { diff.ci.p <- 1 }
    return(list("p.value"=diff.ci.p, "cindex1"=cindex1['C Index'], "cindex2"=cindex2['C Index']))
  }
