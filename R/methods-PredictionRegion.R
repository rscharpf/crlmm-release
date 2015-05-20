setMethod("[", signature(x="PredictionRegion"),
	  function(x, i, j, ..., drop=FALSE){
		  if(missing(i) & missing(j)) return(x)
		  if(!missing(i)){
			  xx <- lapply(x, function(x, i){
				  x$mu <- x$mu[i,,,drop=FALSE]
				  x$cov <- x$cov[i, , , drop=FALSE]
				  return(x)
			  }, i=i)
			  x <- as(xx, "PredictionRegion")
		  }
		  return(x)
	  })
