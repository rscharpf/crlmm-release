

getSigma2 <- function(object){
	phi.a <- phiA(object)
	phi.b <- phiB(object)
	cbind(phi.a, phi.b)
}

getCor <- function(object, gt){
	if(gt =="NULL") return(rep(0, nrow(object)))
	assayDataElement(batchStatistics(object), paste("corr", gt, sep=""))
}

getTau2 <- function(object, gt){
	bs <- batchStatistics(object)
	nma <- paste("tau2A", gt, sep=".")
	nmb <- paste("tau2B", gt, sep=".")
	tau2.a <- assayDataElement(bs, nma)
	tau2.b <- assayDataElement(bs, nmb)
	cbind(tau2.a, tau2.b)
}

getVar <- function(object, nm){
	bs <- batchStatistics(object)
	assayDataElement(bs, nm)
}

setMethod("nuA", signature=signature(object="CNSet"), function(object) nu(object, "A"))
setMethod("nuB", signature=signature(object="CNSet"), function(object) nu(object, "B"))

setMethod("phiA", signature=signature(object="CNSet"), function(object) phi(object, "A"))
setMethod("phiB", signature=signature(object="CNSet"), function(object) phi(object, "B"))
setMethod("phiPrimeA", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "phiPrimeA")
})
setMethod("phiPrimeB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "phiPrimeB")
})
setMethod("tau2A", signature=signature(object="CNSet"), function(object) tau2(object, "A"))
setMethod("tau2B", signature=signature(object="CNSet"), function(object) tau2(object, "B"))

setMethod("Ns", signature=signature(object="CNSet"),
	   function(object, ...){
		   Ns(batchStatistics(object), ...)
	   })
setMethod("corr", signature=signature(object="CNSet"),
	   function(object, ...){
		   corr(batchStatistics(object), ...)
	   })

setMethod("medians", signature=signature(object="CNSet"),
	   function(object, ...){
		   medians(batchStatistics(object), ...)
	   })
setMethod("mads", signature=signature(object="CNSet"),
	   function(object, ...){
		   mads(batchStatistics(object), ...)
	   })
setMethod("tau2", signature=signature(object="CNSet"),
	   function(object, ...){
		   tau2(batchStatistics(object), ...)
	   })
##---------------------------------------------------------------------------
## Number of samples with Genotype AA, AB, or BB by batch
setMethod("N.AA", signature=signature(object="CNSet"), function(object){
	assayDataElement(batchStatistics(object), "N.AA")
})
setMethod("N.AB", signature=signature(object="CNSet"), function(object){
	assayDataElement(batchStatistics(object), "N.AB")
})
setMethod("N.BB", signature=signature(object="CNSet"), function(object){
	assayDataElement(batchStatistics(object), "N.BB")
})

setReplaceMethod("N.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "N.AA", value)
	  })
setReplaceMethod("N.AB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "N.AB", value)
	  })
setReplaceMethod("N.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "N.BB", value)
	  })


##---------------------------------------------------------------------------
##  median intensity by genotype cluster for each allele
##should we update the entire matrix rather than one column...
setMethod("medianA.AA", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "medianA.AA")
})
setMethod("medianA.AB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "medianA.AB")
})
setMethod("medianA.BB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "medianA.BB")
})
setMethod("medianB.AA", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "medianB.AA")
})
setMethod("medianB.AB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "medianB.AB")
})
setMethod("medianB.BB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "medianB.BB")
})
setReplaceMethod("medianA.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "medianA.AA", value)
	  })
setReplaceMethod("medianA.AB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "medianA.AB", value)
	  })
setReplaceMethod("medianA.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "medianA.BB", value)
	  })
setReplaceMethod("medianB.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "medianB.AA", value)
	  })
setReplaceMethod("medianB.AB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "medianB.AB", value)
	  })
setReplaceMethod("medianB.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "medianB.BB", value)
	  })

##---------------------------------------------------------------------------
##  mad intensity by genotype cluster for each allele
setMethod("madA.AA", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "madA.AA")
})
setMethod("madA.AB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "madA.AB")
})
setMethod("madA.BB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "madA.BB")
})
setMethod("madB.AA", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "madB.AA")
})
setMethod("madB.AB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "madB.AB")
})
setMethod("madB.BB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "madB.BB")
})
setReplaceMethod("madA.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "madA.AA", value)
	  })
setReplaceMethod("madA.AB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "madA.AB", value)
	  })
setReplaceMethod("madA.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "madA.BB", value)
	  })
setReplaceMethod("madB.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "madB.AA", value)
	  })
setReplaceMethod("madB.AB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "madB.AB", value)
	  })
setReplaceMethod("madB.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "madB.BB", value)
	  })

##---------------------------------------------------------------------------
##  mad of log(intensities) by genotype cluster for each allele
setMethod("tau2A.AA", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "tau2A.AA")
})
setMethod("tau2A.BB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "tau2A.BB")
})
setMethod("tau2B.AA", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "tau2B.AA")
})
setMethod("tau2B.BB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "tau2B.BB")
})
setReplaceMethod("tau2A.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "tau2A.AA", value)
	  })
setReplaceMethod("tau2A.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "tau2A.BB", value)
	  })
setReplaceMethod("tau2B.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "tau2B.AA", value)
	  })
setReplaceMethod("tau2B.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "tau2B.BB", value)
	  })

## correlation of log2A and log2B within each genotype cluster
setMethod("corrAA", signature=signature(object="CNSet"), function(object){
	assayDataElement(batchStatistics(object), "corrAA")
  })
setMethod("corrAB", signature=signature(object="CNSet"), function(object){
	assayDataElement(batchStatistics(object), "corrAB")
  })
setMethod("corrBB", signature=signature(object="CNSet"), function(object){
	assayDataElement(batchStatistics(object), "corrBB")
  })
setReplaceMethod("corrAA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "corrAA", value)
})

setReplaceMethod("corrAB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "corrAB", value)
})

setReplaceMethod("corrBB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "corrBB", value)
})

setReplaceMethod("phiPrimeA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "phiPrimeA", value)
})

setReplaceMethod("phiPrimeB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "phiPrimeB", value)
})

##setMethod("mad.AA", signature=signature(object="CNSet"), function(object) mads(object, "AA"))
##setMethod("mad.AB", signature=signature(object="CNSet"), function(object) mads(object, "AB"))
##setMethod("mad.BB", signature=signature(object="CNSet"), function(object) mads(object, "BB"))
##

##setReplaceMethod("median.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
##	  function(object, value){
##		  linearParamElementReplace(object, "median.AA", value)
##	  })
##setReplaceMethod("median.AB", signature=signature(object="CNSet", value="ff_or_matrix"),
##	  function(object, value){
##		  linearParamElementReplace(object, "median.AB", value)
##	  })
##setReplaceMethod("median.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
##	  function(object, value){
##		  linearParamElementReplace(object, "median.BB", value)
##	  })
##setReplaceMethod("mad.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
##	  function(object, value){
##		  linearParamElementReplace(object, "mad.AA", value)
##	  })
##setReplaceMethod("mad.AB", signature=signature(object="CNSet", value="ff_or_matrix"),
##	  function(object, value){
##		  linearParamElementReplace(object, "mad.AB", value)
##	  })
##setReplaceMethod("mad.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
##	  function(object, value){
##		  linearParamElementReplace(object, "mad.BB", value)
##	  })

setReplaceMethod("nuA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "nuA", value)
	  })

setReplaceMethod("nuB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "nuB", value)
})

setReplaceMethod("phiA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "phiA", value)
})

setReplaceMethod("phiB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "phiB", value)
})

setReplaceMethod("tau2A", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "tau2A", value)
})

setReplaceMethod("tau2B", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "tau2B", value)
})



setReplaceMethod("flags", signature=signature(object="CNSet", value="ff_or_matrix"),
		 function(object, value){
			 linearParamElementReplace(object, "flags", value)
})
setReplaceMethod("flags", signature=signature(object="CNSet", value="ff_matrix"),
		 function(object, value){
			 linearParamElementReplace(object, "flags", value)
})

setReplaceMethod("snpCall", "CNSet",
                 function(object, value){
			 assayDataElementReplace(object, "call", value)
		 })

setReplaceMethod("snpCallProbability", "CNSet",
                 function(object, value){
			 assayDataElementReplace(object, "callProbability", value)
		 })
