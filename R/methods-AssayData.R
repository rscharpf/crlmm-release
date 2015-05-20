
setMethod("Ns", signature(object="AssayData"),
	  function(object, ...){
		  batchnames <- sampleNames(object)
		  N.AA <- assayDataElement(object, "N.AA")
		  N.AB <- assayDataElement(object, "N.AB")
		  N.BB <- assayDataElement(object, "N.BB")
		  res <- array(NA, dim=c(nrow(N.AA), 3, ncol(N.AA)))
		  dimnames(res)[[2]] <- c("AA", "AB", "BB")
		  dimnames(res)[[3]] <- batchnames
		  res[, "AA", ] <- N.AA[,]
		  res[, "AB", ] <- N.AB[,]
		  res[, "BB", ] <- N.BB[,]
		  return(res)
	  })
setMethod("corr", signature(object="AssayData"),
	  function(object, ...){
		  batchnames <- sampleNames(object)
		  corrAA <- assayDataElement(object, "corrAA")
		  corrAB <- assayDataElement(object, "corrAB")
		  corrBB <- assayDataElement(object, "corrBB")
		  res <- array(NA, dim=c(nrow(corrAA), 3, ncol(corrAA)))
		  dimnames(res)[[2]] <- c("AA", "AB", "BB")
		  dimnames(res)[[3]] <- batchnames
		  res[, "AA", ] <- corrAA[,]
		  res[, "AB", ] <- corrAB[,]
		  res[, "BB", ] <- corrBB[,]
		  return(res)
	  })

setMethod("medians", signature(object="AssayData"),
	  function(object, ...){
		  batchnames <- sampleNames(object)
		  medianA.AA <- assayDataElement(object, "medianA.AA")
		  medianA.AB <- assayDataElement(object, "medianA.AB")
		  medianA.BB <- assayDataElement(object, "medianA.BB")
		  medianB.AA <- assayDataElement(object, "medianB.AA")
		  medianB.AB <- assayDataElement(object, "medianB.AB")
		  medianB.BB <- assayDataElement(object, "medianB.BB")
		  res <- array(NA, dim=c(nrow(medianA.AA), 2, 3, ncol(medianA.AA)))
		  dimnames(res)[[2]] <- c("A", "B")
		  dimnames(res)[[3]] <- c("AA", "AB", "BB")
		  dimnames(res)[[4]] <- batchnames
		  res[, "A", "AA", ] <- medianA.AA[,]
		  res[, "A", "AB", ] <- medianA.AB[,]
		  res[, "A", "BB", ] <- medianA.BB[,]
		  res[, "B", "AA", ] <- medianB.AA[,]
		  res[, "B", "AB", ] <- medianB.AB[,]
		  res[, "B", "BB", ] <- medianB.BB[,]
		  return(res)
})

getMedians <- function(object){
	medianA.AA <- assayDataElement(object, "medianA.AA")
	medianA.AB <- assayDataElement(object, "medianA.AB")
	medianA.BB <- assayDataElement(object, "medianA.BB")
	medianB.AA <- assayDataElement(object, "medianB.AA")
	medianB.AB <- assayDataElement(object, "medianB.AB")
	medianB.BB <- assayDataElement(object, "medianB.BB")
	list(A.AA=medianA.AA,
	     A.AB=medianA.AB,
	     A.BB=medianA.BB,
	     B.AA=medianB.AA,
	     B.AB=medianB.AB,
	     B.BB=medianB.BB)
}

setMethod("mads", signature(object="AssayData"),
	  function(object, ...){
		  batchnames <- sampleNames(object)
		  madA.AA <- assayDataElement(object, "madA.AA")
		  madA.AB <- assayDataElement(object, "madA.AB")
		  madA.BB <- assayDataElement(object, "madA.BB")
		  madB.AA <- assayDataElement(object, "madB.AA")
		  madB.AB <- assayDataElement(object, "madB.AB")
		  madB.BB <- assayDataElement(object, "madB.BB")
		  res <- array(NA, dim=c(nrow(madA.AA), 2, 3, ncol(madA.AA)))
		  dimnames(res)[[2]] <- c("A", "B")
		  dimnames(res)[[3]] <- c("AA", "AB", "BB")
		  dimnames(res)[[4]] <- batchnames
		  res[, "A", "AA", ] <- madA.AA[,]
		  res[, "A", "AB", ] <- madA.AB[,]
		  res[, "A", "BB", ] <- madA.BB[,]
		  res[, "B", "AA", ] <- madB.AA[,]
		  res[, "B", "AB", ] <- madB.AB[,]
		  res[, "B", "BB", ] <- madB.BB[,]
		  return(res)
})



setMethod("tau2", signature(object="AssayData"),
	  function(object, ...){
		  batchnames <- sampleNames(object)
		  tau2A.AA <- assayDataElement(object, "tau2A.AA")
		  tau2A.BB <- assayDataElement(object, "tau2A.BB")
		  tau2B.AA <- assayDataElement(object, "tau2B.AA")
		  tau2B.BB <- assayDataElement(object, "tau2B.BB")
		  res <- array(NA, dim=c(nrow(tau2A.AA), 2, 2, ncol(tau2A.AA)))
		  dimnames(res)[[2]] <- c("A", "B")
		  dimnames(res)[[3]] <- c("AA", "BB")
		  dimnames(res)[[4]] <- batchnames
		  res[, "A", "AA", ] <- tau2A.AA[,]
		  res[, "A", "BB", ] <- tau2A.BB[,]
		  res[, "B", "AA", ] <- tau2B.AA[,]
		  res[, "B", "BB", ] <- tau2B.BB[,]
		  return(res)
})


