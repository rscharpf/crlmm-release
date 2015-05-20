## Method("initialize", "AlleleSet",
##        function(.Object,
##                 assayData = assayDataNew(alleleA=alleleA,
##                                          alleleB=alleleB, ...),
##                 phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
##                 featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
##                 experimentData = new("MIAME"),
##                 annotation = character(),
##                 protocolData = phenoData[,integer(0)],
##                 alleleA = new("matrix"),
##                 alleleB = matrix(numeric(),
## 		                    nrow=nrow(alleleA), ncol=ncol(alleleA),
##                                  dimnames=dimnames(alleleA)),
## 		   chromosome=integer(),
## 		   position=integer(),
## 		   isSnp=integer(),
##                 ...) {
## 		  .Object <- callNextMethod(.Object,
## 					    assayData = assayData,
## 					    phenoData = phenoData,
## 					    featureData = featureData,
## 					    experimentData = experimentData,
## 					    annotation = annotation,
## 					    protocolData = protocolData)
## 		  if(length(annotation) < 1){
## 			  if((length(position) < 1 | length(chromosome) < 1| length(isSnp) < 1)){
## 				  stop("must specify annotation if 'chromosome', 'position', and 'isSnp' are missing")
## 			  } else {
## 				  pData(featureData)$chromosome <- chromosome
## 				  pData(featureData)$position <- position
## 				  pData(featureData)$isSnp <- isSnp
## 			  }
## 		  } else{
## 			  .Object@annotation <- annotation
## 			  if((length(position) < 1 | length(chromosome) < 1| length(isSnp) < 1)){
## 				  if(!isSupportedAnnotation(annotation)){
## 					  stop("The annotation is not supported. Arguments 'chromosome', 'position', and 'isSnp' can be omitted from the initialization only if the annotation is supported (see oligoClasses:::supportedAnnotation()).")
## 				  }
## 			  } else {
## 				  pData(featureData)$chromosome <- chromosome
## 				  pData(featureData)$position <- position
## 				  pData(featureData)$isSnp <- isSnp
## 			  }
## 			  .Object@featureData <- featureData
## 		  }
## 		  ## Do after annotation has been assigned
## 		  if(!(all(c("chromosome", "position", "isSnp") %in% varLabels(featureData))) & isSupportedAnnotation(annotation)){
## 			  ##update the featureData
## 			  .Object@featureData <- addFeatureAnnotation.crlmm(.Object)
## 		  }
## 		  .Object
##        })
## 
## ow to make the initialization platform-specific?
## Method("initialize", "SnpSuperSet",
##        function(.Object,
## 		   call=new("matrix"),		   
##                 callProbability=matrix(NA, nrow(call), ncol(call), dimnames=dimnames(call)),
##                 phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),		   
## 		   featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
## 		   experimentData=new("MIAME"),
## 		   annotation=character(),
## 		   protocolData=phenoData[, integer(0)],
## 		   position=integer(),
## 		   chromosome=integer(),
## 		   isSnp=integer(),...){
## 		  .Object <- callNextMethod(.Object,
## 					    call=call,
## 					    callProbability=callProbability,
## 					    phenoData=phenoData,
## 					    featureData=featureData,
## 					    experimentData=experimentData,
## 					    annotation=annotation,
## 					    protocolData=protocolData,
## 					    position=position,
## 					    chromosome=chromosome,
## 					    isSnp=isSnp, ...)
##        })
##

##setMethod("initialize", "SnpSuperSet",
##          function(.Object,
##		   call=new("matrix"),		   
##                   callProbability=matrix(NA, nrow(call), ncol(call), dimnames=dimnames(call)),
##                   alleleA = new("matrix"),
##                   alleleB = matrix(numeric(),
##		                    nrow=nrow(alleleA), ncol=ncol(alleleA),
##                                    dimnames=dimnames(alleleA)),		   
##                   phenoData = annotatedDataFrameFrom(call, byrow=FALSE),		   
##		   featureData=annotatedDataFrameFrom(call, byrow=TRUE),
##		   experimentData=new("MIAME"),
##		   protocolData=phenoData[, integer(0)],
##		   position=integer(),
##		   chromosome=integer(),
##		   isSnp=integer(),
##		   annotation=character(), ... ){
##		  ##browser()
##		  ##the ... should be additional assayDataElements, if any
##		  .Object <- callNextMethod(.Object,
##					    call=call,
##					    callProbability=callProbability,
##					    alleleA=alleleA,
##					    alleleB=alleleB,
##					    phenoData=phenoData,
##					    featureData=featureData,
##					    experimentData=experimentData,
##					    protocolData=protocolData,
##					    annotation=annotation, ...)
##		  annotation <- .Object@annotation
##		  ##add chromosome, position, isSnp to featureData
##		  if(length(annotation) < 1){
##			  if((length(position) < 1| length(chromosome) < 1 | length(isSnp) < 1)){
##				  stop("must specify annotation if 'chromosome', 'position', and 'isSnp' are missing")
##			  } else {
##				  pData(featureData)$chromosome <- chromosome
##				  pData(featureData)$position <- position
##				  pData(featureData)$isSnp <- isSnp
##			  }
##		  } else{
##			  if((length(position) < 1| length(chromosome) < 1 | length(isSnp) < 1)){
##				  if(!isSupportedAnnotation(annotation)){
##					  stop("The annotation is not supported. Arguments 'chromosome', 'position', and 'isSnp' can be omitted from the initialization only if the annotation is supported (see oligoClasses:::supportedAnnotation()).")
##				  }
##			  } else {
##				  pData(featureData)$chromosome <- chromosome
##				  pData(featureData)$position <- position
##				  pData(featureData)$isSnp <- isSnp
##			  }
##			  .Object@featureData <- featureData
##		  }
##		  ##Do after annotation has been assigned
##		  if(!(all(c("chromosome", "position", "isSnp") %in% varLabels(featureData))) & isSupportedAnnotation(annotation)){
##			  .Object@featureData <- addFeatureAnnotation.crlmm(.Object)
##		  }		  
##		  .Object
##          })



##setMethod("addFeatureAnnotation", "SnpSuperSet", function(object, ...){
##	addFeatureAnnotation.crlmm(object, ...)
##})

##getParam.SnpSuperSet <- function(object, name, batch){
##		  label <- paste(name, batch, sep="_")
##		  colindex <- grep(label, fvarLabels(object))
##		  if(length(colindex) == 1){
##			  param <- fData(object)[, colindex]
##		  }
##		  if(length(colindex) < 1){
##			  param <- NULL
##		  }
##		  if(is.na(colindex)){
##			  stop(paste(label, " not found in object"))
##		  }
##		  if(length(colindex) > 1){
##			  stop(paste(label, " not unique"))
##		  }
##		  return(param)
##	  }



##setMethod("splitByChromosome", "SnpSuperSet", function(object, cnOptions){
##	tmpdir <- cnOptions[["tmpdir"]]
##	outdir <- cnOptions[["outdir"]]	
##	save.it <- cnOptions[["save.it"]]
##	path <- system.file("extdata", package=paste(annotation(object), "Crlmm", sep=""))	
##	load(file.path(path, "snpProbes.rda"))
##	snpProbes <- get("snpProbes")
##	load(file.path(path, "cnProbes.rda"))
##	cnProbes <- get("cnProbes")	
##	k <- grep("chr", colnames(snpProbes))
##	if(length(k) < 1) stop("chr or chromosome not in colnames(snpProbes)")
##	for(CHR in 1:24){
##		cat("Chromosome ", CHR, "\n")
##		snps <- rownames(snpProbes)[snpProbes[, k] == CHR]
##		cnps <- rownames(cnProbes)[cnProbes[, k] == CHR]
##		index <- c(match(snps, featureNames(object)),
##			   match(cnps, featureNames(object)))
##		index <- index[!is.na(index)]
##		callSetPlus <- object[index, ]
##		if(CHR != 24){
##			cnSet <- computeCopynumber(callSetPlus, cnOptions)
##			
##		} else{
##			message("Copy number estimates not available for chromosome Y.  Saving only the 'callSetPlus' object for this chromosome")
##			save(callSetPlus, file=file.path(outdir, paste("callSetPlus_", CHR, ".rda", sep="")))
##		}
##		if(cnOptions[["hiddenMarkovModel"]] & CHR != 24){
##			cnSet <- computeHmm(cnSet, cnOptions)
##		}
##		save(cnSet, file=file.path(outdir, paste("cnSet_", CHR, ".rda", sep="")))
##		saved.objects <- list.files(outdir, pattern="cnSet", full.names=TRUE)
####		} else{ ## save crlmmSet to outdir
####			save(cnSet, file=file.path(outdir, paste("cnSet_", CHR, ".rda", sep="")))
####			saved.objects <- list.files(outdir, pattern="cnSet", full.names=TRUE)			
####		}		
##	}
##	saved.objects
##})

##setMethod("computeCopynumber", "SnpSuperSet",
##	  function(object, cnOptions){
##		  computeCopynumber.SnpSuperSet(object, cnOptions)
##	  })

	
