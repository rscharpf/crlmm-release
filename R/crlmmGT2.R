crlmmGT2 <- function(A, B, SNR, mixtureParams, cdfName, row.names=NULL,
                     col.names=NULL, probs=c(1/3, 1/3, 1/3), DF=6,
                     SNRMin=5, recallMin=10, recallRegMin=1000,
                     gender=NULL, desctrucitve=FALSE, verbose=TRUE,
                     returnParams=FALSE, badSNP=.7,
		     callsGt,
		     callsPr){
	pkgname <- getCrlmmAnnotationName(cdfName)
	stopifnot(require(pkgname, character.only=TRUE, quietly=!verbose))
	open(SNR)
	open(A)
	open(B)
	open(mixtureParams)
	## expect objects to be ff
	keepIndex <- which( SNR[] > SNRMin)
	if(length(keepIndex)==0) stop("No arrays above quality threshold!")
	loader("preprocStuff.rda", .crlmmPkgEnv, pkgname)
	gns <- getVarInEnv("gns", .crlmmPkgEnv)
	if(is.null(rownames(A))){
		stopifnot(nrow(A) == length(gns))
		index <- seq_len(nrow(A))
	} else {
		index <- match(gns, rownames(A))
	}
	snpBatches <- splitIndicesByLength(index, ocProbesets())
	NR <- length(unlist(snpBatches))
	if(verbose) message("Calling ", NR, " SNPs for recalibration... ")
	NC <- ncol(A)
	##
	if(verbose) message("Loading annotations.")
	obj1 <- loader("genotypeStuff.rda", .crlmmPkgEnv, pkgname)
	obj2 <- loader("mixtureStuff.rda", .crlmmPkgEnv, pkgname)
	##
	## this is toget rid of the 'no visible binding' notes
	## variable definitions
	XIndex <- getVarInEnv("XIndex")
	autosomeIndex <- getVarInEnv("autosomeIndex")
	YIndex <- getVarInEnv("YIndex")
	SMEDIAN <- getVarInEnv("SMEDIAN")
	theKnots <- getVarInEnv("theKnots")
	regionInfo <- getVarInEnv("regionInfo")
	params <- getVarInEnv("params")
	rm(list=c(obj1, obj2), envir=.crlmmPkgEnv)
	rm(obj1, obj2)
	## use lexical scope
	imputeGender <- function(XIndex, YIndex){
		if(length(YIndex) > 0){
			a <- log2(as.matrix(A[XIndex,,drop=FALSE]))
			b <- log2(as.matrix(B[XIndex,,drop=FALSE]))
			meds.X <- (apply(a+b, 2, median))/2
			## Y
			a <- log2(as.matrix(A[YIndex,,drop=FALSE]))
			b <- log2(as.matrix(B[YIndex,,drop=FALSE]))
			meds.Y <- (apply(a+b, 2, median))/2
			R <- meds.X - meds.Y
			if(sum(SNR[] > SNRMin) == 1){
				gender <- ifelse(R[SNR[] > SNRMin] > 0.5, 2L, 1L)
			} else{
				gender <- kmeans(R, c(min(R[SNR[]>SNRMin]), max(R[SNR[]>SNRMin])))[["cluster"]]
			}
		} else {
			XMedian <- apply(log2(as.matrix(A[XIndex,,drop=FALSE]))+log2(as.matrix(B[XIndex,, drop=FALSE])), 2, median)/2
			if(sum(SNR > SNRMin) == 1){
				gender <- which.min(c(abs(XMedian-8.9), abs(XMedian-9.5)))
			} else{
				gender <- kmeans(XMedian, c(min(XMedian[SNR[]>SNRMin]), max(XMedian[SNR[]>SNRMin])))[["cluster"]]
			}
		}
		return(gender)
	}
	if(is.null(gender)){
		if(ocProbesets() < length(XIndex)){
			if(verbose) message("Using ", ocProbesets(), " SNPs on chrom X and Y to assign gender.")
			XIndex2 <- sample(XIndex, ocProbesets(), replace=FALSE)
		} else XIndex2 <- XIndex
		if(ocProbesets() < length(YIndex)){
			YIndex2 <- sample(YIndex, ocProbesets(), replace=FALSE)
		} else YIndex2 <- YIndex
		message("Imputing gender")
		gender <- imputeGender(XIndex=XIndex2, YIndex=YIndex2)
		##cnSet$gender[,] <- gender
	}
	Indexes <- list(autosomeIndex, XIndex, YIndex)
	cIndexes <- list(keepIndex,
			 keepIndex[which(gender[keepIndex]==2)],
			 keepIndex[which(gender[keepIndex]==1)])
	## call C
	fIndex <- which(gender==2)
	mIndex <- which(gender==1)
	## different here
	## use gtypeCallerR in batches
	batchSize <- ocProbesets()
	open(A)
	open(B)
	open(mixtureParams)
	process1 <- function(idxBatch){
		snps <- snpBatches[[idxBatch]]
		##rSnps <- range(snps)
		last <- (idxBatch-1)*batchSize
		IndexesBatch <- list(autosomeIndex[autosomeIndex %in% snps]-last,
				     XIndex[XIndex %in% snps]-last,
				     YIndex[YIndex %in% snps]-last)
		IndexesBatch <- lapply(IndexesBatch, as.integer)
		tmpA <- as.matrix(A[snps,])
		tmpB <- as.matrix(B[snps,])
		tmp <- gtypeCallerR(tmpA, tmpB, fIndex, mIndex,
				    params[["centers"]][snps,],
				    params[["scales"]][snps,],
				    params[["N"]][snps,],
				    IndexesBatch, cIndexes,
				    sapply(IndexesBatch, length),
				    sapply(cIndexes, length), SMEDIAN,
				    theKnots, mixtureParams[], DF, probs, 0.025)
		tmp
	}
	## Lexical scope
	gc(verbose=FALSE)
	newparamsBatch <- ocLapply(seq(along=snpBatches), process1, neededPkgs="crlmm")
	gc(verbose=FALSE)
	if(verbose) message("finished process1")
	newparams <- vector("list", 3)
	names(newparams) <- c("centers", "scales", "N")
	newparams[["centers"]] <- do.call("rbind", lapply(newparamsBatch, "[[", 1))
	newparams[["scales"]] <- do.call("rbind", lapply(newparamsBatch, "[[", 2))
	newparams[["N"]] <- do.call("rbind", lapply(newparamsBatch, "[[", 3))
	rm(newparamsBatch); gc(verbose=FALSE)
	if(verbose) message("Done.")
	if(verbose) message("Estimating recalibration parameters.")
	d <- newparams[["centers"]] - params$centers
	##
	##regression
	Index <- intersect(which(pmin(newparams[["N"]][, 1],
				      newparams[["N"]][, 2],
				      newparams[["N"]][, 3]) > recallMin &
				 !apply(regionInfo, 1, any)),
			   autosomeIndex)
	if(length(Index) < recallRegMin){
		warning("Recalibration not possible. Possible cause: small sample size.")
		newparams <- params
		dev <- vector("numeric", nrow(newparams[["centers"]]))
		SS <- matrix(Inf, 3, 3)
		DD <- 0
	}else{
		data4reg <- as.data.frame(newparams[["centers"]][Index,])
		names(data4reg) <- c("AA", "AB", "BB")
		regParams <- cbind(  coef(lm(AA~AB*BB, data=data4reg)),
				   c(coef(lm(AB~AA+BB, data=data4reg)), 0),
				   coef(lm(BB~AA*AB, data=data4reg)))
		rownames(regParams) <- c("intercept", "X", "Y", "XY")
		rm(data4reg)
		##
		minN <- 3
		newparams[["centers"]][newparams[["N"]] < minN] <- NA
		Index <- setdiff(which(rowSums(is.na(newparams[["centers"]]))==1), YIndex)
		if(verbose) message("Filling out empty centers", appendLF=FALSE)
		for(i in Index){
			if(verbose) if(i%%10000==0) message(".", appendLF=FALSE)
			mu <- newparams[["centers"]][i, ]
			j <- which(is.na(mu))
			newparams[["centers"]][i, j] <- c(1, mu[-j], prod(mu[-j]))%*%regParams[, j]
			rm(mu, j)
		}
		##
		##remaing NAs are made like originals
		if(length(YIndex)>0){
			noMoveIndex <- union(setdiff(which(rowSums(is.na(newparams[["centers"]]))>0), YIndex),
					     YIndex[rowSums(is.na(newparams[["centers"]][YIndex, ])>1)])
		}
		snps2ignore <- which(rowSums(is.na(newparams[["centers"]])) > 0)
		snps2keep <- setdiff(autosomeIndex, snps2ignore)
		rm(snps2ignore)
		newparams[["centers"]][is.na(newparams[["centers"]])] <- params[["centers"]][is.na(newparams[["centers"]])]
		if(verbose) cat("\n")
		##
		if(verbose) message("Calculating and standardizing size of shift... ", appendLF=FALSE)
		GG <- DD <- newparams[["centers"]] - params[["centers"]]
		DD <- sweep(DD, 2, colMeans(DD[autosomeIndex, ]))
		SS <- cov(DD[autosomeIndex, ])
		SSI <- solve(SS)
		dev <- vector("numeric", nrow(DD))
		if(length(YIndex)){
			dev[-YIndex] <- apply(DD[-YIndex, ], 1, function(x) x%*%SSI%*%x)
			dev[-YIndex] <- 1/sqrt( (2*pi)^3*det(SS))*exp(-0.5*dev[-YIndex])
			##Now Y (only two params)
			SSY <- SS[c(1, 3), c(1, 3)]
			SSI <- solve(SSY)
			dev[YIndex] <- apply(DD[YIndex, c(1, 3)], 1, function(x) x%*%SSI%*%x)
			dev[YIndex] <- 1/sqrt( (2*pi)^2*det(SSY))*exp(-0.5*dev[YIndex])
		} else {
			dev=apply(DD,1,function(x) x%*%SSI%*%x)
			dev=1/sqrt( (2*pi)^3*det(SS))*exp(-0.5*dev)
		}
		gc(verbose=FALSE)
	}
	if (verbose) message("OK")
	##
	## BC: must keep SD
	params[-2] <- newparams[-2]
	rm(newparams)
	gc(verbose=FALSE)
	##
	if(verbose) message("Calling ", NR, " SNPs... ", appendLF=FALSE)
	##
	## ###################
	## ## MOVE TO C#######
	##
	## running in batches
	callsGt.present <- !missing(callsGt)
	callsPr.present <- !missing(callsPr)
	overwriteAB <- !callsGt.present & !callsPr.present
	if(!overwriteAB){
		open(callsGt)
		open(callsPr)
	}
	process2 <- function(idxBatch){
		snps <- snpBatches[[idxBatch]]
		tmpA <- as.matrix(A[snps,])
		tmpB <- as.matrix(B[snps,])
		rSnps <- range(snps)
		last <- (idxBatch-1)*batchSize
		IndexesBatch <- list(autosomeIndex[autosomeIndex %in% snps]-last,
				     XIndex[XIndex %in% snps]-last,
				     YIndex[YIndex %in% snps]-last)
		IndexesBatch <- lapply(IndexesBatch, as.integer)
		ImNull <- gtypeCallerR2(tmpA, tmpB, fIndex, mIndex,
					params[["centers"]][snps,],
					params[["scales"]][snps,],
					params[["N"]][snps,],
					IndexesBatch, cIndexes,
					sapply(IndexesBatch, length),
					sapply(cIndexes, length),
					SMEDIAN, theKnots, mixtureParams[],
					DF, probs, 0.025,
					which(regionInfo[snps, 2]),
					which(regionInfo[snps, 1]))
		if(overwriteAB){
			A[snps,] <- tmpA
			B[snps,] <- tmpB
		} else {
			callsGt[snps, ] <- tmpA
			callsPr[snps, ] <- tmpB
		}
	}
	gc(verbose=FALSE)
	ocLapply(seq(along=snpBatches), process2, neededPkgs="crlmm")
	close(A)
	close(B)
	close(mixtureParams)
	if(!overwriteAB){
		close(callsGt)
		close(callsPr)
	}
	gc(verbose=FALSE)
	message("Done with process2")
	##  END MOVE TO C#######
	## ##################
	##
	dev <- dev/(dev+1/383)
	if(!is.null(row.names)){ rownames(A) <- rownames(B) <- names(dev) <- row.names}
	if(!is.null(col.names)){ colnames(A) <- colnames(B) <- col.names}
	##
	if(length(Index) >= recallRegMin){
		tmp4batchQC <- DD[autosomeIndex,]/(params[["N"]][autosomeIndex,]+1)
		tmpSnpQc <- dev[autosomeIndex]
		SS <- cov(tmp4batchQC[tmpSnpQc < badSNP,])
		batchQC <- mean(diag(SS))
	}else{
		SS <- matrix(0, 3, 3)
		batchQC <- Inf
	}
	##
	if(verbose) message("Done.")
	if (returnParams){
		return(list(calls=A, confs=B, SNPQC=dev, batchQC=batchQC, params=params, DD=DD, covDD=SS, gender=gender, pkgname=pkgname))
	}else{
		return(list(calls=A, confs=B, SNPQC=dev, batchQC=batchQC, DD=DD, covDD=SS, gender=gender, pkgname=pkgname))
	}
}
