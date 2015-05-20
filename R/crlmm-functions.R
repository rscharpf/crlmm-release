crlmm <- function(filenames, row.names=TRUE, col.names=TRUE,
                  probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5, gender=NULL,
                  save.it=FALSE, load.it=FALSE, intensityFile,
                  mixtureSampleSize=10^5, eps=0.1, verbose=TRUE,
                  cdfName, sns, recallMin=10, recallRegMin=1000,
                  returnParams=FALSE, badSNP=.7){
  if ((load.it | save.it) & missing(intensityFile))
	  stop("'intensityFile' is missing, and you chose either load.it or save.it")
  if (missing(sns)) sns <- basename(filenames)
  if (any(duplicated(sns))) stop("sample identifiers are not unique")
  if (!missing(intensityFile))
    if (load.it & !file.exists(intensityFile)){
      load.it <- FALSE
      message("File ", intensityFile, " does not exist.")
      message("Not loading it, but running SNPRMA from scratch.")
    }
  if (!load.it){
    res <- snprma(filenames, fitMixture=TRUE,
                  mixtureSampleSize=mixtureSampleSize, verbose=verbose,
                  eps=eps, cdfName=cdfName, sns=sns)
    if(save.it){
      t0 <- proc.time()
      save(res, file=intensityFile)
      t0 <- proc.time()-t0
      if (verbose) message("Used ", t0[3], " seconds to save ", intensityFile, ".")
    }
  }else{
    if (verbose) message("Loading ", intensityFile, ".")
    obj <- load(intensityFile)
    if (verbose) message("Done.")
    if (obj != "res")
      stop("Object in ", intensityFile, " seems to be invalid.")
  }
  if(row.names) row.names=res$gns else row.names=NULL
  if(col.names) col.names=res$sns else col.names=NULL
  res2 <- crlmmGT(res[["A"]], res[["B"]], res[["SNR"]],
                  res[["mixtureParams"]], res[["cdfName"]],
                  gender=gender, row.names=row.names,
                  col.names=col.names, recallMin=recallMin,
                  recallRegMin=1000, SNRMin=SNRMin,
                  returnParams=returnParams, badSNP=badSNP,
                  verbose=verbose)

  res2[["SNR"]] <- res[["SNR"]]
  res2[["SKW"]] <- res[["SKW"]]
  return(list2SnpSet(res2, returnParams=returnParams))
}

crlmmGT <- function(A, B, SNR, mixtureParams, cdfName, row.names=NULL,
                    col.names=NULL, probs=c(1/3, 1/3, 1/3), DF=6,
                    SNRMin=5, recallMin=10, recallRegMin=1000,
                    gender=NULL, desctrucitve=FALSE, verbose=TRUE,
                    returnParams=FALSE, badSNP=.7){

  keepIndex <- which(SNR>SNRMin)
  if(length(keepIndex)==0) stop("No arrays above quality threshold!")

  NC <- ncol(A)
  NR <- nrow(A)

  pkgname <- getCrlmmAnnotationName(cdfName)
  if(!require(pkgname, character.only=TRUE, quietly=!verbose)){
    suggCall <- paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep="")
    msg <- paste("If", pkgname,
                 "is installed on an alternative location, please load it manually by using",
                 suggCall)
    message(strwrap(msg))
    stop("Package ", pkgname, " could not be found.")
    rm(suggCall, msg)
  }

  if(verbose) message("Loading annotations.")
  loader("genotypeStuff.rda", .crlmmPkgEnv, pkgname)
  loader("mixtureStuff.rda", .crlmmPkgEnv, pkgname)

  ## this is toget rid of the 'no visible binding' notes
  ## variable definitions
  XIndex <- getVarInEnv("XIndex")
  autosomeIndex <- getVarInEnv("autosomeIndex")
  YIndex <- getVarInEnv("YIndex")
  SMEDIAN <- getVarInEnv("SMEDIAN")
  theKnots <- getVarInEnv("theKnots")
  regionInfo <- getVarInEnv("regionInfo")
  params <- getVarInEnv("params")

  if(is.null(gender)){
	  if(verbose) message("Determining gender.")
	  gender <- imputeGender(A, B, XIndex, YIndex, SNR, SNRMin)
  }

  Indexes <- list(autosomeIndex, XIndex, YIndex)
  cIndexes <- list(keepIndex,
                   keepIndex[which(gender[keepIndex]==2)],
                   keepIndex[which(gender[keepIndex]==1)])

  if(verbose) cat("Calling", NR, "SNPs for recalibration... ")

  ## call C
  fIndex <- which(gender==2)
  mIndex <- which(gender==1)
  newparams <- gtypeCallerR(A, B, fIndex, mIndex,
                            params[["centers"]], params[["scales"]], params[["N"]],
                            Indexes, cIndexes,
                            sapply(Indexes, length), sapply(cIndexes, length),
                            SMEDIAN, theKnots,
                            mixtureParams, DF, probs, 0.025)
  gc(verbose=FALSE)
  names(newparams) <- c("centers", "scales", "N")

  if(verbose) message("Done.")
  if(verbose) message("Estimating recalibration parameters.")
  d <- newparams[["centers"]] - params$centers

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

    minN <- 3
    newparams[["centers"]][newparams[["N"]] < minN] <- NA
    Index <- setdiff(which(rowSums(is.na(newparams[["centers"]]))==1), YIndex)
    if(verbose) cat("Filling out empty centers")
    for(i in Index){
      if(verbose) if(i%%10000==0)cat(".")
      mu <- newparams[["centers"]][i, ]
      j <- which(is.na(mu))
      newparams[["centers"]][i, j] <- c(1, mu[-j], prod(mu[-j]))%*%regParams[, j]
    }

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

    if(verbose) message("Calculating and standardizing size of shift.")
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
  }

  ## BC: must keep SD
  params[-2] <- newparams[-2]

  rm(newparams);gc(verbose=FALSE)
  if(verbose) cat("Calling", NR, "SNPs... ")
  ## ###################
  ## ## MOVE TO C#######
  ImNull <- gtypeCallerR2(A, B, fIndex, mIndex, params[["centers"]],
                          params[["scales"]], params[["N"]], Indexes,
                          cIndexes, sapply(Indexes, length),
                          sapply(cIndexes, length), SMEDIAN, theKnots,
                          mixtureParams, DF, probs, 0.025,
                          which(regionInfo[,2]),
                          which(regionInfo[,1]))
  gc(verbose=FALSE)
  ##  END MOVE TO C#######
  ## ##################

  dev <- dev/(dev+1/383)
  if(!is.null(row.names)){ rownames(A) <- rownames(B) <- names(dev) <- row.names}
  if(!is.null(col.names)){ colnames(A) <- colnames(B) <- col.names}

  if(length(Index) >= recallRegMin){
   tmp4batchQC <- DD[autosomeIndex,]/(params[["N"]][autosomeIndex,]+1)
   tmpSnpQc <- dev[autosomeIndex]
   SS <- cov(tmp4batchQC[tmpSnpQc < badSNP,])
   batchQC <- mean(diag(SS))
  }else{
    SS <- matrix(0, 3, 3)
    batchQC <- Inf
  }

  if(verbose) message("Done.")
  if (returnParams){
    return(list(calls=A, confs=B, SNPQC=dev, batchQC=batchQC, params=params, DD=DD, covDD=SS, gender=gender, pkgname=pkgname))
  }else{
    return(list(calls=A, confs=B, SNPQC=dev, batchQC=batchQC, DD=DD, covDD=SS, gender=gender, pkgname=pkgname))
  }
}


gtypeCallerR <- function(A, B, fIndex, mIndex, theCenters, theScales,
                         theNs, Indexes, cIndexes, nIndexes,
                         ncIndexes, SMEDIAN, knots, params, dft,
                         probs, trim){

  stopifnot(!missing(A), !missing(B), dim(A)==dim(B),
            nrow(A)==nrow(theCenters), nrow(A)==nrow(theScales),
            nrow(A) == nrow(theNs), length(knots)==3, nrow(params)==4,
            ncol(params)==ncol(A), length(trim)==1, length(probs)==3)

  ## make code robust
  ## check types before passing to C

  .Call("gtypeCallerPart1", A, B, as.integer(fIndex),
        as.integer(mIndex), as.numeric(theCenters),
        as.numeric(theScales), as.integer(theNs),
        lapply(Indexes, as.integer), lapply(cIndexes, as.integer),
        as.integer(nIndexes), as.integer(ncIndexes),
        as.numeric(SMEDIAN), as.numeric(knots),
        as.numeric(params), as.integer(dft), as.numeric(probs),
        as.numeric(trim), PACKAGE="crlmm")

}

gtypeCallerR2 <- function(A, B, fIndex, mIndex, theCenters, theScales,
                         theNs, Indexes, cIndexes, nIndexes,
                         ncIndexes, SMEDIAN, knots, params, dft,
                         probs, trim, noTraining, noInfo){

  stopifnot(!missing(A), !missing(B), dim(A)==dim(B),
            nrow(A)==nrow(theCenters), nrow(A)==nrow(theScales),
            nrow(A) == nrow(theNs), length(knots)==3, nrow(params)==4,
            ncol(params)==ncol(A), length(trim)==1, length(probs)==3)

  .Call("gtypeCallerPart2", A, B,
        as.integer(fIndex), as.integer(mIndex),
        as.numeric(theCenters), as.numeric(theScales),
        as.integer(theNs), Indexes, cIndexes, nIndexes, ncIndexes,
        as.numeric(SMEDIAN), as.numeric(knots), as.numeric(params),
        as.integer(dft), as.numeric(probs), as.numeric(trim),
        as.integer(noTraining), as.integer(noInfo), PACKAGE="crlmm")

}

### parallel version
crlmm2 <- function(filenames, row.names=TRUE, col.names=TRUE,
                   probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5, gender=NULL,
                   save.it=FALSE, load.it=FALSE, intensityFile,
                   mixtureSampleSize=10^5, eps=0.1, verbose=TRUE,
                   cdfName, sns, recallMin=10, recallRegMin=1000,
                   returnParams=FALSE, badSNP=.7){
  if ((load.it || save.it) && missing(intensityFile))
    stop("'intensityFile' is missing, and you chose either load.it or save.it")
  if (missing(sns)) sns <- basename(filenames)
  if (!missing(intensityFile))
    if (load.it & !file.exists(intensityFile)){
      load.it <- FALSE
      message("File ", intensityFile, " does not exist.")
      message("Not loading it, but running SNPRMA from scratch.")
    }
  if (!load.it){
    res <- snprma2(filenames, fitMixture=TRUE,
                   mixtureSampleSize=mixtureSampleSize, verbose=verbose,
                   eps=eps, cdfName=cdfName, sns=sns)
    open(res[["A"]])
    open(res[["B"]])
    open(res[["SNR"]])
    open(res[["mixtureParams"]])
    if(save.it){
      t0 <- proc.time()
      save(res, file=intensityFile)
      t0 <- proc.time()-t0
      if (verbose) message("Used ", t0[3], " seconds to save ", intensityFile, ".")
    }
  }else{
    if (verbose) message("Loading ", intensityFile, ".")
    obj <- load(intensityFile)
    if (verbose) message("Done.")
    if (obj != "res")
      stop("Object in ", intensityFile, " seems to be invalid.")
  }
  if(row.names) row.names=res$gns else row.names=NULL
  if(col.names) col.names=res$sns else col.names=NULL
  res2 <- crlmmGT2(res[["A"]], res[["B"]], res[["SNR"]],
                   res[["mixtureParams"]], res[["cdfName"]],
                   gender=gender, row.names=row.names,
                   col.names=col.names, recallMin=recallMin,
                   recallRegMin=1000, SNRMin=SNRMin,
                   returnParams=returnParams, badSNP=badSNP,
                   verbose=verbose)

  res2[["SNR"]] <- res[["SNR"]]
  res2[["SKW"]] <- res[["SKW"]]
  return(list2SnpSet(res2, returnParams=returnParams))
}

imputeGender <- function(A, B, XIndex, YIndex, SNR, SNRMin){
	if(length(YIndex) > 0){
		a <- log2(as.matrix(A[XIndex,,drop=FALSE]))
		b <- log2(as.matrix(B[XIndex,,drop=FALSE]))
		meds.X <- (apply(a+b, 2, median))/2
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
