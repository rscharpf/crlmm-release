setMethod("posteriorMean", signature(object="CNSet"), function(object) assayDataElement(object, "posteriorMean"))
setReplaceMethod("posteriorMean", signature(object="CNSet", value="matrix"), function(object, value) assayDataElementReplace(object, "posteriorMean", value))
##setReplaceMethod("mixtureParams", signature(object="CNSet", value="ANY"), function(object, value) object@mixtureParams <- mixtureParams)

linearParamElementReplace <- function(obj, elt, value) {
    storage.mode <- storageMode(batchStatistics(obj))
    switch(storage.mode,
           "lockedEnvironment" = {
               aData <- copyEnv(batchStatistics(obj))
               if (is.null(value)) rm(list=elt, envir=aData)
               else aData[[elt]] <- value
               Biobase:::assayDataEnvLock(aData)
               batchStatistics(obj) <- aData
           },
           "environment" = {
               if (is.null(value)) rm(list=elt, envir=batchStatistics(obj))
               else batchStatistics(obj)[[elt]] <- value
           },
           list = batchStatistics(obj)[[elt]] <- value)
    obj
}

## parameters
## allele A
##   autosome SNPs
##   autosome NPs
##   chromosome X NPs for women
C1 <- function(object, marker.index, batch.index, sample.index){
	acn <- matrix(NA, nrow=length(marker.index), ncol=length(sample.index))
	for(k in seq_along(batch.index)){
		l <- batch.index[k]
		## calculate cn for all the samples that are in this batch
		jj <- sample.index[as.character(batch(object))[sample.index] == batchNames(object)[l]]
		bg <- nuA(object)[marker.index, l]
		slope <- phiA(object)[marker.index, l]
		I <- as.matrix(A(object)[marker.index, jj])
		acn[, match(jj, sample.index)] <- 1/slope * (I - bg)
	}
	return(as.matrix(acn))
}

## allele B  (treated allele 'A' for chromosome X NPs)
##   autosome SNPs
##   chromosome X for male nonpolymorphic markers
C2 <- function(object, marker.index, batch.index, sample.index, NP.X=FALSE){
	acn <- matrix(NA, nrow=length(marker.index), ncol=length(sample.index))
	for(k in seq_along(batch.index)){
		l <- batch.index[k]
		jj <- sample.index[as.character(batch(object))[sample.index] == batchNames(object)[l]]
		bg <- nuB(object)[marker.index, l]
		slope <- phiB(object)[marker.index, l]
		if(!NP.X){
			I <- as.matrix(B(object)[marker.index, jj])
		} else I <- as.matrix(A(object)[marker.index, jj])
		acn[, match(jj, sample.index)] <- 1/slope * (I - bg)
	}
##	if(length(acn) > 1){
##		acn <- do.call("cbind", acn)
##	} else acn <- acn[[1]]
	return(as.matrix(acn))
}

## Chromosome X SNPs
C3 <- function(object, allele, marker.index, batch.index, sample.index){
##	acn <- vector("list", length(batch.index))
	acn <- matrix(NA, nrow=length(marker.index), ncol=length(sample.index))
	for(k in seq_along(batch.index)){
		l <- batch.index[k]
		##j <- which(as.character(batch(object))[sample.index] == batchNames(object)[l])
		jj <- sample.index[as.character(batch(object))[sample.index] == batchNames(object)[l]]
		##phiA2 and phiB2 are not always estimable  -- need both men and women
		phiA2 <- phiPrimeA(object)[marker.index, l]
		phiB2 <- phiPrimeB(object)[marker.index, l]
		phiA <- phiA(object)[marker.index, l]
		phiB <- phiB(object)[marker.index, l]
		nuA <- nuA(object)[marker.index, l]
		nuB <- nuB(object)[marker.index, l]
		phiA <- phiA(object)[marker.index, l]
		IA <- as.matrix(A(object)[marker.index, jj])
		IB <- as.matrix(B(object)[marker.index, jj])
		## I = nu + acn * phi
		## acn = 1/phi* (I-nu)
		phistar <- phiB2/phiA
		tmp <- (IB - nuB - phistar*IA + phistar*nuA)/phiB
		CB <- tmp/(1-phistar*phiA2/phiB)
		index <- which(is.na(phiB2) | is.na(phiA2))
		if(length(index) > 0){
			cb <- 1/phiB[index] * (IB[index, ] - nuB[index])
			CB[index, ] <- cb
		}
		if(allele == "B"){
			acn[, match(jj, sample.index)] <- CB
			##acn[[k]] <- CB
		}
		if(allele == "A"){
			CA <- (IA-nuA-phiA2*CB)/phiA
			if(length(index) > 0){
				ca <- 1/phiA[index] * (IA[index, ] - nuA[index])
 				CA[index, ] <- ca
			}
			acn[, match(jj, sample.index)] <- CA
		}
	}
##	if(length(acn) > 1){
##		acn <- do.call("cbind", acn)
##	} else acn <- acn[[1]]
	return(as.matrix(acn))
}




ACN <- function(object, allele, i , j){
	if(missing(i) & missing(j)) stop("must specify rows (i) or columns (j)")
	is.ff <- is(calls(object), "ff") | is(calls(object), "ffdf")
	missing.i <- missing(i)
	missing.j <- missing(j)
	if(!missing.i){
		is.ann <- !is.na(chromosome(object)[i])
		is.X <- chromosome(object)[i]==23 & is.ann
		is.auto <- chromosome(object)[i] < 23 & is.ann
		is.snp <- isSnp(object)[i] & is.ann
	} else{
		is.ann <- !is.na(chromosome(object))
		is.X <- chromosome(object)==23 & is.ann
		is.auto <- chromosome(object) < 23 & is.ann
		is.snp <- isSnp(object) & is.ann
		i <- 1:nrow(object)
	}
	## Define batch.index and sample.index
	if(!missing.j) {
		batches <- unique(as.character(batch(object)[j]))
		##batches <- as.character(batch(object)[j])
		batch.index <- match(batches, batchNames(object))
	} else {
		batch.index <- seq_along(batchNames(object))
		j <- 1:ncol(object)
	}
	nr <- length(i)
	nc <- length(j)
	acn <- matrix(NA, nr, nc)
	dimnames(acn) <- list(featureNames(object)[i],
			      sampleNames(object)[j])
	if(allele == "A"){
		if(is.ff){
			open(nuA(object))
			open(phiA(object))
			open(A(object))
		}
		## --
		## 4 types of markers for allele A
		##--
		## 1. autosomal SNPs or autosomal NPs
		if(any(is.auto)){
			auto.index <- which(is.auto)
			marker.index <- i[is.auto]
			acn[auto.index, ] <- C1(object, marker.index, batch.index, j)
		}
		if(any(is.X)){
			##2. CHR X SNPs (men and women)
			if(any(is.snp)){
				if(is.ff) {
					open(phiPrimeA(object))
					open(phiPrimeB(object))
					open(phiB(object))
					open(nuB(object))
					open(B(object))
				}
				marker.index <- i[is.X & is.snp]
				acn.index <- which(is.X & is.snp)
				acn[acn.index, ] <- C3(object, allele="A", marker.index, batch.index, j)
				if(is.ff) {
					close(phiPrimeA(object))
					close(phiPrimeB(object))
					close(phiB(object))
					close(nuB(object))
					close(B(object))
				}
			}
			if(any(!is.snp)){
				marker.index <- i[is.X & !is.snp]
				acn.index <- which(is.X & !is.snp)
				acn[acn.index, ] <- NA
				female.index <- j[object$gender[j] == 2]
				## 3. CHR X NPs: women
				if(length(female.index) > 0){
					female.batch.index <- match(unique(as.character(batch(object))[female.index]), batchNames(object))
					jj <- which(object$gender[j] == 2)
					acn[acn.index, jj] <- C1(object, marker.index, female.batch.index, female.index)
				}
				male.index <- j[object$gender[j] == 1]
				if(length(male.index) > 0){
					if(is.ff){
						open(nuB(object))
						open(phiB(object))
					}
					male.batch.index <- match(unique(as.character(batch(object))[male.index]), batchNames(object))
					jj <- which(object$gender[j] == 1)
					acn[acn.index, jj] <- C2(object, marker.index, male.batch.index, male.index, NP.X=TRUE)
					if(is.ff){
						close(nuB(object))
						close(phiB(object))
					}
				}
			}
		}
		if(is.ff){
			close(nuA(object))
			close(phiA(object))
			close(A(object))
		}
	}
	if(allele == "B"){
		if(is.ff){
			open(nuB(object))
			open(phiB(object))
			open(B(object))
		}
		if(any(!is.snp)){
			acn.index <- which(!is.snp)
			acn[acn.index, ] <- 0
		}
		if(any(is.auto)){
			auto.index <- which(is.auto & is.snp)
			if(length(auto.index) > 0){
				marker.index <- i[auto.index]
				acn[auto.index, ] <- C2(object, marker.index, batch.index, j)
			}
		}
		if(any(is.X)){
			if(is.ff){
				open(phiPrimeA(object))
				open(phiPrimeB(object))
				open(phiA(object))
				open(nuA(object))
				open(A(object))
			}
			marker.index <- i[is.X & is.snp]
			acn.index <- which(is.X & is.snp)
			acn[acn.index, ] <- C3(object, allele="B", marker.index, batch.index, j)
			if(is.ff){
				close(phiPrimeA(object))
				close(phiPrimeB(object))
				close(phiA(object))
				close(nuA(object))
				close(A(object))
			}
			if(any(!is.snp)){
				acn.index <- which(!is.snp)
				marker.index <- i[!is.snp]
				acn[acn.index, ] <- 0
			}
		}
	}
	if(is.ff){
		close(nuB(object))
		close(phiB(object))
		close(B(object))
	}
	return(acn)
}

setMethod("CA", signature=signature(object="CNSet"),
	  function(object, ...){
		  ca <- ACN(object, allele="A", ...)
		  ca[ca < 0] <- 0
		  ca[ca > 5] <- 5
		  return(ca)
	  })
setMethod("CB", signature=signature(object="CNSet"),
	  function(object, ...) {
		  cb <- ACN(object, allele="B", ...)
		  cb[cb < 0] <- 0
		  cb[cb > 5] <- 5
		  return(cb)
	  })

setMethod("totalCopynumber", signature=signature(object="CNSet"),
	  function(object, ...){
		  ca <- CA(object, ...)
		  cb <- CB(object, ...)
		  return(ca+cb)
	  })
rawCopynumber <- totalCopynumber

setMethod("posteriorProbability", signature(object="CNSet"),
	  function(object, predictRegion, copyNumber=0:4, w){
		  if(missing(w)) w <- rep(1/length(copyNumber),length(copyNumber))
		  stopifnot(sum(w)==1)
		  logI <- array(NA, dim=c(nrow(object), ncol(object), 2), dimnames=list(NULL, NULL, LETTERS[1:2]))
		  getIntensity <- function(object){
			  logI[, , 1] <- log2(A(object))
			  logI[, , 2] <- log2(B(object))
			  return(logI)
		  }
		  logI <- getIntensity(object)
		  ##gts <- lapply(as.list(0:4), genotypes)
		  prob <- array(NA, dim=c(nrow(object), ncol(object), length(copyNumber)),
				dimnames=list(NULL,NULL, paste("copynumber",copyNumber, sep="")))
		  bns <- batchNames(object)
		  snp.index <- which(isSnp(object))
		  if(length(snp.index) > 0){
			  for(i in seq_along(copyNumber)){
				  G <- genotypes(copyNumber[i])
				  P <- array(NA, dim=c(length(snp.index), ncol(object), length(G)),
					     dimnames=list(NULL,NULL,G))
				  for(g in seq_along(G)){
					  gt <- G[g]
					  for(j in seq_along(bns)){
						  this.batch <- bns[j]
						  sample.index <- which(batch(object) == this.batch)
						  mu <- predictRegion[[gt]]$mu[snp.index, , this.batch, drop=FALSE]
						  dim(mu) <- dim(mu)[1:2]
						  if(all(is.na(mu))){
							  P[, sample.index, gt] <- NA
						  } else {
							  Sigma <- predictRegion[[gt]]$cov[snp.index, , this.batch, drop=FALSE]
							  dim(Sigma) <- dim(Sigma)[1:2]
							  P[, sample.index, gt] <- dbvn(x=logI[snp.index, sample.index, , drop=FALSE], mu=mu, Sigma=Sigma)
						  }
					  }
				  }
				  ## the marginal probability for the total copy number
				  ##  -- integrate out the probability for the different genotypes
				  for(j in 1:ncol(P)){
					  PP <- P[, j, , drop=FALSE]
					  dim(PP) <- dim(PP)[c(1, 3)]
					  prob[snp.index, j, i] <- rowSums(PP, na.rm=TRUE)
				  }
			  }
		  } ## length(snp.index) > 0
		  ##---------------------------------------------------------------------------
		  ##
		  ##  nonpolymorphic markers
		  ##
		  ##---------------------------------------------------------------------------
		  np.index <- which(!isSnp(object))
		  if(length(np.index) > 0){
			  for(i in seq_along(copyNumber)){
				  G <- genotypes(copyNumber[i], is.snp=FALSE)
				  P <- array(NA, dim=c(length(np.index), ncol(object), length(G)),
					     dimnames=list(NULL,NULL,G))
				  for(g in seq_along(G)){
					  gt <- G[g]
					  for(j in seq_along(bns)){
						  this.batch <- bns[j]
						  sample.index <- which(batch(object) == this.batch)
						  mu <- predictRegion[[gt]]$mu[np.index, 1, this.batch, drop=FALSE]
						  dim(mu) <- dim(mu)[1]
						  if(all(is.na(mu))){
							  P[, sample.index, gt] <- NA
						  } else {
							  Sigma <- predictRegion[[gt]]$cov[np.index, 1, this.batch, drop=FALSE]
							  dim(Sigma) <- dim(Sigma)[1]
							  P[, sample.index, gt] <- dnorm(logI[np.index, sample.index, 1], mean=mu, sd=Sigma)
						  }
					  } ## j in seq_along(bns)
				  } ## g in seq_along(G)
				  ## the marginal probability for the total copy number
				  ##  -- integrate out the probability for the different genotypes
				  prob[np.index, , i] <- P
			  } ## seq_along(copyNumber)
		  } ## length(np.index) > 0
		  ##constant.  Probs across copy number states must
		  ##sum to one.  nc <- apply(prob, c(1,3), sum,
		  ##na.rm=TRUE)
		  wm <- matrix(w, nrow(object), length(copyNumber), byrow=TRUE)
		  nc <- matrix(NA, nrow(object), ncol(object))
		  for(j in seq_len(ncol(object))){
			  pp <- prob[, j, ] * wm
			  nc <- rowSums(pp, na.rm=TRUE)
			  prob[, j, ] <- pp/nc
		  }
		  return(prob)
	  })

setMethod("calculatePosteriorMean", signature(object="CNSet"),
	  function(object, posteriorProb, copyNumber=0:4, ...){
		  stopifnot(dim(posteriorProb)[[3]] == length(copyNumber))
		  pm <- matrix(0, nrow(object), ncol(object))
		  for(i in seq_along(copyNumber)){
			  pm <- pm + posteriorProb[, , i] * copyNumber[i]
		  }
		  return(pm)
	  })

##.bivariateCenter <- function(nu, phi){
##			  ##  lexical scope for mus, CA, CB
##			  if(CA <= 2 & CB <= 2 & (CA+CB) < 4){
##				  mus[,1, ] <- log2(nu[, 1, ] + CA *
##						    phi[, 1, ])
##				  mus[,2, ] <- log2(nu[, 2, ] + CB *
##						    phi[, 2, ])
##			  } else { ## CA > 2
##				  if(CA > 2){
##					  theta <- pi/4*Sigma[,2,]
##					  shiftA <- CA/4*phi[, 1, ] * cos(theta)
##					  shiftB <- CA/4*phi[, 1, ] * sin(theta)
##					  mus[, 1, ] <- log2(nu[, 1, ] + 2 * phi[,1,]+shiftA)
##					  mus[, 2, ] <- log2(nu[, 2, ] + CB *phi[,2,]+shiftB)
##				  }
##				  if(CB > 2){
##					  ## CB > 2
##					  theta <- pi/2-pi/4*Sigma[,2,]
##					  shiftA <- CB/4*phi[, 2, ] * cos(theta)
##					  shiftB <- CB/4*phi[, 2, ] * sin(theta)
##					  mus[, 1, ] <- log2(nu[, 1, ] + CA*phi[,1,]+shiftA)
##					  mus[, 2, ] <- log2(nu[, 2, ]+ 2*phi[,2,]+shiftB)
##				  }
##				  if(CA == 2 & CB == 2){
##					  mus[, 1, ] <- log2(nu[, 1, ] + 1/2*CA*phi[,1,])
##					  mus[, 2, ] <- log2(nu[, 2, ]+ 1/2*CB*phi[,2,])
##				  }
##			  }
##			  mus
##		  }

## for a given copy number, return a named list of bivariate normal prediction regions
##   - elements of list are named by genotype
##   'AAA'  list should have 'mu' and 'Sigma'
##                mu is a R x 2 matrix, R is number of features
##                Sigma is a R x 4 matrix.  Each row can be coerced to a 2 x 2 matrix
##          For nonpolymorphic markers, 2nd column of mu is NA and elements 2-4 of Sigma are NA
##
setMethod("predictionRegion", signature(object="CNSet", copyNumber="integer"),
	  function(object, copyNumber=0:4){
		  ## would it be more efficient to store as an array?
		  ## mu: features x genotype x allele
		  ## Sigma: features x genotype x covariance
		  stopifnot(all(copyNumber %in% 0:4))
		  getNu <- function(object){
			  nu[, 1, ] <- nuA(object)
			  nu[, 2, ] <- nuB(object)
			  nu
		  }
		  getPhi <- function(object){
			  phi[,1,] <- phiA(object)
			  phi[,2,] <- phiB(object)
			  phi
		  }
		  gts <- lapply(as.list(copyNumber), genotypes)
		  nms <- unlist(gts)
		  res <- vector("list", length(nms))
		  ##names(res) <- paste("copyNumber", copyNumber, sep="")
		  names(res) <- nms
		  bnms <- batchNames(object)
		  nu <- array(NA, dim=c(nrow(object), 2, length(bnms)))
		  phi <- array(NA, dim=c(nrow(object), 2, length(bnms)))
		  nu <- getNu(object)
		  phi <- getPhi(object)
		  ##mus <- matrix(NA, nrow(nu), 2, dimnames=list(NULL, LETTERS[1:2]))
		  mus <- array(NA, dim=c(nrow(nu), 2, length(bnms)),
			       dimnames=list(NULL, LETTERS[1:2],
			       bnms))
		  ## Sigma <- matrix(NA, nrow(nu), 3)
		  Sigma <- array(NA, dim=c(nrow(nu), 3, length(bnms)),
				 dimnames=list(NULL, c("varA", "cor", "varB"),
				 bnms))
		  bivariateCenter <- function(nu, phi){
			  ##  lexical scope for mus, CA, CB
			  mus[,1, ] <- log2(nu[, 1, ] + CA *
					    phi[, 1, ])
			  mus[,2, ] <- log2(nu[, 2, ] + CB *
					    phi[, 2, ])
			  return(mus)
		  }
		  np.index <- which(!isSnp(object))
		  for(i in seq_along(copyNumber)){
			  G <- genotypes(copyNumber[i])
			  tmp <- vector("list", length(G))
			  names(tmp) <- G
			  CN <- copyNumber[i]
			  for(g in seq_along(G)){
				  gt <- G[g]
				  CB <- g-1
				  CA <- CN-CB
				  gt.corr <- genotypeCorrelation(gt)
				  nma <- ifelse(CA == 0, "tau2A.BB", "tau2A.AA")
				  nmb <- ifelse(CB == 0, "tau2B.AA", "tau2B.BB")
				  Sigma[, 1, ] <- getVar(object, nma)
				  Sigma[, 3, ] <- getVar(object, nmb)
				  Sigma[, 2, ] <- getCor(object, gt.corr)
				  if(length(np.index) > 0){
					  Sigma[, 1, ] <- getVar(object, "tau2A.AA")
					  Sigma[np.index, 2, ] <- NA
				  }
				  res[[gt]]$mu <- bivariateCenter(nu,phi)
				  ## adjust correlation
##				  if(CA == 0 | CB == 0){
##					  Sigma[,2,] <- 0
##				  }
				  res[[gt]]$cov <- Sigma
			  }
		  }
		  res <- as(res, "PredictionRegion")
		  return(res)
	  })



setMethod("xyplotcrlmm", signature(x="formula", data="CNSet", predictRegion="list"),
	  function(x, data, predictRegion, ...){
		  fns <- featureNames(data)
		  fns <- matrix(fns, nrow(data), ncol(data), byrow=FALSE)
		  fns <- as.character(fns)
		  df <- list(A=as.numeric(log2(A(data))),
			     B=as.numeric(log2(B(data))),
			     gt=as.integer(calls(data)),
			     gt.conf=as.numeric(confs(data)),
			     snpid=fns)#, snp=snpId)
		  df <- as.data.frame(df)
		  df$snpid <- factor(df$snpid, levels=unique(df$snpid), ordered=TRUE)
		  bns <- batchNames(data)
		  predictRegion <- lapply(predictRegion, function(x, bns){
			  batch.index <- match(bns, dimnames(x$mu)[[3]])
			  x$mu <- x$mu[, , batch.index, drop=FALSE]
			  x$cov <- x$cov[, , batch.index, drop=FALSE]
			  return(x)
		  }, bns=bns)
		  ##predictRegion is an argument of ABpanel
		  xyplot(x, df, predictRegion=predictRegion, ...)
	  })

setMethod("xyplot", signature(x="formula", data="CNSet"),
	  function(x, data, ...){
		  if("predictRegion" %in% names(list(...))){
			  xyplotcrlmm(x, data, ...)
		  } else{
			  callNextMethod()
		  }
})

setMethod("calculateRBaf", signature(object="CNSet"),
	  function(object, batch.name, chrom){
		  calculateRBafCNSet(object, batch.name, chrom)
	  })

calculateRBafCNSet <- function(object, batch.name, chrom){
	if(missing(batch.name)) {
		batch.name <- batchNames(object)
		if("grandMean" %in% batch.name)
			batch.name <- batch.name[-length(batch.name)]
	}
	if(missing(chrom)) chrom <- unique(chromosome(object))
	if(!(all(batch.name %in% batchNames(object)))) stop("batch.name must be belong to batchNames(object)")
	chr <- chromosome(object)
	valid.chrs <- chr <= 23 & chr %in% chrom
	index <- which(valid.chrs)
	indexlist <- split(index, chr[index])
	J <- which(batch(object) %in% batch.name)
	sns <- sampleNames(object)[J]
	sampleindex <- split(J, factor(batch(object)[J], levels=unique(batch(object)[J])))
	if(!all(valid.chrs)) warning("Only computing log R ratios and BAFs for autosomes and chr X")
	## if ff package is loaded, these will be ff objects
	chr <- names(indexlist)
	rlist <- blist <- vector("list", length(indexlist))
	path <- ldPath()
	processByChromosome <- function(i, chr, path){
		ldPath(path)
		nr <- length(i)
		##CHR <- names(i)
		bafname <- paste("baf_chr", chr, sep="")
		rname <- paste("lrr_chr", chr, sep="")
		bmatrix <- initializeBigMatrix(bafname, nr=nr, nc=length(sns), vmode="integer")
		rmatrix <- initializeBigMatrix(rname, nr=nr, nc=length(sns), vmode="integer")
		colnames(rmatrix) <- colnames(bmatrix) <- sns
		## put rownames in order of physical position
		ix <- order(position(object)[i])
		i <- i[ix]
		rownames(rmatrix) <- rownames(bmatrix) <- featureNames(object)[i]
		for(j in seq_along(sampleindex)){
			bname <- batch.name[j]
			J <- sampleindex[[j]]
			res <- calculateRTheta(object=object, # crlmm:::
						       batch.name=bname,
						       feature.index=i)
			k <- match(sampleNames(object)[J], sns)
			bmatrix[, k] <- res[["baf"]]
			rmatrix[, k] <- res[["lrr"]]
		}
		list(bmatrix, rmatrix)
	}
	## calcualte R BAF by chromosome
	if(isPackageLoaded("ff")){
		pkgs <- c("oligoClasses", "ff", "Biobase", "crlmm")
	} else pkgs <- c("oligoClasses", "Biobase", "crlmm")
	i <- NULL
	res <- foreach(i=indexlist, chr=names(indexlist), .packages=pkgs) %dopar% {
		processByChromosome(i, chr, path)
	}
	blist <- lapply(res, "[[", 1)
	rlist <- lapply(res, "[[", 2)
	res <- list(baf=blist, lrr=rlist)
	return(res)
}

##setMethod("calculateRBaf", signature(object="CNSet"),
##	  function(object, batch.name){
##		  all.autosomes <- all(chromosome(object) < 23)
##		  if(!all.autosomes){
##			  stop("method currently only defined for chromosomes 1-22")
##		  }
##		  if(missing(batch.name)) batch.name <- batchNames(object)[1]
##		  stopifnot(batch.name %in% batchNames(object))
##		  if(length(batch.name) > 1){
##			  warning("only the first batch in batch.name processed")
##			  batch.name <- batch.name[1]
##		  }
##		  RTheta.aa <- calculateRTheta(object, "AA", batch.name)
##		  RTheta.ab <- calculateRTheta(object, "AB", batch.name)
##		  RTheta.bb <- calculateRTheta(object, "BB", batch.name)
##
##		  J <- which(batch(object) == batch.name)
##
##		  theta.aa <- matrix(RTheta.aa[, "theta"], nrow(object), length(J), byrow=FALSE)
##		  theta.ab <- matrix(RTheta.ab[, "theta"], nrow(object), length(J), byrow=FALSE)
##		  theta.bb <- matrix(RTheta.bb[, "theta"], nrow(object), length(J), byrow=FALSE)
##
##		  a <- A(object)[, J, drop=FALSE]
##		  b <- B(object)[, J, drop=FALSE] ## NA's for b where nonpolymorphic
##		  is.np <- !isSnp(object)
##		  ##b[is.np, ] <- a[is.np, ]
##		  b[is.np, ] <- 0L
##		  dns <- dimnames(a)
##		  dimnames(a) <- dimnames(b) <- NULL
##		  obs.theta <- atan2(b, a)*2/pi
##
##		  lessAA <- obs.theta < theta.aa
##		  lessAB <- obs.theta < theta.ab
##		  lessBB <- obs.theta < theta.bb
##		  grAA <- !lessAA
##		  grAB <- !lessAB
##		  grBB <- !lessBB
##		  not.na <- !is.na(theta.aa)
##		  I1 <- grAA & lessAB & not.na
##		  I2 <- grAB & lessBB & not.na
##
##		  bf <- matrix(NA, nrow(object), ncol(a))
##		  bf[I1] <- 0.5 * ((obs.theta-theta.aa)/(theta.ab-theta.aa))[I1]
##		  bf[I2] <- (.5 * (obs.theta - theta.ab) / (theta.bb - theta.ab))[I2] + 0.5
##		  bf[lessAA] <- 0
##		  bf[grBB] <- 1
##
##		  r.expected <- matrix(NA, nrow(object), ncol(a))
##		  r.aa <- matrix(RTheta.aa[, "R"], nrow(object), length(J), byrow=FALSE)
##		  r.ab <- matrix(RTheta.ab[, "R"], nrow(object), length(J), byrow=FALSE)
##		  r.bb <- matrix(RTheta.bb[, "R"], nrow(object), length(J), byrow=FALSE)
##		  rm(RTheta.aa, RTheta.ab, RTheta.bb); gc()
##		  obs.r <- a+b
##
##		  lessAA <- lessAA & not.na
##		  grBB <- grBB & not.na
##		  tmp <- ((obs.theta - theta.aa) * (r.ab-r.aa)/(theta.ab-theta.aa))[I1] + r.aa[I1]
##		  r.expected[I1] <- tmp
##		  tmp <- ((obs.theta - theta.ab) * (r.bb - r.ab)/(theta.bb-theta.ab))[I2] + r.ab[I2]
##		  r.expected[I2] <- tmp
##		  r.expected[lessAA] <- r.aa[lessAA]
##		  r.expected[grBB] <- r.bb[grBB]
##		  index.np <- which(is.np)
##		  if(length(index.np) > 0){
##			  a.np <- A(object)[index.np, J, drop=FALSE]
##			  meds <- rowMedians(a.np, na.rm=TRUE)
##			  r.expected[index.np, ] <- matrix(meds, length(index.np), ncol(a.np))
##		  }
##		  lrr <- log2(obs.r/r.expected)
##
##		  dimnames(bf) <- dimnames(lrr) <- dns
##		  res <- list(baf=bf,
##			      lrr=lrr)
##		  return(res)
##	  })


setAs("CNSet", "oligoSetList", function(from, to){
	constructOligoSetListFrom(from)
})



setMethod(OligoSetList, "CNSet", function(object,...){
	constructOligoSetListFrom(object, ...)
})
setMethod(BafLrrSetList, "CNSet", function(object,...){
  constructBafLrrSetListFrom(object, ...)
})





constructOligoSetListFrom <- function(object, ...){
  ##row.index <- seq_len(nrow(object))
  ##col.index <- seq_len(ncol(object))
  is.lds <- ifelse(is(calls(object), "ff_matrix") | is(calls(object), "ffdf"), TRUE, FALSE)
  if(is.lds) stopifnot(isPackageLoaded("ff"))
  b.r <- calculateRBaf(object, ...)
  b <- b.r[["baf"]]
  r <- b.r[["lrr"]]
  j <- match(colnames(r[[1]]), sampleNames(object))
  rns <- lapply(r, rownames)
  fDList <- foreach(featureid=rns) %do%{
    featureData(object)[match(featureid, featureNames(object)), ]
  }
  names(fDList) <- sapply(fDList, function(x) chromosome(x)[1])
  gtPlist <- gtlist <- vector("list", length(r))
  for(i in seq_along(r)){
    gtlist[[i]] <- initializeBigMatrix("call", nr=nrow(r[[i]]), nc=length(j), vmode="integer")
    gtPlist[[i]] <- initializeBigMatrix("callPr", nr=nrow(r[[i]]), nc=length(j), vmode="integer")
    featureid <- rownames(r[[i]])
    ix <- match(featureid, featureNames(object))
    rownames(gtPlist[[i]]) <- rownames(gtlist[[i]]) <- featureid
    colnames(gtPlist[[i]]) <- colnames(gtlist[[i]]) <- colnames(r[[i]])
    for(k in seq_along(j)){
      gtlist[[i]][, k] <- calls(object)[ix, j[k]]
      gtPlist[[i]][, k] <- snpCallProbability(object)[ix, j[k]]
    }
  }
  ad <- AssayDataList(baf=b, copyNumber=r, call=gtlist, callProbability=gtPlist)
  object <- new("oligoSetList",
                assayDataList=ad,
                featureDataList=fDList,
                chromosome=names(fDList),
                phenoData=phenoData(object)[j, ],
                annotation=annotation(object),
                genome=genomeBuild(object))
  return(object)
}



constructBafLrrSetListFrom <- function(object, ...){
	is.lds <- ifelse(is(calls(object), "ff_matrix") | is(calls(object), "ffdf"), TRUE, FALSE)
	if(is.lds) stopifnot(isPackageLoaded("ff"))
	b.r <- calculateRBaf(object, ...)
	b <- b.r[["baf"]]
	r <- b.r[["lrr"]]
	j <- match(colnames(r[[1]]), sampleNames(object))
	rns <- lapply(r, rownames)
	featureid <- NULL
	fDList <- foreach(featureid=rns) %do%{
		featureData(object)[match(featureid, featureNames(object)), ]
	}
	names(fDList) <- sapply(fDList, function(x) chromosome(x)[1])
	ad <- AssayDataList(baf=b, lrr=r)
	object <- new("BafLrrSetList",
		      assayDataList=ad,
		      featureDataList=fDList,
		      chromosome=names(fDList),
		      phenoData=phenoData(object)[j, ],
		      annotation=annotation(object),
		      genome=genomeBuild(object))
	return(object)
}

computeBR <- constructBafLrrSetListFrom
