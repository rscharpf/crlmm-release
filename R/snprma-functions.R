snprma <- function(filenames, mixtureSampleSize=10^5, fitMixture=FALSE,
                   eps=0.1, verbose=TRUE, seed=1, cdfName, sns){
  if (missing(sns)) sns <- basename(filenames)
  if (missing(cdfName))
    cdfName <- read.celfile.header(filenames[1])$cdfName
  pkgname <- getCrlmmAnnotationName(cdfName)

  if(!require(pkgname, character.only=TRUE, quietly=!verbose)){
    suggCall <- paste("library(", pkgname,
                      ", lib.loc='/Altern/Lib/Loc')",
                      sep="")
    msg <- paste("If", pkgname,
                 "is installed on an alternative location,",
                 "please load it manually by using", suggCall)
    message(strwrap(msg))
    stop("Package ", pkgname, " could not be found.")
  }

  if(verbose) message("Loading annotations and mixture model parameters.")
  loader("preprocStuff.rda", .crlmmPkgEnv, pkgname)
  loader("genotypeStuff.rda", .crlmmPkgEnv, pkgname)
  loader("mixtureStuff.rda", .crlmmPkgEnv, pkgname)
  autosomeIndex <- getVarInEnv("autosomeIndex")
  pnsa <- getVarInEnv("pnsa")
  pnsb <- getVarInEnv("pnsb")
  fid <- getVarInEnv("fid")
  reference <- getVarInEnv("reference")
  aIndex <- getVarInEnv("aIndex")
  bIndex <- getVarInEnv("bIndex")
  SMEDIAN <- getVarInEnv("SMEDIAN")
  theKnots <- getVarInEnv("theKnots")
  gns <- getVarInEnv("gns")

  ##We will read each cel file, summarize, and run EM one by one
  ##We will save parameters of EM to use later
  mixtureParams <- matrix(0, 4, length(filenames))
  SNR <- vector("numeric", length(filenames))
  SKW <- vector("numeric", length(filenames))
##   mixtureParams <- initializeBigMatrix("crlmmMixt-", 4, length(filenames))
##   SNR <- initializeBigVector("crlmmSNR-", length(filenames), "numeric")
##   SKW <- initializeBigVector("crlmmSKW-", length(filenames), "numeric")

  ## This is the sample for the fitting of splines
  ## BC: I like better the idea of the user passing the seed,
  ##     because this might intefere with other analyses
  ##     (like what happened to GCRMA)
  set.seed(seed)
  idx <- sort(sample(autosomeIndex, mixtureSampleSize))
  ##S will hold (A+B)/2 and M will hold A-B
  ##NOTE: We actually dont need to save S. Only for pics etc...
  ##f is the correction. we save to avoid recomputing
  A <- matrix(as.integer(0), length(pnsa), length(filenames))
  B <- matrix(as.integer(0), length(pnsb), length(filenames))

  if(verbose){
    message("Processing ", length(filenames), " files.")
    pb <- txtProgressBar(min=0, max=length(filenames), style=3)
  }

  ##for skewness. no need to do everything
  idx2 <- sample(length(fid), 10^5)

  ##We start looping throug cel files
  for(i in seq(along=filenames)){
    y <- as.matrix(read.celfile(filenames[i], intensity.means.only=TRUE)[["INTENSITY"]][["MEAN"]][fid])
    x <- log2(y[idx2])
    SKW[i] <- mean((x-mean(x))^3)/(sd(x)^3)
    rm(x)
    y <- normalize.quantiles.use.target(y, target=reference)
    A[, i] <- intMedianSummaries(y[aIndex, 1, drop=FALSE], pnsa)
    B[, i] <- intMedianSummaries(y[bIndex, 1, drop=FALSE], pnsb)
    ##Now to fit the EM
    if(fitMixture){
      S <- (log2(A[idx, i])+log2(B[idx, i]))/2 - SMEDIAN
      M <- log2(A[idx, i])-log2(B[idx, i])

      ##we need to test the choice of eps.. it is not the max diff between funcs
      tmp <- fitAffySnpMixture56(S, M, theKnots, eps=eps)

      mixtureParams[, i] <- tmp[["coef"]]
      SNR[i] <- tmp[["medF1"]]^2/(tmp[["sigma1"]]^2+tmp[["sigma2"]]^2)
    }
    if (verbose) setTxtProgressBar(pb, i)
  }
  if (verbose) close(pb)
  if (!fitMixture) SNR <- mixtureParams <- NA
  ## gns comes from preprocStuff.rda
  list(A=A, B=B, sns=sns, gns=gns, SNR=SNR, SKW=SKW, mixtureParams=mixtureParams, cdfName=cdfName)
}

fitAffySnpMixture56 <- function(S, M, knots, probs=rep(1/3, 3), eps=.01, maxit=10, verbose=FALSE){
  ##56 stands for 5 and 6 arrays but will also work for Illumina
  ##Note the unfortunate choice of numbering:
  ##1 is BB, 2 AB, and 3 AA. Opposite to everything else!
  ##this is legacy code I decided not to change.
  ## this why at the end we report -coefs: F1 is the negative f
  mus <- append(quantile(M, c(1, 5)/6, names=FALSE), 0, 1)
  sigmas <- rep(mad(c(M[M<mus[1]]-mus[1], M[M>mus[3]]-mus[3])), 3)
  sigmas[2] <- sigmas[2]/2

  weights <- apply(cbind(mus, sigmas), 1, function(p) dnorm(M, p[1], p[2]))
  previousF1 <- -Inf
  change <- eps+1
  it <- 0

  if(verbose) message("Max change must be under ", eps, ".")
  matS <- stupidSplineBasis(S, knots)
  while (change > eps & it < maxit){
    it <- it+1
    ## E
    z <- sweep(weights, 2, probs, "*")
    LogLik <- rowSums(z)
    z <- sweep(z, 1, LogLik, "/")
    probs <- colMeans(z)

    ## M
    fit1 <- crossprod(chol2inv(chol(crossprod(sweep(matS, 1, z[, 1], FUN="*"), matS))), crossprod(matS, z[, 1]*M))

    fit2 <- sum(z[, 2]*M)/sum(z[, 2])
    F1 <- matS%*%fit1
    sigmas[c(1, 3)] <- sqrt(sum(z[, 1]*(M-F1)^2)/sum(z[, 1]))
    sigmas[2] <- sqrt(sum(z[, 2]*(M-fit2)^2)/sum(z[, 2]))

    weights[, 1] <- dnorm(M, F1, sigmas[1])
    weights[, 2] <- dnorm(M, fit2, sigmas[2])
    weights[, 3] <- dnorm(M, -F1, sigmas[3])

    change <- max(abs(F1-previousF1))
    previousF1 <- F1
    if(verbose) message("Iter ", it, ": ", change, ".")
  }
  medF1 <- median(-F1)
 return(list(coef= -fit1, medF1=medF1, sigma1=sigmas[1], sigma2=sigmas[2]))
}

snprma2 <- function(filenames, mixtureSampleSize=10^5, fitMixture=FALSE,
                    eps=0.1, verbose=TRUE, seed=1, cdfName, sns){
  if (missing(sns)) sns <- basename(filenames)
  if (missing(cdfName))
    cdfName <- read.celfile.header(filenames[1])[["cdfName"]]
  pkgname <- getCrlmmAnnotationName(cdfName)
  stopifnot(require(pkgname, character.only=TRUE, quietly=!verbose))

  if(verbose) message("Loading annotations and mixture model parameters.")
  obj <- loader("preprocStuff.rda", .crlmmPkgEnv, pkgname)
  pnsa <- getVarInEnv("pnsa")
  pnsb <- getVarInEnv("pnsb")
  gns <- getVarInEnv("gns")
  rm(list=obj, envir=.crlmmPkgEnv)
  rm(obj)

  ##We will read each cel file, summarize, and run EM one by one
  ##We will save parameters of EM to use later
  if(verbose) message("Initializing objects.")
  mixtureParams <- initializeBigMatrix("crlmmMixt-", 4, length(filenames), "double")
  SNR <- initializeBigVector("crlmmSNR-", length(filenames), "double")
  SKW <- initializeBigVector("crlmmSKW-", length(filenames), "double")

  ## This is the sample for the fitting of splines
  ## BC: I like better the idea of the user passing the seed,
  ##     because this might intefere with other analyses
  ##     (like what happened to GCRMA)
  ##S will hold (A+B)/2 and M will hold A-B
  ##NOTE: We actually dont need to save S. Only for pics etc...
  ##f is the correction. we save to avoid recomputing
  A <- initializeBigMatrix("crlmmA-", length(pnsa), length(filenames), "integer")
  B <- initializeBigMatrix("crlmmB-", length(pnsb), length(filenames), "integer")

  sampleBatches <- splitIndicesByNode(seq(along=filenames))

  if(verbose) message("Processing ", length(filenames), " files.")

  ocLapply(sampleBatches, processCEL, filenames=filenames,
           fitMixture=fitMixture, A=A, B=B, SKW=SKW, SNR=SNR,
           mixtureParams=mixtureParams, eps=eps, seed=seed,
           mixtureSampleSize=mixtureSampleSize, pkgname=pkgname,
           neededPkgs=c("crlmm", pkgname))
  close(mixtureParams)
  close(SNR)
  close(SKW)
  close(A)
  close(B)

  list(A=A, B=B, sns=sns, gns=gns, SNR=SNR, SKW=SKW, mixtureParams=mixtureParams, cdfName=cdfName)
}


processCEL <- function(i, filenames, fitMixture, A, B, SKW, SNR,
                       mixtureParams, eps, seed, mixtureSampleSize,
                       pkgname){
  obj1 <- loader("preprocStuff.rda", .crlmmPkgEnv, pkgname)
  obj2 <- loader("genotypeStuff.rda", .crlmmPkgEnv, pkgname)
  obj3 <- loader("mixtureStuff.rda", .crlmmPkgEnv, pkgname)
  autosomeIndex <- getVarInEnv("autosomeIndex")
  pnsa <- getVarInEnv("pnsa")
  pnsb <- getVarInEnv("pnsb")
  fid <- getVarInEnv("fid")
  reference <- getVarInEnv("reference")
  aIndex <- getVarInEnv("aIndex")
  bIndex <- getVarInEnv("bIndex")
  SMEDIAN <- getVarInEnv("SMEDIAN")
  theKnots <- getVarInEnv("theKnots")
  gns <- getVarInEnv("gns")
  rm(list=c(obj1, obj2, obj3), envir=.crlmmPkgEnv)
  rm(obj1, obj2, obj3)

  ## for mixture
  set.seed(seed)
  idx <- sort(sample(autosomeIndex, mixtureSampleSize))
  ##for skewness. no need to do everything
  idx2 <- sample(length(fid), 10^5)

  open(A)
  open(B)
  open(SKW)
  open(mixtureParams)
  open(SNR)

  for (k in i){
    y <- as.matrix(read.celfile(filenames[k], intensity.means.only=TRUE)[["INTENSITY"]][["MEAN"]][fid])
    x <- log2(y[idx2])
    SKW[k] <- mean((x-mean(x))^3)/(sd(x)^3)
    rm(x)
    y <- normalize.quantiles.use.target(y, target=reference)
    A[, k] <- intMedianSummaries(y[aIndex, 1, drop=FALSE], pnsa)
    B[, k] <- intMedianSummaries(y[bIndex, 1, drop=FALSE], pnsb)
    rm(y)

    if(fitMixture){
      S <- (log2(A[idx,k])+log2(B[idx, k]))/2 - SMEDIAN
      M <- log2(A[idx, k])-log2(B[idx, k])
      tmp <- fitAffySnpMixture56(S, M, theKnots, eps=eps)
      rm(S, M)
      mixtureParams[, k] <- tmp[["coef"]]
      SNR[k] <- tmp[["medF1"]]^2/(tmp[["sigma1"]]^2+tmp[["sigma2"]]^2)
      rm(tmp)
    } else {
      mixtureParams[, k] <- NA
      SNR[k] <- NA
    }
  }
  close(A)
  close(B)
  close(SKW)
  close(mixtureParams)
  close(SNR)
  rm(list=ls())
  gc(verbose=FALSE)
  TRUE
}
