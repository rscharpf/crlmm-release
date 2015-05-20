# function below works OK provided all .idat files are in the current working directory
# - could add an option to allow files in Illumina directory structure to be handled
# or to use the optional 'Path' column in sampleSheet
# - there is a lot of header information that is currently discarded - could try and store this somewhere in the resulting NChannelSet

readIdatFiles = function(sampleSheet=NULL,
			  arrayNames=NULL,
			  ids=NULL,
			  path=".",
			  arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
			  highDensity=FALSE,
			  sep="_",
			  fileExt=list(green="Grn.idat", red="Red.idat"),
			  saveDate=FALSE, verbose=FALSE) {
	verbose <- FALSE
       if(!is.null(arrayNames)) {
               pd = new("AnnotatedDataFrame", data = data.frame(Sample_ID=arrayNames))
       }
       if(!is.null(sampleSheet)) { # get array info from Illumina's sample sheet
	       if(is.null(arrayNames)){
		       if(!is.null(arrayInfoColNames$barcode) && (arrayInfoColNames$barcode %in% colnames(sampleSheet))) {
			       barcode = sampleSheet[,arrayInfoColNames$barcode]
			       arrayNames=barcode
		       }
		       if(!is.null(arrayInfoColNames$position) && (arrayInfoColNames$position %in% colnames(sampleSheet))) {
			       position = sampleSheet[,arrayInfoColNames$position]
			       if(is.null(arrayNames))
				       arrayNames=position
			       else
				       arrayNames = paste(arrayNames, position, sep=sep)
			       if(highDensity) {
				       hdExt = list(A="R01C01", B="R01C02", C="R02C01", D="R02C02")
				       for(i in names(hdExt))
					       arrayNames = sub(paste(sep, i, sep=""), paste(sep, hdExt[[i]], sep=""), arrayNames)
			       }
		       }
	       }
	       pd = new("AnnotatedDataFrame", data = sampleSheet)
	       sampleNames(pd) = make.unique(basename(arrayNames))
       }
       if(is.null(arrayNames)) {
               arrayNames = gsub(paste(sep, fileExt$green, sep=""), "", dir(pattern=fileExt$green, path=path))
               if(!is.null(sampleSheet)) {
                      sampleSheet=NULL
                      cat("Could not find required info in \'sampleSheet\' - ignoring.  Check \'sampleSheet\' and/or \'arrayInfoColNames\'\n")
               }
               pd = new("AnnotatedDataFrame", data = data.frame(Sample_ID=arrayNames))
       }
       narrays = length(arrayNames)
       grnfiles = paste(arrayNames, fileExt$green, sep=sep)
       redfiles = paste(arrayNames, fileExt$red, sep=sep)
       if(length(grnfiles)==0 || length(redfiles)==0)
	       stop("Cannot find .idat files")
       if(length(grnfiles)!=length(redfiles))
	       stop("Cannot find matching .idat files")
       if(path[1] != "." & path[1] != ""){
	       grnidats = file.path(path, grnfiles)
	       redidats = file.path(path, redfiles)
       }  else {
	       message("path arg not set.  Assuming files are in local directory, or that complete path is provided")
	       grnidats = grnfiles
	       redidats = redfiles
       }
       if(!all(file.exists(grnidats))) stop("Missing some of the *Grn.idat files")
       if(!all(file.exists(redidats))) stop("Missing some of the *Red.idat files")
       headerInfo = list(nProbes = rep(NA, narrays),
                         Barcode = rep(NA, narrays),
                         ChipType = rep(NA, narrays),
                         Manifest = rep(NA, narrays), # not sure about this one - sometimes blank
                         Position = rep(NA, narrays)) # this may also vary a bit
       dates = list(decode=rep(NA, narrays),
                    scan=rep(NA, narrays))
       ## read in the data
       for(i in seq_along(arrayNames)) {
	       if(verbose) {
	       ## RS
		       ##cat("reading", arrayNames[i], "\t")
		       cat("reading", basename(arrayNames[i]), "\t")
		       cat(paste(sep, fileExt$green, sep=""), "\t")
	       }
	       idsG = idsR = G = R = NULL
	       G = readIDAT(grnidats[i])
	       idsG = rownames(G$Quants)
	       headerInfo$nProbes[i] = G$nSNPsRead
	       headerInfo$Barcode[i] = G$Barcode
	       headerInfo$ChipType[i] = G$ChipType
	       headerInfo$Manifest[i] = G$Unknown$MostlyNull
	       headerInfo$Position[i] = G$Unknowns$MostlyA
               if(headerInfo$nProbes[i]>(headerInfo$nProbes[1]+headerInfo$nProbes[1]*0.04) || headerInfo$nProbes[i]<(headerInfo$nProbes[1]-headerInfo$nProbes[1]*0.04)) {
		       warning("Chips are not of the same type.  Skipping ", basename(grnidats[i]), " and ", basename(redidats[i]))
		       next()
	       }
               if(saveDate) {
                      if(nrow(G$RunInfo)>=2) {
	              dates$decode[i] = G$RunInfo[1, 1]
	              dates$scan[i] = G$RunInfo[2, 1]
                      }
               }
	       if(i==1) {
		       if(is.null(ids) && !is.null(G)){
			       ids = idsG
		       } # else stop("Could not find probe IDs")
		       nprobes = length(ids)
		       narrays = length(arrayNames)
		       RG = new("NChannelSet",
		                 R=matrix(0, nprobes, narrays),
		                 G=matrix(0, nprobes, narrays),
		                 zero=matrix(1, nprobes, narrays),
				 annotation=headerInfo$Manifest[1],
				 phenoData=pd, storage.mode="environment")
		       featureNames(RG) = ids
		       if(!is.null(sampleSheet) && !is.null(sampleSheet$Sample_ID)){
		            sampleNames(RG) = make.unique(sampleSheet$Sample_ID)
		       } else  sampleNames(RG) = make.unique(basename(arrayNames))
		       gc(verbose=FALSE)
	       }
	       if(length(ids)==length(idsG)) {
		       if(sum(ids==idsG)==nprobes) {
			       RG@assayData$G[,i] = G$Quants[, "Mean"]
			       zeroG = G$Quants[, "NBeads"]==0
		       }
	       } else {
		       indG = match(ids, idsG)
                       nasG = is.na(indG)
		       RG@assayData$G[!nasG,i] = G$Quants[indG[!nasG], "Mean"]
		       zeroG = G$Quants[indG[!nasG], "NBeads"]==0
	       }
	       rm(G)
	       gc(verbose=FALSE)
	       if(verbose) {
                      cat(paste(sep, fileExt$red, sep=""), "\n")
	       }
	       R = readIDAT(redidats[i])
	       idsR = rownames(R$Quants)

	       if(length(ids)==length(idsG)) {
		       if(sum(ids==idsR)==nprobes) {
			       RG@assayData$R[,i] = R$Quants[ ,"Mean"]
		               zeroR = R$Quants[ ,"NBeads"]==0
                               RG@assayData$zero[,i] = zeroG | zeroR
		       }
	       } else {
		       indR = match(ids, idsR)
	               nasR = is.na(indR)
		       RG@assayData$R[!nasR,i] = R$Quants[indR[!nasR], "Mean"]
		       zeroR = R$Quants[indR[!nasR], "NBeads"]==0
                       RG@assayData$zero[!nasR,i] = zeroG | zeroR
	       }
#	       RG@assayData$zero[,i] = zeroG | zeroR
	       rm(R, zeroG, zeroR)
	       gc(verbose=FALSE)
       }
       if(saveDate) {
	       protocolData(RG)[["ScanDate"]] = dates$scan
       }
       storageMode(RG) = "lockedEnvironment"
       RG
}


getNumberOfSNPs = function(afile, path){
    fullfilename = file.path(path, afile)
    headerSection = readLines(fullfilename, n=15)

    headerLine = headerSection[10][1]
    delimiterList = c(",", "\t")

    headers = unlist(strsplit(headerLine, delimiterList[1]))
    if (length(headers)!=1) {
        delimiterIndex = 1
    }
    if (length(headers) == 1) {
        headers = unlist(strsplit(headerLine, delimiterList[2]))
        if (length(headers) != 1) {
            delimiterIndex = 2
        }
        if (length(headers) == 1) {
            stop("Input file ", fullfilename, " is not delimited by either comm(,) or tab(\\t)")
        }
    }

    SNPLine = headerSection[5][1]
    elements = unlist(strsplit(SNPLine, delimiterList[delimiterIndex]))
    numSNP = as.integer(elements[2])
    return(numSNP)
}

checkNumberOfSNPs = function(filenames, path){
    numSNP = getNumberOfSNPs(filenames[1], path)
    if (length(filenames) > 1) {
		for (i in 2:length(filenames)){
			if (getNumberOfSNPs(filenames[i], path) != numSNP){
				return(FALSE)
			}
		}
    }
    return(TRUE)
}


getNumberOfSamples = function(filenames, path, numSNP){
    sampleCount = rep(0, length(filenames))
    for (i in 1:length(filenames)){
    	# number of sample in input file line 7 is not reliable, calculate from number of lines and number of SNPs
    	fullfilename = file.path(path, filenames[i])
	LineCount = .Call("countFileLines", fullfilename)
  	if (((LineCount - 10) %% numSNP) != 0){
           stop("Please check input file: ", fullfilename, " Line count is not a multiple of number of SNPs")
    	}
	sampleCount[i] = LineCount %/% numSNP
    }	
    return(sampleCount)
}

processOneGenCallFile = function(afile, numSNP, numSample,
    colnames=list("SampleID"="Sample ID", "SNPID"="SNP Name", "XRaw"="X Raw", "YRaw"="Y Raw"),
    verbose=FALSE) {
  
    headerSection = readLines(afile, n=15)

    headerLine = headerSection[10][1]
    delimiterList = c(",", "\t")

    headers = unlist(strsplit(headerLine, delimiterList[1]))
    if (length(headers)!=1) {
        delimiterIndex = 1
    }
    if (length(headers) == 1) {
        headers = unlist(strsplit(headerLine, delimiterList[2]))
        if (length(headers) != 1) {
            delimiterIndex = 2
        }
        if (length(headers) == 1) {
            stop("Input file is not delimited by either comm(,) or tab(\\t)")
        }
    }
     
    if(sum(is.na(match(colnames, headers))) != 0)
	stop("Cannot find required columns: ", colnames[is.na(match(colnames, headers))], " in ", file, "\nPlease check whether this data was exported.")
    
    SNPIDPos = which(headers == colnames$SNPID)
    sampleIDPos = which(headers == colnames$SampleID)
    XValuePos = which(headers == colnames$XRaw)
    YValuePos = which(headers == colnames$YRaw)   
        
    if(verbose) {
        message("Number of SNPs in file: ", afile, " is ", numSNP, " and number of samples is ", numSample)
        if (delimiterIndex == 1) message("File is comma-seperated. ")
        if (delimiterIndex == 2) message("File is tab-seperated. ")              
    }

    values = .Call("readGenCallOutputCFunc", afile, numSNP, numSample, SNPIDPos, sampleIDPos, XValuePos, YValuePos, delimiterIndex)
             
    return(values)
}

readGenCallOutput = function(filenames, path=".", cdfName,
    colnames=list("SampleID"="Sample ID", "SNPID"="SNP Name", "XRaw"="X Raw", "YRaw"="Y Raw"),
    type=list("SampleID"="character", "SNPID"="character", "XRaw"="integer", "YRaw"="integer"), verbose=FALSE) {

    if(!identical(names(type), names(colnames)))
       stop("The arguments 'colnames' and 'type' must have consistent names")
    if(missing(cdfName)) stop("must specify cdfName")

    if (!checkNumberOfSNPs(filenames, path)){
       stop("Number of SNPs in each file must be identical to form one output NChannelSet object.")
    }

    if (verbose) message("Checking number of samples and features in each file. ")
    numSNP = getNumberOfSNPs(filenames[1], path)
    sampleCounts = getNumberOfSamples(filenames, path, numSNP)
    numSample = sum(sampleCounts)

    X = initializeBigMatrix(name = "X", nr = numSNP, nc = numSample, vmode = "integer")
    Y = initializeBigMatrix(name = "Y", nr = numSNP, nc = numSample, vmode = "integer")
    zero = initializeBigMatrix(name = "zero", nr = numSNP, nc = numSample, vmode = "integer")
    
    totSampleNames = rep(NA, numSample)
   
    baseIndex = 1
    if (verbose) message("Start processing ", length(filenames), " input file(s)")
    for (i in 1:length(filenames)){
    	fullfilename = file.path(path, filenames[i])
    	valuesThisFile = processOneGenCallFile(fullfilename, numSNP, sampleCounts[i], colnames, verbose)
	
	if (i == 1){
	    totSNPNames = rownames(valuesThisFile$Xvalues)
	} else {
	    # matching on SNP names? now assume they come in order
	}
	maxIndex = baseIndex + sampleCounts[i] - 1
	X[, baseIndex:maxIndex] = valuesThisFile$Xvalues
	Y[, baseIndex:maxIndex] = valuesThisFile$Yvalues
        zero[, baseIndex:maxIndex] = (X[, baseIndex:maxIndex] == 0) || (Y[, baseIndex:maxIndex] == 0)
	totSampleNames[baseIndex:(baseIndex + sampleCounts[i] - 1)] = colnames(valuesThisFile$Xvalues)
	rm(valuesThisFile)
	baseIndex = baseIndex + sampleCounts[i]
    }

    if(verbose) message("Creating NChannelSet object\n")

    XY = new("NChannelSet", X=X, Y=Y, zero=zero, annotation=cdfName, storage.mode = "environment")
    sampleNames(XY) = totSampleNames
    featureNames(XY) = totSNPNames
    
    if(verbose)
    cat("Done\n")

    XY
}



RGtoXY = function(RG, chipType, verbose=TRUE) {
  needToLoad <- !all(sapply(c('addressA', 'addressB', 'base'), isLoaded))
  if(needToLoad){
	  chipList = c("human1mv1c",# 1M
	  "human370v1c",            # 370CNV
	  "human650v3a",            # 650Y
	  "human610quadv1b",        # 610 quad
	  "human660quadv1a",        # 660 quad
	  "human370quadv3c",        # 370CNV quad
	  "human550v3b",            # 550K
	  "human1mduov3b",          # 1M Duo
	  "humanomni1quadv1b",      # Omni1 quad
	  "humanomni25quadv1b",     # Omni2.5 quad
	  "humanomni258v1a",        # Omni2.5 8 v1 A
          "humanomni258v1p1b",      # Omni2.5 8 v1.1 B
          "humanomni5quadv1b",      # Omni5 quad  
	  "humanomniexpress12v1b",  # Omni express 12
	  "humanimmuno12v1b",       # Immuno chip 12
          "humancytosnp12v2p1h",    # CytoSNP 12
          "humanexome12v1p2a",      # Exome 12 v1.2 A
          "humanomniexpexome8v1p1b") # Omni Express Exome 8 v1.1b
	  ## RS: added cleancdfname()
	  if(missing(chipType)){
		  chipType = match.arg(cleancdfname(annotation(RG)), chipList)
	  } else chipType = match.arg(cleancdfname(chipType), chipList)
	  pkgname = getCrlmmAnnotationName(chipType)
	  if(!require(pkgname, character.only=TRUE, quietly=!verbose)){
		  suggCall = paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep="")
		  msg = paste("If", pkgname, "is installed on an alternative location, please load it manually by using", suggCall)
		  message(strwrap(msg))
		  stop("Package ", pkgname, " could not be found.")
		  rm(suggCall, msg)
	  }
	  if(verbose) message("Loading chip annotation information.")
	  loader("address.rda", .crlmmPkgEnv, pkgname)
  }
  aids = getVarInEnv("addressA") # comes from AddressA_ID or Address column in manifest
  bids = getVarInEnv("addressB") # comes from AddressB_ID or Address2 column in manifest
  snpbase = getVarInEnv("base")
##  loader(‘file.rda’)
##  x = getVarInEnv(‘x’)
##  y = getVarInEnv(‘y’)
##
##  I’d consider using something like:
##
##	  needToLoad = !all(sapply(c(‘x’, ‘y’), isLoaded))
##  if (needToLoad){
##	  loader(‘file.rda’)
##	  x = getVarInEnv(‘x’)
##	  y = getVarInEnv(‘y’)
##  }
  ids = names(aids)
  nsnps = length(aids)
  narrays = ncol(RG)
#  aidcol = match("AddressA_ID", colnames(annot))
#  if(is.na(aidcol))
#    aidcol = match("Address", colnames(annot))
#  bidcol = match("AddressB_ID", colnames(annot))
#  if(is.na(bidcol))
#    bidcol = match("Address2", colnames(annot))
#  aids = annot[, aidcol]
#  bids = annot[, bidcol]
#  snpids = annot[,"Name"]
#  snpbase = annot[,"SNP"]
  infI = !is.na(bids) & bids!=0
  aord = match(aids, featureNames(RG)) # NAs are possible here
  bord = match(bids, featureNames(RG)) # and here
#  argrg = aids[rrgg]
#  brgrg = bids[rrgg]
  xyPhenoData = AnnotatedDataFrame(data=RG@phenoData@data,varMetadata=RG@phenoData@varMetadata)
  levels(xyPhenoData@varMetadata$channel) = c("X","Y","zero","_ALL_")
  XY <- new("NChannelSet",
        assayDataNew(X=matrix(0, nsnps, narrays),Y=matrix(0, nsnps, narrays),zero=matrix(0, nsnps, narrays)),
        phenoData=xyPhenoData, protocolData=RG@protocolData, annotation=chipType)
   storageMode(XY) = "environment"
#  XY <- new("NChannelSet",
#	    X=matrix(0, nsnps, narrays),
#	    Y=matrix(0, nsnps, narrays),
#	    zero=matrix(0, nsnps, narrays),
#	    annotation=chipType, phenoData=RG@phenoData,
#	    protocolData=RG@protocolData, storage.mode="environment")
  featureNames(XY) = ids
  sampleNames(XY) = sampleNames(RG)
  ## RS
  rm(list=c("bord", "infI", "aids", "bids", "ids", "snpbase"))
#  print(gc())
  gc(verbose=FALSE)
  # Need to initialize - matrices filled with NAs to begin with
#  XY@assayData$X[1:nsnps,] = 0
#  XY@assayData$Y[1:nsnps,] = 0
#  XY@assayData$zero[1:nsnps,] = 0

  # First sort out Infinium II SNPs, X -> R (allele A)  and Y -> G (allele B) from the same probe
  ## RS added
  not.na <- !is.na(aord)
  not.na.aord <- aord[not.na]
  ## RS substitution  !is.na(aord) -> not.na
  ##r <- as.matrix(exprs(channel(RG, "R"))[not.na.aord,]) # mostly red
  r <- as.matrix(assayData(RG)[["R"]][not.na.aord, ])
  XY@assayData$X[not.na,] <- r
  rm(r);gc(verbose=FALSE)
  g <- as.matrix(assayData(RG)[["G"]][not.na.aord,]) # mostly green
  XY@assayData$Y[not.na,] <- g
  rm(g); gc(verbose=FALSE)
  ##z <- as.matrix(exprs(channel(RG, "zero"))[not.na.aord,]) # mostly green
  z <- as.matrix(assayData(RG)[["zero"]][not.na.aord, ])
  XY@assayData$zero[not.na,] <- z
  rm(z); gc(verbose=FALSE)
  ##RS added
  rm(RG)
#  print(gc())
  gc(verbose=FALSE)
  ## Warning - not 100% sure that the code below is correct - could be more complicated than this
#  infIRR = infI & (snpbase=="[A/T]" | snpbase=="[T/A]" | snpbase=="[a/t]" | snpbase=="[t/a]")

#  X[infIRR,] = exprs(channel(RG, "R"))[aord[infIRR],] # mostly red
#  Y[infIRR,] = exprs(channel(RG, "R"))[bord[infIRR],] # mostly green

  # Finally Infinium I where X -> G from allele A probe and Y -> G from allele B probe
#  infIGG = infI & (snpbase=="[C/G]" | snpbase=="[G/C]" | snpbase=="[g/c]" | snpbase=="[c/g]")

#  X[infIGG,] = exprs(channel(RG, "G"))[aord[infIGG],] # mostly red
#  Y[infIGG,] = exprs(channel(RG, "G"))[bord[infIGG],] # mostly green

  #  For now zero out Infinium I probes
#  XY@assayData$X[infI,] = 0
#  XY@assayData$Y[infI,] = 0
#  XY@assayData$zero[infI,] = 0
#  gc(verbose=FALSE)
  XY
}


stripNormalize = function(XY, useTarget=TRUE, verbose=TRUE, quantile.method="between") {
        if(quantile.method=="between") {
	  if(useTarget){
		objectsNeeded <- c("stripnum", "reference")
	  } else objectsNeeded <- "stripnum"
	  needToLoad <- !all(sapply(objectsNeeded, isLoaded))
	  if(needToLoad){
		pkgname = getCrlmmAnnotationName(annotation(XY))
		if(!require(pkgname, character.only=TRUE, quietly=!verbose)){
			suggCall = paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep="")
			msg = paste("If", pkgname, "is installed on an alternative location, please load it manually by using", suggCall)
			message(strwrap(msg))
			stop("Package ", pkgname, " could not be found.")
			rm(suggCall, msg)
		}
		if(verbose) message("Loading strip and reference normalization information.")
		loader("preprocStuff.rda", .crlmmPkgEnv, pkgname)
	  }
	  stripnum = getVarInEnv("stripnum")
	  if(useTarget)
		targetdist = getVarInEnv("reference")
	  if(verbose){
		message("Quantile normalizing ", ncol(XY), " arrays by ", max(stripnum), " strips.")
		if (getRversion() > '2.7.0') pb = txtProgressBar(min=0, max=max(stripnum), style=3)
	  }
          for(s in 1:max(stripnum)) {
            if(verbose) {
              if (getRversion() > '2.7.0') setTxtProgressBar(pb, s)
              else cat(".")
            }
            sel = stripnum==s
            ##RS: faster to access data directly
            ##subX = as.matrix(exprs(channel(XY, "X"))[sel,])
            ##subY = as.matrix(exprs(channel(XY, "Y"))[sel,])
            subX <- as.matrix(assayData(XY)[["X"]][sel, ])
            subY <- as.matrix(assayData(XY)[["Y"]][sel, ])
            if(useTarget)
              tmp = normalize.quantiles.use.target(cbind(subX, subY), targetdist[[s]])
            else
              tmp = normalize.quantiles(cbind(subX, subY))
            XY@assayData$X[sel,] = matrix(as.integer(tmp[,1:(ncol(tmp)/2)]+16), nrow(tmp), ncol(tmp)/2)
            XY@assayData$Y[sel,] = matrix(as.integer(tmp[,(ncol(tmp)/2+1):ncol(tmp)]+16), nrow(tmp), ncol(tmp)/2)
            rm(subX, subY, tmp, sel)
            gc(verbose=FALSE)
          }
          if(verbose)
            cat("\n")
        }
        if(quantile.method=="within") {  # ignore strip information
	  if(useTarget){
		objectsNeeded <- c("Xref", "Yref")
	  } else objectsNeeded <- ""
	  needToLoad <- !all(sapply(objectsNeeded, isLoaded))
	  if(needToLoad){
		pkgname = getCrlmmAnnotationName(annotation(XY))
		if(!require(pkgname, character.only=TRUE, quietly=!verbose)){
			suggCall = paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep="")
			msg = paste("If", pkgname, "is installed on an alternative location, please load it manually by using", suggCall)
			message(strwrap(msg))
			stop("Package ", pkgname, " could not be found.")
			rm(suggCall, msg)
		}
		if(verbose) message("Loading reference normalization information.")
		loader("targetXY.rda", .crlmmPkgEnv, pkgname)
	  }
	  if(useTarget) {
                Xref = getVarInEnv("Xref")
	        Yref = getVarInEnv("Yref")
          } else{
                Xref = normalize.quantiles.determine.target(as.matrix(assayData(XY)[["X"]]))
                gc(verbose=FALSE)
                Yref = normalize.quantiles.determine.target(as.matrix(assayData(XY)[["Y"]]))
                gc(verbose=FALSE)
          }
	  if(verbose){
		message("Quantile normalizing ", ncol(XY), " arrays, one at a time.")
		if (getRversion() > '2.7.0') pb = txtProgressBar(min=0, max=ncol(XY), style=3)
	  }
          for(s in 1:ncol(XY)) {
            if(verbose) {
              if (getRversion() > '2.7.0') setTxtProgressBar(pb, s)
              else cat(".")
            }
            ##RS: faster to access data directly
            ##subX = as.matrix(exprs(channel(XY, "X"))[sel,])
            ##subY = as.matrix(exprs(channel(XY, "Y"))[sel,])
            subX <- as.matrix(assayData(XY)[["X"]][,s])
            subY <- as.matrix(assayData(XY)[["Y"]][,s])
            tmpX = normalize.quantiles.use.target(subX, as.integer(Xref))
            tmpY = normalize.quantiles.use.target(subY, as.integer(Yref))              
            XY@assayData$X[,s] = as.integer(tmpX+16)
            XY@assayData$Y[,s] = as.integer(tmpY+16)
            rm(subX, subY, tmpX, tmpY)
          }
          if(verbose)
            cat("\n")
        }
  XY
}


preprocessInfinium2 = function(XY, mixtureSampleSize=10^5,
				fitMixture=TRUE,
				eps=0.1,
				verbose=TRUE,
				seed=1,
				cdfName,
				sns,
				stripNorm=TRUE,
				useTarget=TRUE,
                                quantile.method="between") { #,
#               outdir=".") {
#				save.it=FALSE,
#				snpFile,
#				cnFile) {
	if(stripNorm)
		XY = stripNormalize(XY, useTarget=useTarget, verbose=verbose, quantile.method=quantile.method)
	## MR: the code below is mostly straight from snprma.R
	if (missing(sns)) sns = sampleNames(XY) #$X
	if(missing(cdfName))
		cdfName = annotation(XY)
	needToLoad <- !all(sapply(c("autosomeIndex",
				    "SMEDIAN",
				    "theKnots",
				    "npProbesFid"), isLoaded))
	if(needToLoad){
		pkgname = getCrlmmAnnotationName(cdfName)
		if(!require(pkgname, character.only=TRUE, quietly=!verbose)){
			suggCall = paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep="")
			msg = paste("If", pkgname, "is installed on an alternative location, please load it manually by using", suggCall)
			message(strwrap(msg))
			stop("Package ", pkgname, " could not be found.")
			rm(suggCall, msg)
		}
		if(verbose) message("Loading snp annotation and mixture model parameters.")
		loader("genotypeStuff.rda", .crlmmPkgEnv, pkgname)
                if(fitMixture)
 		  loader("mixtureStuff.rda", .crlmmPkgEnv, pkgname)
		loader("snpProbesFid.rda", .crlmmPkgEnv, pkgname)
		loader("npProbesFid.rda", .crlmmPkgEnv, pkgname)
	}
	autosomeIndex = getVarInEnv("autosomeIndex")
        if(fitMixture) {
          SMEDIAN = getVarInEnv("SMEDIAN")
          theKnots = getVarInEnv("theKnots")
        }
        npIndex = getVarInEnv("npProbesFid")
	snpIndex = getVarInEnv("snpProbesFid")
	narrays = ncol(XY)
	nprobes = length(npIndex)
	if(length(nprobes)>0) {
		## RS: channel creates an expression set.  This is much slower than directly accessing the data
		A <- matrix(as.integer(assayData(XY)[["X"]][npIndex, ]), nprobes, narrays)
		##system.time(A <- matrix(as.integer(exprs(channel(XY, "X"))[npIndex,]), nprobes, narrays))
		B <- matrix(as.integer(assayData(XY)[["Y"]][npIndex, ]), nprobes, narrays)
		##B = matrix(as.integer(exprs(channel(XY, "Y"))[npIndex,]), nprobes, narrays)
		## new lines below - useful to keep track of zeroed out probes
		zero <- matrix(as.integer(assayData(XY)[["zero"]][npIndex, ]), nprobes, narrays)
		##zero = matrix(as.integer(exprs(channel(XY, "zero"))[npIndex,]), nprobes, narrays)
		
		colnames(A) = colnames(B) = colnames(zero) = sns
		rownames(A) = rownames(B) = rownames(zero) = names(npIndex)
		cnAB = list(A=A, B=B, zero=zero, sns=sns, gns=names(npIndex), cdfName=cdfName)
		##      t0 <- proc.time()
		##      save(cnAB, file=cnFile)
		##      t0 <- proc.time()-t0
		##      if(verbose) message("Used ", round(t0[3],1), " seconds to save ", cnFile, ".")
		rm(A, B, zero)
#		print(gc())
		gc(verbose=FALSE)
	}
  # next process snp probes
  nprobes = length(snpIndex)
  ##We will read each cel file, summarize, and run EM one by one
  ##We will save parameters of EM to use later
  mixtureParams = matrix(NA, 4, narrays)
  SNR = rep(NA, narrays)
  SKW = rep(NA, narrays)
  ## This is the sample for the fitting of splines
  ## BC: I like better the idea of the user passing the seed,
  ##     because this might intefere with other analyses
  ##     (like what happened to GCRMA)
  set.seed(seed)
  idx = sort(sample(autosomeIndex, mixtureSampleSize))
  idx2 = sample(nprobes, 10^5)
  ##S will hold (A+B)/2 and M will hold A-B
  ##NOTE: We actually dont need to save S. Only for pics etc...
  ##f is the correction. we save to avoid recomputing

  A = matrix(0, nprobes, narrays)
  B = matrix(0, nprobes, narrays)
  zero = matrix(NA, nprobes, narrays)
  if(verbose && fitMixture){
     message("Calibrating ", narrays, " arrays.")
     if (getRversion() > '2.7.0') pb = txtProgressBar(min=0, max=narrays, style=3)
  }
  for(i in 1:narrays){
	  ## RS: faster to access data directly without using channel method
	  ##A[,i] = as.integer(exprs(channel(XY, "X"))[snpIndex,i])
	  A[, i] <- as.integer(assayData(XY)[["X"]][snpIndex, i])
	  B[, i] <- as.integer(assayData(XY)[["Y"]][snpIndex, i])
	  zero[, i] <- as.integer(assayData(XY)[["zero"]][snpIndex, i])
	  ##B[,i] = as.integer(exprs(channel(XY, "Y"))[snpIndex,i])
	  ##zero[,i] = as.integer(exprs(channel(XY, "zero"))[snpIndex,i])
	  SKW[i] = mean((A[idx2,i]-mean(A[idx2,i]))^3)/(sd(A[idx2,i])^3)
	  if(fitMixture){
		  S = (log2(A[idx,i])+log2(B[idx,i]))/2 - SMEDIAN
		  M = log2(A[idx,i])-log2(B[idx,i])
		  ##we need to test the choice of eps.. it is not the max diff between funcs
		  tmp = fitAffySnpMixture56(S, M, theKnots, eps=eps)
		  mixtureParams[, i] = tmp[["coef"]]
		  SNR[i] = tmp[["medF1"]]^2/(tmp[["sigma1"]]^2+tmp[["sigma2"]]^2)
	          if(verbose) {
		    if (getRversion() > '2.7.0') setTxtProgressBar(pb, i)
		    else cat(".")
                  }
	  }
	  ## run garbage collection every now and then
	  if(i %% 100 == 0) gc(verbose=FALSE);
  }
  if (verbose && fitMixture) {
	  if (getRversion() > '2.7.0') close(pb)
	  else cat("\n")
  }
  if (!fitMixture) SNR = mixtureParams = NA
  ## gns comes from preprocStuff.rda

#  if(save.it & !missing(snpFile)) {
#    t0 <- proc.time()
#    save(res, file=snpFile)
#    t0 <- proc.time()-t0
#    if(verbose) message("Used ", round(t0[3],1), " seconds to save ", snpFile, ".")
#  }

  res = list(A=A, B=B,
             zero=zero, sns=sns, gns=names(snpIndex), SNR=SNR, SKW=SKW,
             mixtureParams=mixtureParams, cdfName=cdfName, cnAB=cnAB)

#  open(res[["A"]])
#  open(res[["B"]])
#  open(res[["zero"]])
#  open(res[["SNR"]])
#  open(res[["mixtureParams"]])

#  if(save.it & !missing(snpFile)) {
#    t0 <- proc.time()
#    save(res, file=snpFile)
#    t0 <- proc.time()-t0
#    if(verbose) message("Used ", round(t0[3],1), " seconds to save ", snpFile, ".")
#  }

#  close(res[["A"]])
#  close(res[["B"]])
#  close(res[["zero"]])
#  close(res[["SNR"]])
#  close(res[["mixtureParams"]])

  return(res)
}


crlmmIllumina <- function(RG, XY, stripNorm=TRUE, useTarget=TRUE,
			  row.names=TRUE, col.names=TRUE,
			  probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5, gender=NULL,
			  seed=1,
			  mixtureSampleSize=10^5, eps=0.1, verbose=TRUE,
			  cdfName, sns, recallMin=10, recallRegMin=1000,
			  returnParams=FALSE, badSNP=.7) {
	if(missing(cdfName)) {
		if(!missing(RG))
			cdfName = annotation(RG)
		if(!missing(XY))
			cdfName = annotation(XY)
	}
	if(!isValidCdfName(cdfName))
		stop("cdfName not valid.  see validCdfNames")
	if(!missing(RG)) {
		if(missing(XY))
			XY = RGtoXY(RG, chipType=cdfName)
		else
			stop("Both RG and XY specified - please use one or the other")
	}
	if (missing(sns)) sns = sampleNames(XY)
	res <- preprocessInfinium2(XY, mixtureSampleSize=mixtureSampleSize,
				   fitMixture=TRUE, verbose=verbose,
				   seed=seed, eps=eps, cdfName=cdfName, sns=sns,
				   stripNorm=stripNorm, useTarget=useTarget) #,
	if(row.names) row.names=res$gns else row.names=NULL
	if(col.names) col.names=res$sns else col.names=NULL
	res2 <- crlmmGT(A=res[["A"]],
			B=res[["B"]],
			SNR=res[["SNR"]],
			mixtureParams=res[["mixtureParams"]],
			cdfName=cdfName,
			row.names=row.names,
			col.names=col.names,
			probs=probs,
			DF=DF,
			SNRMin=SNRMin,
			recallMin=recallMin,
			recallRegMin=recallRegMin,
			gender=gender,
			verbose=verbose,
			returnParams=returnParams,
			badSNP=badSNP)
	res2[["SNR"]] = res[["SNR"]]
	res2[["SKW"]] = res[["SKW"]]
	rm(res); gc(verbose=FALSE)
	return(list2SnpSet(res2, returnParams=returnParams))
}


## MR: Below is a more memory efficient version of crlmmIllumina() which
## reads in the .idats and genotypes in the one function and removes objects
## after they have been used
crlmmIlluminaV2 = function(sampleSheet=NULL,
			  arrayNames=NULL,
			  ids=NULL,
			  path=".",
			  arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
			  highDensity=FALSE,
			  sep="_",
			  fileExt=list(green="Grn.idat", red="Red.idat"),
			  saveDate=FALSE,
			  stripNorm=TRUE,
			  useTarget=TRUE,
 #                         outdir=".",
			  row.names=TRUE,
			  col.names=TRUE,
			  probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5, gender=NULL,
                          seed=1,  # save.it=FALSE, snpFile, cnFile,
                          mixtureSampleSize=10^5, eps=0.1, verbose=TRUE,
                          cdfName, sns, recallMin=10, recallRegMin=1000,
                          returnParams=FALSE, badSNP=.7) {

    if(missing(cdfName)) stop("must specify cdfName")
    if(!isValidCdfName(cdfName)) stop("cdfName not valid.  see validCdfNames")

#    is.lds = ifelse(isPackageLoaded("ff"), TRUE, FALSE)
    is.lds=FALSE
    RG = readIdatFiles(sampleSheet=sampleSheet, arrayNames=arrayNames,
                       ids=ids, path=path, arrayInfoColNames=arrayInfoColNames,
                       highDensity=highDensity, sep=sep, fileExt=fileExt, saveDate=saveDate)


    XY = RGtoXY(RG, chipType=cdfName)
#    if(is.lds) {
#      open(RG@assayData$R); open(RG@assayData$G); open(RG@assayData$zero)
#      delete(RG@assayData$R); delete(RG@assayData$G); delete(RG@assayData$zero)
#    }
    rm(RG); gc(verbose=FALSE)

    if (missing(sns)) { sns = sampleNames(XY)
    }
    res = preprocessInfinium2(XY, mixtureSampleSize=mixtureSampleSize, fitMixture=TRUE, verbose=verbose,
                               seed=seed, eps=eps, cdfName=cdfName, sns=sns, stripNorm=stripNorm, useTarget=useTarget) #, outdir=outdir) #,
#                               save.it=save.it, snpFile=snpFile, cnFile=cnFile)

#    if(is.lds) {
#      open(XY@assayData$X); open(XY@assayData$Y); open(XY@assayData$zero)
#      delete(XY@assayData$X); delete(XY@assayData$Y); delete(XY@assayData$zero)
#    }
    rm(XY); gc(verbose=FALSE)

    if(row.names) row.names=res$gns else row.names=NULL
    if(col.names) col.names=res$sns else col.names=NULL
##
##  if(is.lds){
##	  res2 <- crlmmGT2(A=res[["A"]],
##			   B=res[["B"]],
##			   SNR=res[["SNR"]],
##			   mixtureParams=res[["mixtureParams"]],
##			   cdfName=cdfName,
##			   row.names=row.names,
##			   col.names=col.names,
##			   probs=probs,
##			   DF=DF,
##			   SNRMin=SNRMin,
##			   recallMin=recallMin,
##			   recallRegMin=recallRegMin,
##			   gender=gender,
##			   verbose=verbose,
##			   returnParams=returnParams,
##			   badSNP=badSNP)
##  } else {
	  res2 <- crlmmGT(A=res[["A"]],
			   B=res[["B"]],
			   SNR=res[["SNR"]],
			   mixtureParams=res[["mixtureParams"]],
			   cdfName=cdfName,
			   row.names=row.names,
			   col.names=col.names,
			   probs=probs,
			   DF=DF,
			   SNRMin=SNRMin,
			   recallMin=recallMin,
			   recallRegMin=recallRegMin,
			   gender=gender,
			   verbose=verbose,
			   returnParams=returnParams,
			   badSNP=badSNP)
##  }

##    FUN = ifelse(is.lds, "crlmmGT2", "crlmmGT")
##    ## genotyping
##    crlmmGTfxn = function(FUN,...){
##		switch(FUN,
##		       crlmmGT2=crlmmGT2(...),
##		       crlmmGT=crlmmGT(...))
##              }
##    res2 = crlmmGTfxn(FUN,
##                     A=res[["A"]],
##                     B=res[["B"]],
##                     SNR=res[["SNR"]],
##                     mixtureParams=res[["mixtureParams"]],
##                     cdfName=cdfName,
##                     row.names=row.names,
##                     col.names=col.names,
##                     probs=probs,
##                     DF=DF,
##                     SNRMin=SNRMin,
##                     recallMin=recallMin,
##                     recallRegMin=recallRegMin,
##                     gender=gender,
##                     verbose=verbose,
##                     returnParams=returnParams,
##                     badSNP=badSNP)

#    if(is.lds) {
#      open(res[["SNR"]]); open(res[["SKW"]])
#    }
    res2[["SNR"]] = res[["SNR"]]
    res2[["SKW"]] = res[["SKW"]]
 #  if(is.lds) {
 #    delete(res[["A"]]); delete(res[["B"]])
 #    delete(res[["SKW"]]); delete(res[["SNR"]]); delete(res[["mixtureParams"]])
 #  }
    rm(res); gc(verbose=FALSE)
    return(list2SnpSet(res2, returnParams=returnParams))
}

# Functions analogous to Rob's Affy functions to set up container
getProtocolData.Illumina = function(filenames, sep="_", fileExt="Grn.idat", verbose=FALSE) {
       narrays = length(filenames)

       headerInfo = list(nProbes = rep(NA, narrays),
                         Barcode = rep(NA, narrays),
                         ChipType = rep(NA, narrays),
                         Manifest = rep(NA, narrays),
                         Position = rep(NA, narrays))

       scanDates = data.frame(ScanDate=rep(NA, narrays), DecodeDate=rep(NA, narrays))
       rownames(scanDates) <- make.unique(gsub(paste(sep, fileExt, sep=""), "", filenames))
       ## read in the data
       for(i in seq_along(filenames)) {
               if(verbose)
	               cat("reading", filenames[i], "\n")
	       idsG = G = NULL
	       G = readIDAT(filenames[i])
	       idsG = rownames(G$Quants)
	       headerInfo$nProbes[i] = G$nSNPsRead
	       headerInfo$Barcode[i] = G$Barcode
	       headerInfo$ChipType[i] = G$ChipType
	       headerInfo$Manifest[i] = G$Unknown$MostlyNull
	       headerInfo$Position[i] = G$Unknowns$MostlyA
               if(headerInfo$nProbes[i]>(headerInfo$nProbes[1]+10000) || headerInfo$nProbes[i]<(headerInfo$nProbes[1]-10000)) {
		       warning("Chips are not of the same type.  Skipping ", basename(filenames[i]))
		       next()
	       }
               if(nrow(G$RunInfo)>=2) {
   	           scanDates$ScanDate[i] = G$RunInfo[1, 1]
	           scanDates$DecodeDate[i] = G$RunInfo[2, 1]
               }
	       rm(G)
	       gc(verbose=FALSE)
       }
       protocoldata = new("AnnotatedDataFrame",
			    data=scanDates,
			    varMetadata=data.frame(labelDescription=colnames(scanDates),
			                           row.names=colnames(scanDates)))
       return(protocoldata)
}

getAvailableIlluminaGenomeBuild <- function(path){
	snp.file <- list.files(path, pattern="snpProbes_hg")
	if(length(snp.file) > 1){
		## use hg19
		message("genome build not specified. Using build hg19 for annotation.")
		snp.file <- snp.file[1]
	}
	if(length(snp.file) == 1)
		genome <- gsub(".rda", "", strsplit(snp.file, "snpProbes_")[[1]][[2]])
	if(length(snp.file)==0) genome <- ""
	genome
}


constructInf <- function(sampleSheet=NULL,
			 arrayNames=NULL,
			 path=".",
			 arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
			 highDensity=FALSE,
			 sep="_",
			 fileExt=list(green="Grn.idat", red="Red.idat"),
                         XY,
			 cdfName,
			 verbose=FALSE,
			 batch=NULL, #fns,
			 saveDate=TRUE) { #, outdir="."){
        if(is.null(XY)) {
  	  if(!is.null(arrayNames)) {
		pd = new("AnnotatedDataFrame", data = data.frame(Sample_ID=arrayNames))
	  }
	  if(!is.null(sampleSheet)) { # get array info from Illumina's sample sheet
		if(is.null(arrayNames)) {
			if(!is.null(arrayInfoColNames$barcode) && (arrayInfoColNames$barcode %in% colnames(sampleSheet))) {
				barcode = sampleSheet[,arrayInfoColNames$barcode]
				arrayNames=barcode
			}
			if(!is.null(arrayInfoColNames$position) && (arrayInfoColNames$position %in% colnames(sampleSheet))) {
				position = sampleSheet[,arrayInfoColNames$position]
				if(is.null(arrayNames))
					arrayNames=position
				else
					arrayNames = paste(arrayNames, position, sep=sep)
				if(highDensity) {
					hdExt = list(A="R01C01", B="R01C02", C="R02C01", D="R02C02")
					for(i in names(hdExt))
						arrayNames = sub(paste(sep, i, sep=""), paste(sep, hdExt[[i]], sep=""), arrayNames)
				}
			}
		  }
		  pd = new("AnnotatedDataFrame", data = sampleSheet)
		  sampleNames(pd) = make.unique(basename(arrayNames))
	  }
	  if(is.null(arrayNames)) {
		arrayNames = gsub(paste(sep, fileExt$green, sep=""), "", dir(pattern=fileExt$green, path=path))
		if(!is.null(sampleSheet)) {
			sampleSheet=NULL
			cat("Could not find required info in \'sampleSheet\' - ignoring.  Check \'sampleSheet\' and/or \'arrayInfoColNames\'\n")
		}
		pd = new("AnnotatedDataFrame", data = data.frame(Sample_ID=arrayNames))
  	  }
	  narrays = length(arrayNames)

	  if(!is.null(batch)) {
		stopifnot(length(batch) == narrays)
	  }
 	  if(is.null(batch)) {
                batch = rep("batch1", narrays) # assume only one batch stop("Must specify 'batch'")
	  }
          if(is(batch, "factor")) batch = as.character(batch)

	  grnfiles = paste(arrayNames, fileExt$green, sep=sep)
	  redfiles = paste(arrayNames, fileExt$red, sep=sep)
	  if(length(grnfiles)==0 || length(redfiles)==0)
		stop("Cannot find .idat files")
	  if(length(grnfiles)!=length(redfiles))
		stop("Cannot find matching .idat files")
	  if(path[1] != "." & path[1] != ""){
		grnidats = file.path(path, grnfiles)
		redidats = file.path(path, redfiles)
	  }  else {
		message("path arg not set.  Assuming files are in local directory, or that complete path is provided")
		grnidats = grnfiles
		redidats = redfiles
    	  }
	  if(!all(file.exists(grnidats))) stop("Missing some of the *Grn.idat files")
	  if(!all(file.exists(redidats))) stop("Missing some of the *Red.idat files")

	  if(verbose) message("Initializing container for genotyping and copy number estimation")
	  pkgname <- getCrlmmAnnotationName(cdfName)
	  path <- system.file("extdata", package=pkgname)
	  genome <- getAvailableIlluminaGenomeBuild(path)
	  featureData = getFeatureData(cdfName, copynumber=TRUE, genome=genome)
	  nr = nrow(featureData); nc = narrays
	  sns <-  if (!is.null(sampleSheet) && !is.null(sampleSheet$Sample_ID)) {
		make.unique(sampleSheet$Sample_ID)
          } else{
		make.unique(basename(arrayNames))
          }
	  biga <- initializeBigMatrix(name="A", nr, nc)
	  bigb <- initializeBigMatrix(name="B", nr, nc)
	  bigc <- initializeBigMatrix(name="call", nr, nc)
	  bigd <- initializeBigMatrix(name="callPr", nr,nc)
	  colnames(biga) <- colnames(bigb) <- colnames(bigc) <- colnames(bigd) <- sns
	  cnSet <- new("CNSet",
		     alleleA=biga,
		     alleleB=bigb,
		     call=bigc,
		     callProbability=bigd,
		     annotation=cdfName,
		     featureData=featureData,
		     batch=batch,
		     genome=genome)
##        if (!is.null(sampleSheet) && !is.null(sampleSheet$Sample_ID)) {
##		sampleNames(cnSet) = make.unique(sampleSheet$Sample_ID)
##        } else{
##		sampleNames(cnSet) <- make.unique(basename(arrayNames))
##	}
  	  if(saveDate){
		protocolData = getProtocolData.Illumina(grnidats, sep=sep, fileExt=fileExt$green, verbose=verbose)
	  } else{
		protocolData = annotatedDataFrameFrom(A(cnSet), byrow=FALSE)
	  }
	  rownames(pData(protocolData)) = sampleNames(cnSet)
	  protocolData(cnSet) = protocolData
	##featureData(cnSet) = featureData
	  featureNames(cnSet) = featureNames(featureData)
	  cnSet$gender <- initializeBigVector("gender", ncol(cnSet), vmode="integer")
	  cnSet$SNR = initializeBigVector("crlmmSNR-", ncol(cnSet), "double")
	  cnSet$SKW = initializeBigVector("crlmmSKW-", ncol(cnSet), "double")
	##sampleNames(cnSet) <- basename(sampleNames(cnSet))
        } else { # if XY specified, easier set-up of cnSet
          narrays = ncol(XY)
          if(verbose) message("Initializing container for genotyping and copy number estimation")           
          if(!is.null(batch)) {
              stopifnot(length(batch) == narrays)
          }
          if(is.null(batch)) {
              batch = rep("batch1", narrays) # assume only one batch stop("Must specify 'batch'")
          }
          if(is(batch, "factor")) batch = as.character(batch)
          pkgname <- getCrlmmAnnotationName(cdfName)
          path <- system.file("extdata", package=pkgname)
          genome <- getAvailableIlluminaGenomeBuild(path)
          featureData = getFeatureData(cdfName, copynumber=TRUE, genome=genome)
          nr = nrow(featureData); nc = narrays
          sns <-  sampleNames(XY)
          biga <- initializeBigMatrix(name="A", nr, nc)
          bigb <- initializeBigMatrix(name="B", nr, nc)
          bigc <- initializeBigMatrix(name="call", nr, nc)
          bigd <- initializeBigMatrix(name="callPr", nr,nc)
          colnames(biga) <- colnames(bigb) <- colnames(bigc) <- colnames(bigd) <- sns
          cnSet <- new("CNSet",
                       alleleA=biga,
                       alleleB=bigb,
                       call=bigc,
                       callProbability=bigd,
                       annotation=cdfName,
                       featureData=featureData,
                       batch=batch,
                       genome=genome)
          protocolData = annotatedDataFrameFrom(A(cnSet), byrow=FALSE)
          rownames(pData(protocolData)) = sampleNames(cnSet)
          protocolData(cnSet) = protocolData
          featureNames(cnSet) = featureNames(featureData)
          cnSet$gender = initializeBigVector("gender", ncol(cnSet), vmode="integer")
          cnSet$SNR = initializeBigVector("crlmmSNR-", ncol(cnSet), "double")
          cnSet$SKW = initializeBigVector("crlmmSKW-", ncol(cnSet), "double")
        }
	return(cnSet)
}
construct.Illumina <- constructInf

preprocessInf <- function(cnSet,
		       sampleSheet=NULL,
		       arrayNames=NULL,
		       ids=NULL,
		       path=".",
		       arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
		       highDensity=TRUE,
		       sep="_",
		       fileExt=list(green="Grn.idat", red="Red.idat"),
                       XY,
		       saveDate=TRUE,
		       stripNorm=TRUE,
		       useTarget=TRUE,
		       mixtureSampleSize=10^5,
		       fitMixture=TRUE,
                       quantile.method="between",
		       eps=0.1,
		       verbose=TRUE,
		       seed=1, cdfName){
	stopifnot(all(c("gender", "SNR", "SKW") %in% varLabels(cnSet)))
	open(A(cnSet))
	open(B(cnSet))
	sns <- sampleNames(cnSet)
	pkgname = getCrlmmAnnotationName(annotation(cnSet))
	is.snp = isSnp(cnSet)
	snp.index = which(is.snp)
	narrays = ncol(cnSet)
	sampleBatches <- splitIndicesByLength(seq_len(ncol(cnSet)), ocSamples()/getDoParWorkers())
	mixtureParams = initializeBigMatrix("crlmmMixt-", 4, narrays, "double")
	ocLapply(seq_along(sampleBatches),
		 processIDAT, # crlmm:::
		 sampleBatches=sampleBatches,
		 sampleSheet=sampleSheet,
		 arrayNames=arrayNames,
		 ids=ids,
		 path=path,
		 arrayInfoColNames=arrayInfoColNames,
		 highDensity=highDensity,
		 sep=sep,
		 fileExt=fileExt,
                 XY=XY,
		 saveDate=saveDate,
		 verbose=verbose,
		 mixtureSampleSize=mixtureSampleSize,
		 fitMixture=fitMixture,
		 eps=eps,
		 seed=seed,
		 cdfName=cdfName,
		 sns=sns,
		 stripNorm=stripNorm,
		 useTarget=useTarget,
                 quantile.method=quantile.method,                 
		 A=A(cnSet),
		 B=B(cnSet),
		 GT=calls(cnSet),
		 GTP=snpCallProbability(cnSet),
		 SKW=cnSet$SKW,
		 SNR=cnSet$SNR,
		 mixtureParams=mixtureParams,
		 is.snp=is.snp,
		 neededPkgs=c("crlmm", pkgname)) # outdir=outdir,
	return(mixtureParams)
}
preprocess <- preprocessInf

genotypeInf <- function(cnSet, mixtureParams, probs=rep(1/3,3),
			SNRMin=5,
			recallMin=10,
			recallRegMin=1000,
			verbose=TRUE,
			returnParams=TRUE,
			badSNP=0.7,
			gender=NULL,
			DF=6,
			cdfName,
                        call.method="crlmm",
                        trueCalls=NULL){
	is.snp = isSnp(cnSet)
	snp.index = which(is.snp)
##	narrays = ncol(cnSet)
##	open(A(cnSet))
##	open(B(cnSet))
##	open(snpCall(cnSet))
##	open(snpCallProbability(cnSet))
##	## crlmmGT2 overwrites the normalized intensities with calls and confidenceScores
##	message("Writing to call and callProbability slots")
##	for(j in 1:ncol(cnSet)){
##		snpCall(cnSet)[snp.index, j] <- as.integer(A(cnSet)[snp.index, j])
##		snpCallProbability(cnSet)[snp.index, j] <- as.integer(B(cnSet)[snp.index, j])
##	}
##	message("Writing complete.  Begin genotyping...")
##	close(A(cnSet))
##	close(B(cnSet))
        if(call.method=="crlmm") {
  	  tmp <- crlmmGT2(A=A(cnSet),
			B=B(cnSet),
			SNR=cnSet$SNR,
			mixtureParams=mixtureParams,
			cdfName=cdfName,
			col.names=sampleNames(cnSet),
			probs=probs,
			DF=DF,
			SNRMin=SNRMin,
			recallMin=recallMin,
			recallRegMin=recallRegMin,
			gender=gender,
			verbose=verbose,
			returnParams=returnParams,
			badSNP=badSNP,
			callsGt=calls(cnSet),
			callsPr=snpCallProbability(cnSet))#,
	##RS:  I changed the API for crlmmGT2 to be consistent with crlmmGT
           open(cnSet$gender)
           cnSet$gender[,] = tmp$gender
           close(cnSet$gender)
        }
        if(call.method=="krlmm")
          tmp <- krlmm(cnSet, cdfName, gender=gender, trueCalls=trueCalls, verbose=verbose) # new function required...  cnSet, cdfName and gender are probably the only arguments you need...
	## snp.names=featureNames(cnSet)[snp.index])
	if(verbose) message("Genotyping finished.") # Updating container with genotype calls and confidence scores.")
	TRUE
}


genotype.Illumina <- function(sampleSheet=NULL,
			      arrayNames=NULL,
			      ids=NULL,
			      path=".",
			      arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
			      highDensity=FALSE,
			      sep="_",
			      fileExt=list(green="Grn.idat", red="Red.idat"),
                              XY=NULL,
                              call.method="crlmm",
                              trueCalls=NULL,
			      cdfName,
			      copynumber=TRUE,
			      batch=NULL,
			      saveDate=FALSE,
			      stripNorm=TRUE,
			      useTarget=TRUE,
                              quantile.method="between",                             
			      mixtureSampleSize=10^5,
			      fitMixture=TRUE,
			      eps=0.1,
			      verbose=TRUE,
			      seed=1,
			      sns,
			      probs=rep(1/3, 3),
			      DF=6,
			      SNRMin=5,
			      recallMin=10,
			      recallRegMin=1000,
			      gender=NULL,
			      returnParams=TRUE,
			      badSNP=0.7) {
        krlmm.supported = c("humanomni1quadv1b",      # Omni1 quad
	                    "humanomni25quadv1b",     # Omni2.5 quad
	                    "humanomni258v1a",        # Omni2.5 8 v1 A
                            "humanomni258v1p1b",      # Omni2.5 8 v1.1 B
                            "humanomni5quadv1b",      # Omni5 quad
			    "humanexome12v1p2a",      # Exome 12 v1.2 A
                            "humanomniexpexome8v1p1b") # Omni Express Exome 8 v1.1b
        crlmm.supported = c("human1mv1c",             # 1M
                            "human370v1c",            # 370CNV
	                    "human650v3a",            # 650Y
	                    "human610quadv1b",        # 610 quad
	                    "human660quadv1a",        # 660 quad
	                    "human370quadv3c",        # 370CNV quad
	                    "human550v3b",            # 550K
	                    "human1mduov3b",          # 1M Duo
                            "humanomni1quadv1b",      # Omni1 quad
#			    "humanomni258v1a",        # Omni2.5 8 v1 A
#                           "humanomni258v1p1b",      # Omni2.5 8 v1.1 B
	                    "humanomniexpress12v1b",  # Omni express 12
	                    "humanimmuno12v1b",       # Immuno chip 12
                            "humancytosnp12v2p1h",    # CytoSNP 12
                            "humanomniexpexome8v1p1b") # Omni Express Exome 8 v1.1b
        if(call.method=="krlmm") {
          if(!any(cdfName==krlmm.supported))
             stop(cdfName, " platform not supported by krlmm.  Consider setting call.method=\'crlmm\'")           
          fitMixture=FALSE
          quantile.method="within"
        }
        if(call.method=="crlmm") {
          if(!any(cdfName==crlmm.supported))
             stop(cdfName, " platform not supported by crlmm.  Consider setting call.method=\'krlmm\'")
          fitMixture=TRUE
          # quantile.method="between"
        }
	is.lds = ifelse(isPackageLoaded("ff"), TRUE, FALSE)
        if (!(is.lds))
            stop("Package ff not loaded")
        if(!is.null(XY) && missing(cdfName))
          cdfName = getCrlmmAnnotationName(annotation(cnSet))
	if(missing(cdfName)) stop("must specify cdfName")
	if(!isValidCdfName(cdfName)) stop("cdfName not valid.  see validCdfNames")
#	stopifnot(!missing(batch))
#	if(is(batch, "factor")) batch <- as.character(batch)
        pkgname <- getCrlmmAnnotationName(cdfName)
#        if(is.null(cnSet)) {
  	  message("Instantiate CNSet container.")
	  cnSet <- constructInf(sampleSheet=sampleSheet,
				    arrayNames=arrayNames,
				    path=path,
				    arrayInfoColNames=arrayInfoColNames,
				    highDensity=highDensity,
				    sep=sep,
				    fileExt=fileExt,
                                    XY=XY,
				    cdfName=cdfName,
				    verbose=verbose,
				    batch=batch,
				    saveDate=saveDate)
#        }
        if(call.method=="krlmm" && ncol(cnSet)<8)
           stop(paste("To run krlmm you need at least 8 samples (in general, the more the better).  You currently have", ncol(cnSet)))
	mixtureParams <- preprocessInf(cnSet=cnSet,
				    sampleSheet=sampleSheet,
				    arrayNames=arrayNames,
				    ids=ids,
				    path=path,
				    arrayInfoColNames=arrayInfoColNames,
				    highDensity=highDensity,
				    sep=sep,
				    fileExt=fileExt,
				    saveDate=saveDate,
                                    XY=XY,
				    stripNorm=stripNorm,
				    useTarget=useTarget,
				    mixtureSampleSize=mixtureSampleSize,
				    fitMixture=fitMixture,
                                    quantile.method=quantile.method,
				    eps=eps,
				    verbose=verbose,
				    seed=seed,
				    cdfName=cdfName)
	message("Preprocessing complete.  Begin genotyping...")
	genotypeInf(cnSet=cnSet, mixtureParams=mixtureParams,
		    probs=probs,
		    SNRMin=SNRMin,
		    recallMin=recallMin,
		    recallRegMin=recallRegMin,
		    verbose=verbose,
		    returnParams=returnParams,
		    badSNP=badSNP,
		    gender=gender,
		    DF=DF,
		    cdfName=cdfName,
                    call.method=call.method,
                    trueCalls=trueCalls)
	return(cnSet)
}


processIDAT <- function(stratum, sampleBatches, sampleSheet=NULL,
			arrayNames=NULL,
			ids=NULL,
			path=".",
			arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
			highDensity=FALSE,
			sep="_",
			fileExt=list(green="Grn.idat", red="Red.idat"),
                        XY,
			saveDate=FALSE,
			verbose=TRUE,
			mixtureSampleSize=10^5,
			fitMixture=TRUE,
			eps=0.1,
			seed=1,
			cdfName,
			sns,
			stripNorm=TRUE,
			useTarget=TRUE,
                        quantile.method=quantile.method,                        
			A, B,
			GT,
			GTP,
			SKW, SNR, mixtureParams, is.snp) { #, outdir=".") {
	message("Processing sample stratum ", stratum, " of ", length(sampleBatches))
	sel <- sampleBatches[[stratum]]
        if(length(path)>= length(sel)) path = path[sel]
        if(is.null(XY)) {
          RG = readIdatFiles(sampleSheet=sampleSheet[sel,], arrayNames=arrayNames[sel],
                       ids=ids, path=path, arrayInfoColNames=arrayInfoColNames,
                       highDensity=highDensity, sep=sep, fileExt=fileExt, saveDate=saveDate, verbose=verbose)
          XY = RGtoXY(RG, chipType=cdfName)
          rm(RG)
          gc(verbose=FALSE)
          if (missing(sns) || length(sns)!=ncol(XY)) sns = sampleNames(XY)
          res = preprocessInfinium2(XY, mixtureSampleSize=mixtureSampleSize, fitMixture=TRUE, verbose=verbose,
                               seed=seed, eps=eps, cdfName=cdfName, sns=sns, stripNorm=stripNorm, useTarget=useTarget,
                               quantile.method=quantile.method) #, outdir=outdir)
          rm(XY)
        }else{ #XY already available
          if (missing(sns) || length(sns)!=ncol(XY)) sns = sampleNames(XY)  
          res = preprocessInfinium2(XY[,sel], mixtureSampleSize=mixtureSampleSize, fitMixture=TRUE, verbose=verbose,
                                           seed=seed, eps=eps, cdfName=cdfName, sns=sns[sel], stripNorm=stripNorm, useTarget=useTarget,
                                           quantile.method=quantile.method) 
        }
        gc(verbose=FALSE)
	if(verbose) message("Finished preprocessing.")
        snp.index = which(is.snp)
	np.index = which(!is.snp)
	open(A); open(B);
	open(GT); open(GTP)
	Amatrix <- res[["A"]]
	Bmatrix <- res[["B"]]

	## Amatrix and Bmatrix are ordinary matrices--not ff objects.
	## By writing columns of a ordinary matrix to GT and GTP, we
	## save one read step later on.  GT and GTP will be updated to
	## calls and call probabilities by the crlmmGT2 function. The A
	## and B slots will retain the normalized intensities for the
	## A and B alleles
	for(j in seq_along(sel)){
		jj <- sel[j]
		A[snp.index, jj] <- Amatrix[, j]
		GT[snp.index, jj] <- Amatrix[, j]
		B[snp.index, jj] <- Bmatrix[, j]
		GTP[snp.index, jj] <- Bmatrix[, j]
	}
	if(length(np.index)>0) {
		cnAmatrix <- res[["cnAB"]]$A
		cnBmatrix <- res[["cnAB"]]$B
		for(j in seq_along(sel)){
			jj <- sel[j]
			A[np.index, jj] <- cnAmatrix[, j]
			GT[np.index, jj] <- cnAmatrix[, j]
			B[np.index, jj] <- cnBmatrix[, j]
			GTP[np.index, jj] <- cnBmatrix[, j]
		}
        }
	open(SKW); open(SNR); open(mixtureParams)
	SKW[sel] = res[["SKW"]][]
	SNR[sel] = res[["SNR"]][]
	mixtureParams[, sel] = res[["mixtureParams"]][]
        close(A); close(B)
	close(GT); close(GTP)
        close(SNR); close(SKW)
        close(mixtureParams)
        rm(res)
	gc(verbose=FALSE)
        TRUE
      }
