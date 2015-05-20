krlmm <- function(cnSet, cdfName, gender=NULL, trueCalls=NULL, verbose=TRUE) {
  
    pkgname <- getCrlmmAnnotationName(cdfName) # crlmm:::
    
### pre-processing, output M    
    M = computeLogRatio(cnSet, verbose)
    message("leaving out novariant SNPs")
    
    # For SNPs with less than 3 distinct data point, exclude them from downstream analysis
    uniqueCount = apply(M[,], 1, function(x){length(unique(x))})
    SNPtoProcessIndex = uniqueCount >= 3
    noVariantSNPIndex = uniqueCount < 3
    M = M[SNPtoProcessIndex, ]

    numSNP = nrow(M)
    numSample = ncol(M)
    
    calls = oligoClasses::initializeBigMatrix(name="calls", numSNP, numSample, vmode = "integer")
    scores = oligoClasses::initializeBigMatrix(name="scores", numSNP, numSample, vmode = "double")
    open(calls)
    open(scores)

    rownames(calls) = rownames(M)
    rownames(scores) = rownames(M)
    colnames(calls) = colnames(M)
    colnames(scores) = colnames(M)  
    
    priormeans = calculatePriorValues(M, numSNP, verbose)
    VGLMparameters = calculateParameters(M, priormeans, numSNP, verbose)
  
### retrieve or calculate coefficients
    krlmmCoefficients = getKrlmmVGLMCoefficients(pkgname, trueCalls, VGLMparameters, verbose, numSample, colnames(M))
    
### do VGLM fit, to predict k for each SNP
    kPrediction <- predictKwithVGLM(VGLMparameters, krlmmCoefficients, verbose);
    rm(VGLMparameters)
    
### assign calls
    assignCalls(calls, M, kPrediction, priormeans, numSNP, numSample, verbose);
    rm(kPrediction)

### assign confidence scores
    computeCallPr(scores, M, calls, numSNP, numSample, verbose)
    
### add back SNPs excluded before
    AddBackNoVarianceSNPs(cnSet, calls, scores, numSNP, numSample, SNPtoProcessIndex, noVariantSNPIndex)

    loader("genotypeStuff.rda", .crlmmPkgEnv, pkgname)
    YIndex <- getVarInEnv("YIndex")
    loader("annotation.rda", .crlmmPkgEnv, pkgname)
    annotation <- getVarInEnv("annot")
    
### impute gender if gender information not provided
    if (is.null(gender)) {
        gender = krlmmImputeGender(cnSet, annotation, YIndex, verbose)
    }

### double-check ChrY SNP cluster, impute gender if gender information not provided
    if (!(is.null(gender))) {
        verifyChrYSNPs(cnSet, M, calls, gender, annotation, YIndex, priormeans, verbose)
    }

    close(calls)
    close(scores)
    rm(calls)
    rm(scores)
    rm(M)
    rm(annotation)
    rm(YIndex)
    TRUE
}

#######################################################################################################

getNumberOfCores <- function(){
    defaultCores = min(detectCores(), 8)
    return(getOption("krlmm.cores", defaultCores))    
}


#######################################################################################################

computeLogRatio <- function(cnSet, verbose, blockSize = 500000){
    # compute log-ratio in blocksize of 500,000 by default
    message("Start computing log ratio")
    A <- A(cnSet)
    B <- B(cnSet)
    open(A)
    open(B)
    
    numSNP <- nrow(A)
    numSample <- ncol(A)
    
    M <- oligoClasses::initializeBigMatrix(name="M", numSNP, numSample, vmode = "double")
    rownames(M) = rownames(A)
    colnames(M) = colnames(A)
  
    numBlock = ceiling(numSNP / blockSize)
    for (i in 1:numBlock){
        if (verbose) message(" -- Processing segment ", i, " out of ", numBlock)
        subsetA = as.matrix(A[((i-1) * (blockSize) + 1):min(i * blockSize, numSNP), ])
        subsetB = as.matrix(B[((i-1) * (blockSize) + 1):min(i * blockSize, numSNP), ])                     
        subsetM = .Call("krlmmComputeM", subsetA, subsetB, PACKAGE="crlmm")
        M[((i-1) * (blockSize) + 1):min(i * blockSize, numSNP),] = subsetM[1:nrow(subsetM), ]
        rm(subsetA, subsetB, subsetM)
    }
    
    close(A)
    close(B)
    if (verbose) message("Done computing log ratio")
    return(M)    
}

computeAverageLogIntensity <- function(cnSet, verbose, blockSize = 500000){
    # compute average log intensity in blocksize of 500,000 by default
    message("Start computing average log intensity")
    A <- A(cnSet)
    B <- B(cnSet)
    open(A)
    open(B)
    
    numSNP <- nrow(A)
    numSample <- ncol(A)
    
    S <- oligoClasses::initializeBigMatrix(name="S", numSNP, numSample, vmode = "double")
    rownames(S) = rownames(A)
    colnames(S) = colnames(A)
      
    numBlock = ceiling(numSNP / blockSize)
    for (i in 1:numBlock){
        if (verbose) message(" -- Processing segment ", i, " out of ", numBlock)
        subsetA = as.matrix(A[((i-1) * (blockSize) + 1):min(i * blockSize, numSNP), ])
        subsetB = as.matrix(B[((i-1) * (blockSize) + 1):min(i * blockSize, numSNP), ])                     
        subsetS = .Call("krlmmComputeS", subsetA, subsetB, PACKAGE="crlmm")
        S[((i-1) * (blockSize) + 1):min(i * blockSize, numSNP),] = subsetS[1:nrow(subsetS), ]
        rm(subsetA, subsetB, subsetS)
    }
    
    close(A)
    close(B)
    if (verbose) message("Done computing average log intensity")
    return(S)
  
}


#######################################################################################################

calculatePriorValues <- function(M, numSNP, verbose) {

    calculateOneKmeans <- function(x) {
        tmp = kmeans(x, 3, nstart=45)
        return(sort(tmp$centers, decreasing = T))
    }

    if (verbose) message("Start calculating Prior Means")
    cl <- makeCluster(getNumberOfCores())
    centers <- parRapply(cl, M, calculateOneKmeans)
    stopCluster(cl) 
    centers <- matrix(centers, numSNP, 3, byrow = TRUE)
    priormeans = apply(centers, 2, FUN="median", na.rm=TRUE)
    if(abs(sum(priormeans))>1) {
      checksymmetric= apply(centers,1,function(x){abs(sum(x))})<1
      priormeans=apply(centers[checksymmetric,],2, FUN="median", na.rm=TRUE)
    }
    if (verbose) message("Done calculating Prior Means")
    return(priormeans)
}

#######################################################################################################

calculateKrlmmCoefficients <- function(trueCalls, params, numSample, samplenames){
#    if (!require(VGAM)) {
#        message("VGAM package not found, fall back to use defined platform-specific coefficients")
#        return(NULL)
#    }
    amatch = match(rownames(params), rownames(trueCalls))
    amatch = amatch[!is.na(amatch)]
    trueCalls = trueCalls[amatch,]

    amatch = match(rownames(trueCalls), rownames(params))
    params = params[amatch, ]
    
    amatch = match(samplenames, colnames(trueCalls))
    amatch = amatch[!is.na(amatch)]    
    trueCalls = trueCalls[, amatch]
        
    if ((numSample <= 40) && (ncol(trueCalls) < round(numSample/2))){
        message("Sample size is ", numSample, ", KRLMM requires at least trueCalls of ", round(numSample/2), " samples to calculate coefficients")
        return(NULL)
    }
    if ((numSample > 40) && (ncol(trueCalls) < 20)){
        message("At least trueCalls of 20 samples are required to calculate coefficients")
        return(NULL)
    }
    if (nrow(trueCalls) < 1000){
        message("At lease trueCalls of 1000 SNPs are required to calculate coefficients")
        return(NULL)
    }

    getna = apply(trueCalls, 1, function(x){sum(is.na(x))>=10})
    truek = apply(trueCalls, 1, function(x){length(unique(x[!is.na(x)]))})
   
    params1 = params[!getna,]
    truek1 = truek[!getna]

    xx = data.frame(params1)
    t = as.factor(as.numeric(truek1)) 
    xx = data.frame(xx, t)
    fit = suppressWarnings(vglm(t~., multinomial(refLevel=1), xx))
    coe = coefficients(fit) # VGAM::
    return(coe)    
}

getKrlmmVGLMCoefficients <- function(pkgname, trueCalls, params, verbose, numSample, samplenames){
    if (!is.null(trueCalls)) {
        coe = calculateKrlmmCoefficients(trueCalls, params, numSample, samplenames)
        if (!is.null(coe)) {
            if (verbose) message ("Done calculating platform-specific coefficients to predict number of clusters")
            return(coe)
        }
    }
    if (!is.null(trueCalls)) message("Fall back to use defined platform-specific coefficients")
    if (verbose) message ("Retrieving defined platform-specific coefficients to predict number of clusters")
    loader("krlmmVGLMCoefficients.rda", .crlmmPkgEnv, pkgname)
    return(getVarInEnv("krlmmCoefficients"))      
}


predictKwithVGLM <- function(data, coe, verbose){
    if (verbose) message("Start predicting number of clusters")
    logit1 <- rep(0, nrow(data))
    logit2 <- coe[1]+coe[3]*data[,1]+coe[5]*data[,2]+coe[7]*data[,3]+coe[9]*data[,4]+coe[11]*data[,5]+coe[13]*data[,6]+coe[15]*data[,7]+coe[17]*data[,8]
    logit23 <- coe[2]+coe[4]*data[,1]+coe[6]*data[,2]+coe[8]*data[,3]+coe[10]*data[,4]+coe[12]*data[,5]+coe[14]*data[,6]+coe[16]*data[,7]+coe[18]*data[,8]
   
    logits <- cbind(logit1, logit2, logit23)
    rm(logit1)
    rm(logit2)
    rm(logit23)
    p.unscaled <- exp(logits)
    rm(logits)   
    p <- p.unscaled / rowSums(p.unscaled)
    clusterPrediction = apply(p, 1, function(x){which.max(x)})
    rm(p.unscaled)
    rm(p)
    if (verbose) message("Done predicting number of clusters")
    return(clusterPrediction)
}

#######################################################################################################

assignCalls3Cluster <- function(intensities, priormeans){
    prior31 = c(priormeans[1]/2,priormeans[2],priormeans[3])
    prior32 = c(priormeans[1],priormeans[2],priormeans[3]/2)	

    emp <- rep(NA, length(intensities))
	ansCalls <- emp
    tmp <- try(kmeans(intensities, priormeans, nstart=1),TRUE)
    if(class(tmp) == "try-error") {
        tmp1 <- try(kmeans(intensities, prior31, nstart=1),TRUE)
        if(class(tmp1) == "try-error"){
            tmp2 <- try(kmeans(intensities, prior32, nstart=1),TRUE)
            if(class(tmp2) == "try-error"){
                ansCalls = emp
            }else{
                ansCalls = tmp2$cluster
            }
        }else{
            ansCalls = tmp1$cluster
        }
    }else{
        ansCalls = tmp$cluster
    }
    rm(prior31, prior32)
    return(ansCalls)
 
}

#######################################################################################################

assignCalls2Cluster <- function(intensities, priormeans){
    closest <- rep(NA, 3)	
    for(j in 1:3){
        distance <- as.vector(abs(priormeans[j]-intensities))
        closest[j] <- intensities[which.min(distance)]
    }
    prior2 <- priormeans[priormeans!=priormeans[which.max(abs(closest-priormeans))]]
	
    emp <- rep(NA, length(intensities))
    ansCalls <- emp    
    tmp <- try(kmeans(intensities, prior2, nstart=1), TRUE)
    if(class(tmp) == "try-error") {
        ansCalls <- emp
    }else{
        if(prior2[1]==priormeans[2] && prior2[2]==priormeans[3]){
            mp <- abs(tmp$centers[1]-tmp$centers[2])
            pp <- abs(priormeans[2]-priormeans[3])
            if((mp/pp)<=0.25){
                ansCalls <- emp
            }else{
                ansCalls <- tmp$cluster+1
            }
        }
        if(prior2[1]==priormeans[1] && prior2[2]==priormeans[2]){
            mp=abs(tmp$centers[1]-tmp$centers[2])
            pp=abs(priormeans[1]-priormeans[2])
            if((mp/pp)<=0.25){
                ansCalls = emp
            }else{
                ansCalls <- tmp$cluster
            }
        }
        if(prior2[1]==priormeans[1] && prior2[2]==priormeans[3]){
            mp <- abs(tmp$centers[1]-tmp$centers[2])
            pp <- abs(priormeans[1]-priormeans[3])
            if ((mp/pp) <=0.25){
                ansCalls <- emp
            }else{
                ansCalls[tmp$cluster==1] <- 1
                ansCalls[tmp$cluster==2] <- 3
            }
        }
    }
    rm(tmp)
    return(ansCalls)
}

#######################################################################################################

assignCalls1Cluster <- function(intensities, priormeans){
    closest <- rep(NA, 3)
    for(j in 1:3){
        distance <- as.vector(abs(priormeans[j]-intensities))
        closest[j]=intensities[which.min(distance)]
    }
    clusterindex <- which.min(abs(closest-priormeans))
    return(rep(clusterindex, length(intensities)))
}
  
#######################################################################################################

assignCallsOneSNP <- function(x, priormeans, numSample){
    tolerance = 1e-3
       
    k = x[numSample + 1]
    values = x[1:numSample]  
    
    if (abs(k - 2) <= tolerance) {
        tmp <- try(kmeans(values, priormeans, nstart=1),TRUE)
        if (!(class(tmp)=="try-error")) {
            k <- 3;
        }
    }

    successful <- FALSE;
    if (abs(k - 3) <= tolerance){
        SNPCalls <- assignCalls3Cluster(values, priormeans)
        successful <- !is.na(SNPCalls[1])
        if (!successful) { 
            k <- 2
        }
    }
		
    if ( (abs(k - 2) <= tolerance) && (!successful)){
        SNPCalls <- assignCalls2Cluster(values, priormeans)
        successful <- !is.na(SNPCalls[1])
        if (!successful) { 
            k <- 1
        }			
    }

    if ( (abs(k - 1) <= tolerance) && (!successful)){
        SNPCalls <- assignCalls1Cluster(values, priormeans)
    }
    return(SNPCalls)
}


assignCalls <- function(callsGt, M, a, priormeans, numSNP, numSample, verbose, blockSize=500000){
    # process by block size of 500,000 by default
    message("Start assign calls")
       
    numBlock = ceiling(numSNP / blockSize)

    for (i in 1:numBlock){
        if (verbose) message(" -- Processing segment ", i, " out of ", numBlock)
        subsetM = as.matrix(M[((i-1) * (blockSize) + 1):min(i * blockSize, numSNP), ])
        thisnumSNP = nrow(subsetM)
        subseta = a[((i-1) * (blockSize) + 1):min(i * blockSize, numSNP)]
        subseta = matrix(subseta, thisnumSNP, 1)

        testValues = cbind(subsetM, subseta)

        cl <- makeCluster(getNumberOfCores())
        callAns <- parRapply(cl, testValues, assignCallsOneSNP, priormeans, numSample)
        stopCluster(cl)

        callAnsAllValues = unlist(callAns)

        subsetcalls <- matrix(callAnsAllValues, thisnumSNP, numSample, byrow = TRUE)
        
        callsGt[((i-1) * (blockSize) + 1):min(i * blockSize, numSNP),] = subsetcalls[1:thisnumSNP, ]
        rm(subsetM, subseta, subsetcalls)
    }
    message("Done assign calls")
}

#######################################################################################################


calculateParameters <- function(M, priormeans, numSNP, verbose) {
    if (verbose) message("Start calculating 3-clusters parameters")
    params3cluster <- calculateParameters3Cluster(M, priormeans, numSNP, verbose);
    if (verbose) message("Done calculating 3-cluster parameters")

    if (verbose) message("Start calculating 2-cluster parameters")
    params2cluster <- calculateParameters2Cluster(M, priormeans, numSNP, verbose);
    if (verbose) message("Done calculating 2-cluster parameters")

    if (verbose) message("Start calculating 1-cluster parameters")
    params1cluster <- calculateParameters1Cluster(M, priormeans, numSNP, verbose);    
    if (verbose) message("Done calculating 1-cluster parameters")

    parameters <- cbind(as.matrix(params1cluster$elemh), as.matrix(params1cluster$eless), as.matrix(params2cluster$elemh), as.matrix(params2cluster$eless),
                        as.matrix(params3cluster$elemh), as.matrix(params3cluster$eless), as.matrix(params2cluster$elehw), as.matrix(params3cluster$elehw));
    
    rownames(parameters) = rownames(M)
    rm(params3cluster)
    rm(params2cluster)
    rm(params1cluster)
    return(parameters)
}

hardyweinberg <- function(x, minN = 8){
   if (length(x) < minN){
       return(NA) 
    } else {
       result = .Call("krlmmHardyweinberg", x)
       return(result)
    }
}

calculateOneParams3Cluster <- function(x, priors, priorDown, priorUp){

    tmp = try(kmeans(x, priors, nstart=1), TRUE)
    if(class(tmp) == "kmeans") {
        flag = 1
     } else {
        tmp = try(kmeans(x, priorDown, nstart=1), TRUE)
        if(class(tmp) == "kmeans") {
            flag = 2
        } else {
            tmp = try(kmeans(x, priorUp, nstart=1), TRUE)
            if (class(tmp) == "kmeans") {
                flag = 3
            } else {
                tmp = kmeans(x, 3, nstart=35)
                flag = 4
            }
        }  
    }
    ss = sum(tmp$withinss)
    hw = hardyweinberg(tmp$cluster) 
    centers = sort(tmp$centers, decreasing = TRUE)

    return(c(centers, ss, hw, flag))
}

calculateParameters3Cluster <- function(M, priormeans, numSNP, verbose) {
    Apriors = priormeans
    ApriorDown = c(Apriors[1]/2, Apriors[2], Apriors[3]) # shift-down
    ApriorUp = c(Apriors[1], Apriors[2], Apriors[3]/2) # shift-up

    cl <- makeCluster(getNumberOfCores())
    clusterEvalQ(cl, library(crlmm))
    parAns <- parRapply(cl, M, calculateOneParams3Cluster, Apriors, ApriorDown, ApriorUp)

    stopCluster(cl)
    parAnsAllValues = unlist(parAns)
    parameters <- matrix(parAnsAllValues, numSNP, 6, byrow = TRUE)
   
    centers <- parameters[, 1:3]
    ss <- parameters[, 4]
    hw <- parameters[, 5]
    flag <- parameters[, 6]

    rm(parAns)
    rm(parameters)
   
    sigma=solve(cov(centers, use="complete.obs"))
    mh = calculateMahalDist3Cluster(centers, sigma, flag, Apriors, ApriorDown, ApriorUp, numSNP)

    rm(sigma)
    rm(centers)

    gc()
    return(list(eless = ss, elemh = mh, elehw = hw))
}


calculateMahalDist3Cluster <- function(centers, sigma, flag, priors, priorDown, priorUp, numSNP){
    mahaldist = rep(NA, numSNP)

    tolerance = 1e-3
    for (i in 1:numSNP) {
        if ((abs(flag[i] - 1) <= tolerance) || (abs(flag[i]- 4) <= tolerance)) difference = centers[i, ] - priors
        if (abs(flag[i] - 2) <= tolerance) difference = centers[i, ] - priorDown
        if (abs(flag[i] - 3) <= tolerance) difference = centers[i, ] - priorUp         
        mahaldist[i] = as.vector(difference)%*%sigma%*%as.vector(difference)
    }
    return(mahaldist)
}


calculateOneParams2Cluster <- function(x, priors){
    aa = rep(NA, 3)   
    for(j in 1:3){
        dist = as.vector(abs(priors[j]-x))
        aa[j]=x[which.min(dist)]
    } 
    prior2 = priors[priors!=priors[which.max(abs(aa-priors))]]
     

    tmp = try(kmeans(x, prior2, nstart=1), TRUE)
    if (class(tmp)=="kmeans") {
        centers = tmp$centers
        hw = hardyweinberg(tmp$cluster)
    }
    rm(tmp)
    tmp = kmeans(x, 2, nstart = 35)
    if (tmp$centers[1] < tmp$centers[2]){
        centers = c(tmp$centers[2], tmp$centers[1])
        hw = hardyweinberg(3 - tmp$cluster)   
    } else {
        centers = tmp$centers
        hw = hardyweinberg(tmp$cluster)
    }
    ss = sum(tmp$withinss)     
    return(c(centers, prior2, ss, hw))
}


calculateParameters2Cluster <- function(M, priormeans, numSNP, verbose) {
    Apriors = priormeans
   
    cl <- makeCluster(getNumberOfCores())
    clusterEvalQ(cl, library(crlmm))
    parAns <- parRapply(cl, M, calculateOneParams2Cluster, Apriors)
    stopCluster(cl)

    parAnsAllValues = unlist(parAns)
    parameters <- matrix(parAnsAllValues, numSNP, 6, byrow = TRUE)

    centers <- parameters[, 1:2]
    priors2cluster <- parameters[, 3:4]
    ss <- parameters[, 5]
    hw <- parameters[, 6]
    
    rm(parAns)
    rm(parameters)
    sigma=solve(cov(centers, use="complete.obs"))

    mh = calculateMahalDist2Cluster(centers, sigma, priors2cluster, numSNP)

    rm(sigma)
    rm(centers)

    gc()
    return(list(eless = ss, elemh = mh, elehw = hw))
}


calculateMahalDist2Cluster <- function(centers, sigma, priors, numSNP){
    mahaldist = rep(NA, numSNP)
    
    for (i in 1:numSNP) {
        difference <- centers[i,] - priors[i,]
        mahaldist[i] = as.vector(difference)%*%sigma%*%as.vector(difference)        
    }
    return(mahaldist)
}

calculateOneParams1Cluster <- function(x, priors){
    center = mean(x)
    diff <- x - center
    diffsquare <- diff^2
    ss = sum(diffsquare)

    closestPrior = priors[which.min(abs(priors - center))]

    return(c(center, closestPrior, ss))
}


calculateParameters1Cluster <- function(M, priormeans, numSNP, verbose) {
    Apriors = priormeans
   
    cl <- makeCluster(getNumberOfCores())
    clusterEvalQ(cl, library(crlmm))
    parAns <- parRapply(cl, M, calculateOneParams1Cluster, Apriors)
    stopCluster(cl)

    parAnsAllValues = unlist(parAns)
    parameters <- matrix(parAnsAllValues, numSNP, 3, byrow = TRUE)

    centers <- matrix(parameters[, 1], numSNP, 1)
    prior1cluster <- matrix(parameters[, 2], numSNP, 1)
    ss <- parameters[, 3]
    
    rm(parAns)
    rm(parameters)
    sigma=solve(cov(centers, use="complete.obs"))

    mh = calculateMahalDist1Cluster(centers, sigma, prior1cluster, numSNP)

    rm(sigma)
    rm(centers)

    gc()
    return(list(eless = ss, elemh = mh))
}


calculateMahalDist1Cluster <- function(centers, sigma, priors, numSNP){
    mahaldist = rep(NA, numSNP)
    
    for(i in 1:numSNP) {
        difference <- as.vector(centers[i, 1] - priors[i, 1])
	mahaldist[i]=difference%*%sigma%*%difference
    }
    return(mahaldist)
}
   
#############################################


computeCallPr <- function(callsPr, M, calls, numSNP, numSample, verbose, blockSize = 500000){
    # compute log-ratio in blocksize of 500,000 by default
    if (verbose) message("Start computing confidence score")
    
    numBlock = ceiling(numSNP / blockSize)
    for (i in 1:numBlock){
        if (verbose) message(" -- Processing segment ", i, " out of ", numBlock)
        subsetM = as.matrix(M[((i-1) * (blockSize) + 1):min(i * blockSize, numSNP), ])
        subsetCalls = as.matrix(calls[((i-1) * (blockSize) + 1):min(i * blockSize, numSNP), ])                     
        subsetCallProb = .Call("krlmmConfidenceScore", subsetM, subsetCalls, PACKAGE="crlmm");
        callsPr[((i-1) * (blockSize) + 1):min(i * blockSize, numSNP),] = subsetCallProb[1:nrow(subsetM), ]
        rm(subsetM, subsetCalls, subsetCallProb)
    }

    if (verbose) message("Done computing confidence score")
}

#############################################

AddBackNoVarianceSNPs <- function(cnSet, calls, scores, numSNP, numSample, variantSNPIndex, noVariantSNPIndex){
    callsGt = calls(cnSet)
    callsPr = snpCallProbability(cnSet)
    open(callsGt)
    open(callsPr)

    callsGt[variantSNPIndex, ] = calls[,]
    callsPr[variantSNPIndex, ] = scores[,]  
    
    callsGt[noVariantSNPIndex, ] = NA
    callsPr[noVariantSNPIndex, ] = 0     

    close(callsGt)
    close(callsPr)    
}

#############################################

krlmmImputeGender <- function(cnSet, annotation, YIndex, verbose){
    if (verbose) message("Start imputing gender")    
    S = computeAverageLogIntensity(cnSet, verbose) 

    # S is calculated and saved in original SNP order. 
    matchy = match(annotation[YIndex, 2], rownames(S))
    matchy = matchy[!is.na(matchy)]
    if (length(matchy) <= 10){
        predictedGender = rep(NA, ncol(A))
    }
    Sy = S[matchy,]

    uniqueDataPoint = apply(Sy, 1, function(x){length(unique(x))})
    validYSNPs = uniqueDataPoint >= 2
    SyValid = Sy[validYSNPs, ]
    
    if (nrow(SyValid) < 20){
        message("Not enough ChrY SNPs, skipping gender prediction step");
        predictedGender = rep(NA, ncol(Sy))
    }
   
    rm(S)
    numYChrSNP = nrow(SyValid)

    allS = matrix(NA, numYChrSNP, 2)

    for (i in 1:numYChrSNP) {
        tmp = kmeans(SyValid[i,] ,2, nstart=45)
        allS[i,] = sort(tmp$centers, decreasing=F)
    }
    priorS = apply(allS, 2, FUN="median", na.rm=TRUE)

    if (abs(priorS[1] - priorS[2]) <= 1.6) {
        message("Skipping gender prediction step");
        predictedGender = rep(NA, ncol(Sy))        
    }
    
    meanmatrix = apply(Sy, 2, median)

    Sy1 = abs(meanmatrix - priorS[1])
    Sy2 = abs(meanmatrix - priorS[2])

    # output male - 1, female - 2, female S-values are smaller than male S-values. 
    test = cbind(Sy2, Sy1)
    predictedGender = apply(test, 1, which.min)

    open(cnSet$gender)
    cnSet$gender[,] = predictedGender
    close(cnSet$gender)
    
    if (verbose) message("Done imputing gender")
    return(predictedGender)   
}


#############################################

verifyChrYSNPs <- function(cnSet, M, calls, gender, annotation, YIndex, priormeans, verbose){
    if (verbose) message("Start verifying SNPs on Chromosome Y")
    callsGt = calls(cnSet)
    callsPr = snpCallProbability(cnSet)
    open(callsGt)
    open(callsPr)
       
    matchy = match(annotation[YIndex, 2], rownames(M))
    matchy = matchy[!is.na(matchy)]
   
    MChrY = M[matchy,]
    callsChrY = calls[matchy,]
  
    male = gender == 1
    female = gender == 2    
    
    checkK = apply(callsChrY[, male], 1, function(x){ length(unique(x[!is.na(x)])) } )
    
    for(i in 1:nrow(MChrY)){
        # Chromosome Y SNPs, no calls for female, k = 1 or 2 permitted for male samples
        thisChrYSNPorigPosition = match(rownames(callsChrY)[i], rownames(callsGt))
        callsGt[thisChrYSNPorigPosition, female] = NA
        callsPr[thisChrYSNPorigPosition, female] = 0

        # early exit for k = 1 or 2 cases. Only re-assign calls to male samples if we previouly assigned
        # male samples to 3 clusters by mistake. 
        if (checkK[i] < 3) next;
          
        if (class(try(kmeans(MChrY[i, male], c(priormeans[1], priormeans[3]), nstart=1), TRUE)) != "try-error"){
           
            maleSampleCalls = kmeans(MChrY[i, male],c(priormeans[1], priormeans[3]), nstart=45)$cluster
            callsGt[thisChrYSNPorigPosition, male][maleSampleCalls == 1] = 1
            callsGt[thisChrYSNPorigPosition, male][maleSampleCalls == 2] = 3
        } else {
                    
            distanceToPrior1 = mean(abs(MChrY[i, male] - priormeans[1]))
            distanceToPrior3 = mean(abs(MChrY[i, male] - priormeans[3]))                
            callsGt[thisChrYSNPorigPosition, male] = ifelse(distanceToPrior1 < distanceToPrior3, 1, 3)
        }

        MMaleSamples = MChrY[i, male]
        callsMaleSample = callsGt[thisChrYSNPorigPosition, male]
        scoresMaleSample = .Call("krlmmConfidenceScore", t(as.matrix(MMaleSamples)), t(as.matrix(callsMaleSample)), PACKAGE="crlmm");

        callsPr[thisChrYSNPorigPosition, male] = scoresMaleSample       
    }

    close(callsGt)
    close(callsPr)
    
    if (verbose) message("Done verifying SNPs on Chromosome Y")
}
