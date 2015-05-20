test_crlmm <- function(){
	library(oligoClasses)
	if (require(genomewidesnp6Crlmm) & require(hapmapsnp6)){
		path <- system.file("celFiles", package="hapmapsnp6")
		cels <- list.celfiles(path, full.names=TRUE)
		(crlmmOutput <- crlmm(cels))
	}
}

test_duplicates <- function(){
	library(crlmm);library(RUnit); library(oligoClasses)
	if (require(genomewidesnp6Crlmm) & require(hapmapsnp6)){
		path <- system.file("celFiles", package="hapmapsnp6")
		cels <- list.celfiles(path, full.names=TRUE)
		## Quickly fails because sample identifiers are not unique
		checkException(crlmmOutput <- crlmm(cels[c(1,1, 2)]), silent=FALSE)

	}
	if(FALSE){
		library2(ff)
		datadir <- "/thumper/ctsa/snpmicroarray/illumina/IDATS/370k"
		## read in your samplesheet
		samplesheet <- read.csv(file.path(datadir, "HumanHap370Duo_Sample_Map.csv"), header=TRUE, as.is=TRUE)
		samplesheet <- samplesheet[-c(28:46,61:75,78:79), ]
		arrayNames <- file.path(datadir, unique(samplesheet[, "SentrixPosition"]))
		index <- c(1,1,2,2,3,4)
		any(duplicate(arrayNames[index]))
		arrayNames <- arrayNames[index]
		ss <- samplesheet[index, ]
		arrayInfo <- list(barcode=NULL, position="SentrixPosition")
		cnSet <- genotype.Illumina(sampleSheet=ss,
					   arrayNames=arrayNames,
					   arrayInfoColNames=arrayInfo,
					   cdfName="human370v1c",
					   batch=rep("1", nrow(ss)))
		checkTrue(validObject(cnSet))
	}

}
