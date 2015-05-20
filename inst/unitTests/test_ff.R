genotypingTest <- function(){
	require(genomewidesnp6Crlmm) & require(hapmapsnp6)
	path <- system.file("celFiles", package="hapmapsnp6")
	cels <- list.celfiles(path, full.names=TRUE)
	ocProbesets(50e3)
	batch <- as.factor(rep("A", length(cels)))
	(cnSet <- genotype(cels, cdfName="genomewidesnp6", batch=batch))
	library(ff)
	ldPath(tempdir())
	(cnSet2 <- genotype(cels, cdfName="genomewidesnp6", batch=batch))
	checkTrue(all.equal(calls(cnSet), calls(cnSet2)[,]))
}

genotypingTestIllumina <- function(){
	setwd("/thumper/ctsa/snpmicroarray/illumina/IDATS/370k/")
	library(crlmm)
	library(ff)
	ldPath(tempdir())

	samples <- read.csv("samples370k.csv", as.is=TRUE)
	RG <- readIdatFiles(sampleSheet=samples,
			    arrayInfoColNames=list(barcode=NULL, position="SentrixPosition"),
			    saveDate=TRUE)

	crlmmResult <-  crlmmIllumina(RG=RG,
				      cdfName="human370v1c",
				      returnParams=TRUE)
	checkTrue(is(calls(crlmmResult)[1:5,1], "integer"))

	crlmmResult2 <- crlmmIlluminaV2(sampleSheet=samples,
					arrayInfoColNames=list(barcode=NULL, position="SentrixPosition"),
					saveDate=TRUE, cdfName="human370v1c", returnParams=TRUE)
	checkTrue(identical(calls(crlmmResult)[1:5, ]),
		  identical(calls(crlmmResult2)[1:5, ]))

	crlmmResult3 <- genotype.Illumina(sampleSheet=samples,
					  arrayInfoColNames=list(barcode=NULL,
					  position="SentrixPosition"),
					  saveDate=TRUE, cdfName="human370v1c",
					  batch = as.factor(rep(1, nrow(samples))))

	checkTrue(identical(calls(crlmmResult)[1:5, ]),
		  identical(calls(crlmmResult3)[1:5, ]))

}
