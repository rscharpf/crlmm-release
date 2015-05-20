## Need to submit from a fat node!
all:	illumina affy 

affy:	copynumber.Rnw
	echo "Stangle(\"copynumber.Rnw\")" | R-devel --no-save --no-restore;	
	cat ~/bin/cluster.template | perl -pe "s/Rprog/copynumber.R/" > copynumber.R.sh
	qsub -m e -r y -cwd -l mem_free=12G,h_vmem=16G copynumber.R.sh

illumina:	illumina_copynumber.Rnw
	echo "Stangle(\"illumina_copynumber.Rnw\")" | R --no-save --no-restore;	
	cat ~/bin/cluster.templateRelease | perl -pe "s/Rprog/illumina_copynumber.R/" > illumina_copynumber.R.sh
	qsub -m e -r y -cwd -l mem_free=12G,h_vmem=16G illumina_copynumber.R.sh

texcopynumber:	copynumber.tex
	texi2dvi --pdf --clean --quiet copynumber.tex

texillumina_copynumber:	illumina_copynumber.tex
	texi2dvi --pdf --clean --quiet illumina_copynumber.tex		

