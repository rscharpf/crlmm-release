#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include "utils.h"

static double mydt(double x, int df){
  return(pow(1.0+pow(x, 2.0)/ (double)df, -((double)df+1.0)/2.0));
}

SEXP gtypeCallerPart1(SEXP A, SEXP B, SEXP fIndex, SEXP mIndex,
		      SEXP theCenters, SEXP theScales, SEXP theNs,
		      SEXP Indexes, SEXP cIndexes, SEXP nIndexes,
		      SEXP ncIndexes, SEXP SMEDIAN,
		      SEXP knots, SEXP mixtureParams, SEXP df,
		      SEXP probs, SEXP trim){
  /*
    ARGUMENTS
    ---------
    A: intensity matrix for allele A
    B: intensity matrix for allele B
    fIndex: indexes for females (columns in A/B for females)
    mIndex: indexes for males (columns in A/B for males)
    theCenters: matrix with SNP-specific centers (3 columns: AA/AB/BB)
    theScales: matrix with SNP-specific scales (3 columns: AA/AB/BB)
    theNs: matrix with SNP-specific counts (3 columns: AA/AB/BB)
    Indexes: list with 3 elements (autosomeIndex, XIndex, YIndex) for SNPs
    cIndexes: list with 3 elements (keepIndex, keepIndexFemale, keepIndexMale) for arrays
    SMEDIAN: scalar (median S)
    knots: knots for mixture
    mixtureParams: mixture parameters
    probs: genotype priors (1/3) for *ALL* SNPs. It's a vector of length 3
    trim: drop rate to estimate means

    ASSUMPTIONS
    -----------
    A and B have the same dimensions
    fIndex and mIndex are in a valid range for A/B
    Number of rows of theCenters, theScales, theNs match the number of rows in A/B
    Length of Indexes and cIndexes is 3
    3 knots
    4xNSAMPLES parameters
    priors has length 3 and is the same for *ALL* SNPs
    trim in (0, .5)

    INTERNAL VARIABLES
    -------- ---------
    likelihood: matrix (nsample rows x 3 columns)
    rowsAB, colsAB: dimensions of A and B
    lenLists = 3, length of Indexes and cIndexes
    h, i: iteration
    nIndex: length of Indexes[[h]]
    ii: particular value of Indexes
    M: log-ratio for a particular SNP
    S: adjusted average log-intensity for a particular SNP
    f: f for a particular SNP

    TODO
    ----
    - Get length of a vector from a list within C (adding nIndexes and ncIndexes for the moment)
  */

  /*
    ==========================================
              VARIABLE DECLARATION
    ==========================================
  */
  // Organizing variables
  int nSnpClasses, k;
  nSnpClasses = GET_LENGTH(Indexes);
  int nSnpsPerClass[nSnpClasses];
  for (k=0; k < nSnpClasses; k++)
    nSnpsPerClass[k] = GET_LENGTH(VECTOR_ELT(Indexes, k));
				  
  // General Variables
  int rowsAB, colsAB, h, i, ii, j, elem, nMales, nFemales;
  rowsAB = INTEGER(getAttrib(A, R_DimSymbol))[0];
  colsAB = INTEGER(getAttrib(A, R_DimSymbol))[1];
  double likelihood[colsAB*3], M[colsAB], S[colsAB], f[colsAB];

  // Constants
  //const int lenLists=3;

  // Buffers
  int intbuffer, ibv1[colsAB], ib2;
  double buffer;

  // All pointers appear here
  int *ptr2nIndexes, *ptr2A, *ptr2B, *ptr2Ns, *ptr2df, *ptr2mIndex, *ptr2fIndex; //*ptr2ncIndexes;
  double *ptr2Smedian, *ptr2knots, *ptr2params, *ptr2Centers, *ptr2Scales, *ptr2probs, *ptr2trim;
  ptr2nIndexes = INTEGER_POINTER(AS_INTEGER(nIndexes));
  ptr2A = INTEGER_POINTER(AS_INTEGER(A));
  ptr2B = INTEGER_POINTER(AS_INTEGER(B));
  ptr2Ns = INTEGER_POINTER(AS_INTEGER(theNs));
  ptr2df = INTEGER_POINTER(AS_INTEGER(df));
  ptr2mIndex = INTEGER_POINTER(AS_INTEGER(mIndex));
  ptr2fIndex = INTEGER_POINTER(AS_INTEGER(fIndex));
 // ptr2ncIndexes = INTEGER_POINTER(AS_INTEGER(ncIndexes));
  ptr2Smedian = NUMERIC_POINTER(AS_NUMERIC(SMEDIAN));
  ptr2knots = NUMERIC_POINTER(AS_NUMERIC(knots));
  ptr2params = NUMERIC_POINTER(AS_NUMERIC(mixtureParams));
  ptr2Centers = NUMERIC_POINTER(AS_NUMERIC(theCenters));
  ptr2Scales = NUMERIC_POINTER(AS_NUMERIC(theScales));
  ptr2probs = NUMERIC_POINTER(AS_NUMERIC(probs));
  ptr2trim = NUMERIC_POINTER(AS_NUMERIC(trim));
  // End pointers

  // These will be returned to R
  double *ptr2e1, *ptr2e2, *ptr2e3;
  SEXP estimates1, estimates2, estimates3, output;
  PROTECT(estimates1 = allocMatrix(REALSXP, rowsAB, 3));
  PROTECT(estimates2 = allocMatrix(REALSXP, rowsAB, 3));
  PROTECT(estimates3 = allocMatrix(REALSXP, rowsAB, 3));
  
  ptr2e1 = NUMERIC_POINTER(estimates1);
  ptr2e2 = NUMERIC_POINTER(estimates2);
  ptr2e3 = NUMERIC_POINTER(estimates3);
  /*
    ==========================================
            END VARIABLE DECLARATION
    ==========================================
  */
  nMales = GET_LENGTH(mIndex);
  nFemales = GET_LENGTH(fIndex);

  for (h=0; h < nSnpClasses; h++){
    ib2 = GET_LENGTH(VECTOR_ELT(cIndexes, h));
    double dbv[ib2];
    int ibv[ib2];
    if (nSnpsPerClass[h] > 0)
      for (i=0; i < ptr2nIndexes[h]; i++){
	/* if (i%100000 == 0) Rprintf("+");*/
	// Substract 1, as coming from R it is 1-based and C is 0-based.
	ii=INTEGER(VECTOR_ELT(Indexes, h))[i] - 1;
	for (j=0; j < colsAB; j++){
	  // j is an index for vectors whose length is number of samples (SAMPLE)
	  // elem is an index for A and B **only** (or objs with SNP rows and SAMPLE columns)
	  //Rprintf("J %d I %d Rows %d\n", j, ii, rowsAB);

	  elem = rowsAB * j + ii;
	  //Rprintf("\nElemt %d ", elem);
	  
	  M[j] = (log2(ptr2A[elem])-log2(ptr2B[elem]));
	  S[j] = (log2(ptr2A[elem])+log2(ptr2B[elem]))/2 - ptr2Smedian[0];
	  buffer = fmax(fmin(S[j], ptr2knots[2]), ptr2knots[0]);
	  f[j] = ptr2params[j*4+0]+ptr2params[j*4+1]*buffer+ptr2params[j*4+2]*pow(buffer, 2.0)+ptr2params[j*4+3]*pow(fmax(0, buffer-ptr2knots[1]), 2.0);
	  //Rprintf("M %f S %f f %f ", M[j], S[j], f[j]);

	  // buffer here is sigma
	  // likelihood for AA
	  // All likelihoods already multiplied by prior to save time
	  buffer = ptr2Scales[ii] * sdCorrection(&ptr2Ns[ii]);
	  likelihood[j] = mydt( ((M[j]-f[j])-ptr2Centers[ii])/buffer, ptr2df[0])*ptr2probs[0];

	  //Rprintf("L1 %2.4f ", likelihood[j]);

	  // likelihood for AB
	  buffer = ptr2Scales[ii+rowsAB] * sdCorrection(&ptr2Ns[ii+rowsAB]);
	  likelihood[j+colsAB] = mydt( (M[j]-ptr2Centers[ii+rowsAB])/buffer, ptr2df[0])*ptr2probs[1];
	  
	  // intbuffer (here) is the subject ID as in R (ie. 1-based)
	  intbuffer = j+1;
	  if (nMales > 0)
	    if (h > 0)
	      if (intInSet(&intbuffer, ptr2mIndex, &nMales) > 0)
		likelihood[j+colsAB] = 0;

	  //Rprintf("L2 %2.4f ", likelihood[j+colsAB]);


	  // likelihood for BB
	  buffer = ptr2Scales[ii+2*rowsAB] * sdCorrection(&ptr2Ns[ii+2*rowsAB]);
	  likelihood[j+2*colsAB] = mydt( ((M[j]+f[j])-ptr2Centers[ii+2*rowsAB])/buffer, ptr2df[0])*ptr2probs[2];
	  
	  // Females on Y: 1 to avoid NAs. Later made 0 (RI)
	  // To save some time: 1*priors = priors
	  if (nFemales > 0)
	    if (h == 2)
	      if (intInSet(&intbuffer, ptr2fIndex, &nFemales) >0){
		likelihood[j] = ptr2probs[2];
		likelihood[j+colsAB] = ptr2probs[2];
		likelihood[j+2*colsAB] = ptr2probs[2];
	      }

	  //Rprintf("L3 %2.4f ", likelihood[j+2*colsAB]);


	  // Compute simple posterior
	  buffer = likelihood[j]+likelihood[j+colsAB]+likelihood[j+2*colsAB];
	  likelihood[j]/=buffer;
	  likelihood[j+colsAB]/=buffer;
	  likelihood[j+2*colsAB]/=buffer;

	  if (nFemales > 0)
	    if (h == 2)
	      if (intInSet(&intbuffer, ptr2fIndex, &nFemales) >0){
		likelihood[j] = 0;
		likelihood[j+colsAB] = 0;
		likelihood[j+2*colsAB] = 0;
	      }

	  ibv1[j] = genotypeCall(&likelihood[j], &likelihood[j+colsAB], &likelihood[j+2*colsAB]);
	}
      
	for (j=0; j < ib2; j++){
	  intbuffer = INTEGER(VECTOR_ELT(cIndexes, h))[j]-1;
	  dbv[j]=M[intbuffer]-f[intbuffer];
	  ibv[j]=ibv1[intbuffer];
	}
	trimmed_mean(dbv, ibv, 1, ptr2trim[0], GET_LENGTH(VECTOR_ELT(cIndexes, h)), rowsAB, ptr2e1, ptr2e2, ptr2e3, ii);
	for (j=0; j < ib2; j++){
	  intbuffer = INTEGER(VECTOR_ELT(cIndexes, h))[j]-1;
	  dbv[j]=M[intbuffer];
	}
	trimmed_mean(dbv, ibv, 2, ptr2trim[0], GET_LENGTH(VECTOR_ELT(cIndexes, h)), rowsAB, ptr2e1, ptr2e2, ptr2e3, ii);
	for (j=0; j < ib2; j++){
	  intbuffer = INTEGER(VECTOR_ELT(cIndexes, h))[j]-1;
	  dbv[j] = M[intbuffer]+f[intbuffer];
	}
	trimmed_mean(dbv, ibv, 3, ptr2trim[0], GET_LENGTH(VECTOR_ELT(cIndexes, h)), rowsAB, ptr2e1, ptr2e2, ptr2e3, ii);
      } /* for Snp */
  } /* for SnpClass */
  PROTECT(output = allocVector(VECSXP,3));
  SET_VECTOR_ELT(output, 0, estimates1);
  SET_VECTOR_ELT(output, 1, estimates2);
  SET_VECTOR_ELT(output, 2, estimates3);
  UNPROTECT(4);
  return(output);
}

SEXP gtypeCallerPart2(SEXP A, SEXP B, SEXP fIndex, SEXP mIndex,
		      SEXP theCenters, SEXP theScales, SEXP theNs,
		      SEXP Indexes, SEXP cIndexes, SEXP nIndexes,
		      SEXP ncIndexes, SEXP SMEDIAN,
		      SEXP knots, SEXP mixtureParams, SEXP df,
		      SEXP probs, SEXP trim, SEXP noTraining, SEXP noInfo){
  /*
    WARNING!!! REMEMBER TO MODIFY MY TWIN TOO!

    ARGUMENTS
    ---------
    A: intensity matrix for allele A
    B: intensity matrix for allele B
    fIndex: indexes for females (columns in A/B for females)
    mIndex: indexes for males (columns in A/B for males)
    theCenters: matrix with SNP-specific centers (3 columns: AA/AB/BB)
    theScales: matrix with SNP-specific scales (3 columns: AA/AB/BB)
    theNs: matrix with SNP-specific counts (3 columns: AA/AB/BB)
    Indexes: list with 3 elements (autosomeIndex, XIndex, YIndex) for SNPs
    cIndexes: list with 3 elements (keepIndex, keepIndexFemale, keepIndexMale) for arrays
    SMEDIAN: scalar (median S)
    knots: knots for mixture
    mixtureParams: mixture parameters
    probs: genotype priors (1/3) for *ALL* SNPs. It's a vector of length 3
    trim: drop rate to estimate means

    ASSUMPTIONS
    -----------
    A and B have the same dimensions
    fIndex and mIndex are in a valid range for A/B
    Number of rows of theCenters, theScales, theNs match the number of rows in A/B
    Length of Indexes and cIndexes is 3
    3 knots
    4xNSAMPLES parameters
    priors has length 3 and is the same for *ALL* SNPs
    trim in (0, .5)
    Indexes and cIndexes have the same length.

    INTERNAL VARIABLES
    -------- ---------
    likelihood: matrix (nsample rows x 3 columns)
    rowsAB, colsAB: dimensions of A and B
    lenLists = 3, length of Indexes and cIndexes
    h, i: iteration
    nIndex: length of Indexes[[h]]
    ii: particular value of Indexes
    M: log-ratio for a particular SNP
    S: adjusted average log-intensity for a particular SNP
    f: f for a particular SNP

    TODO
    ----
    - Get length of a vector from a list within C (adding nIndexes and ncIndexes for the moment)
  */

  /*
    ==========================================
              VARIABLE DECLARATION
    ==========================================
  */
  // Organizing variables
  int nSnpClasses, k;
  nSnpClasses = GET_LENGTH(Indexes);
  int nSnpsPerClass[nSnpClasses];
  for (k=0; k < nSnpClasses; k++)
    nSnpsPerClass[k] = GET_LENGTH(VECTOR_ELT(Indexes, k));

  // General Variables
  int rowsAB, colsAB, h, i, ii, j, elem, nMales, nFemales;
  rowsAB = INTEGER(getAttrib(A, R_DimSymbol))[0];
  colsAB = INTEGER(getAttrib(A, R_DimSymbol))[1];
  double likelihood[colsAB*3], M[colsAB], S[colsAB], f[colsAB];

  // Buffers
  int intbuffer, ib2, ib3, ibSnpLevel1=0, ibSnpLevel2=0;
  double buffer;

  ib2 = GET_LENGTH(noTraining);
  ib3 = GET_LENGTH(noInfo);

  // All pointers appear here
  int *ptr2nIndexes, *ptr2A, *ptr2B, *ptr2Ns, *ptr2df, *ptr2mIndex, *ptr2fIndex; // *ptr2ncIndexes, 
  int *ptr2noTraining, *ptr2noInfo;
  double *ptr2Smedian, *ptr2knots, *ptr2params, *ptr2Centers, *ptr2Scales, *ptr2probs; // *ptr2trim;
  ptr2nIndexes = INTEGER_POINTER(AS_INTEGER(nIndexes));
  ptr2A = INTEGER_POINTER(AS_INTEGER(A));
  ptr2B = INTEGER_POINTER(AS_INTEGER(B));
  ptr2Ns = INTEGER_POINTER(AS_INTEGER(theNs));
  ptr2df = INTEGER_POINTER(AS_INTEGER(df));
  ptr2mIndex = INTEGER_POINTER(AS_INTEGER(mIndex));
  ptr2fIndex = INTEGER_POINTER(AS_INTEGER(fIndex));
 // ptr2ncIndexes = INTEGER_POINTER(AS_INTEGER(ncIndexes));
  ptr2noTraining = INTEGER_POINTER(AS_INTEGER(noTraining));
  ptr2noInfo = INTEGER_POINTER(AS_INTEGER(noInfo));

  ptr2Smedian = NUMERIC_POINTER(AS_NUMERIC(SMEDIAN));
  ptr2knots = NUMERIC_POINTER(AS_NUMERIC(knots));
  ptr2params = NUMERIC_POINTER(AS_NUMERIC(mixtureParams));
  ptr2Centers = NUMERIC_POINTER(AS_NUMERIC(theCenters));
  ptr2Scales = NUMERIC_POINTER(AS_NUMERIC(theScales));
  ptr2probs = NUMERIC_POINTER(AS_NUMERIC(probs));
  // ptr2trim = NUMERIC_POINTER(AS_NUMERIC(trim));

  // End pointers

  /*
    ==========================================
            END VARIABLE DECLARATION
    ==========================================
  */
  nMales = GET_LENGTH(mIndex);
  nFemales = GET_LENGTH(fIndex);

  for (h=0; h < nSnpClasses; h++){
    if (nSnpsPerClass[h] > 0)
      for (i=0; i < ptr2nIndexes[h]; i++){
	/* if (i%100000 == 0) Rprintf("+"); */
	// Substract 1, as coming from R it is 1-based and C is 0-based.
	ii=INTEGER(VECTOR_ELT(Indexes, h))[i] - 1;
	intbuffer = ii+1;
	if (intInSet(&intbuffer, ptr2noTraining, &ib2) > 0) ibSnpLevel1 = 1;
	if (intInSet(&intbuffer, ptr2noInfo, &ib3) > 0) ibSnpLevel2 = 1;
	for (j=0; j < colsAB; j++){
	  // j is an index for vectors whose length is number of samples (SAMPLE)
	  // elem is an index for A and B **only** (or objs with SNP rows and SAMPLE columns)
	  elem = rowsAB * j + ii;
	  M[j] = (log2(ptr2A[elem])-log2(ptr2B[elem]));
	  S[j] = (log2(ptr2A[elem])+log2(ptr2B[elem]))/2 - ptr2Smedian[0];
	  buffer = fmax(fmin(S[j], ptr2knots[2]), ptr2knots[0]);
	  f[j] = ptr2params[j*4+0]+ptr2params[j*4+1]*buffer+ptr2params[j*4+2]*pow(buffer, 2.0)+ptr2params[j*4+3]*pow(fmax(0, buffer-ptr2knots[1]), 2.0);
	  
	  // buffer here is sigma
	  // likelihood for AA
	  // All likelihoods already multiplied by prior to save time
	  buffer = ptr2Scales[ii] * sdCorrection(&ptr2Ns[ii]);
	  likelihood[j] = mydt( ((M[j]-f[j])-ptr2Centers[ii])/buffer, ptr2df[0])*ptr2probs[0];

	  // likelihood for AB
	  buffer = ptr2Scales[ii+rowsAB] * sdCorrection(&ptr2Ns[ii+rowsAB]);
	  likelihood[j+colsAB] = mydt( (M[j]-ptr2Centers[ii+rowsAB])/buffer, ptr2df[0])*ptr2probs[1];
	  
	  // intbuffer (here) is the subject ID as in R (ie. 1-based)
	  // h > 0 is chr X or Y
	  intbuffer = j+1;
	  if (nMales > 0)
	    if (h > 0 && intInSet(&intbuffer, ptr2mIndex, &nMales) > 0) likelihood[j+colsAB] = 0;

	  // likelihood for BB
	  buffer = ptr2Scales[ii+2*rowsAB] * sdCorrection(&ptr2Ns[ii+2*rowsAB]);
	  likelihood[j+2*colsAB] = mydt( ((M[j]+f[j])-ptr2Centers[ii+2*rowsAB])/buffer, ptr2df[0])*ptr2probs[2];

	  // Females on Y: 1 to avoid NAs. Later made 0 (RI)
	  // To save some time: 1*priors = priors
	  if (nFemales > 0)
	    if (h == 2 && intInSet(&intbuffer, ptr2fIndex, &nFemales) > 0)
	      likelihood[j] = likelihood[j+colsAB] = likelihood[j+2*colsAB] = ptr2probs[2];
	  
	  // Compute simple posterior
	  buffer = likelihood[j]+likelihood[j+colsAB]+likelihood[j+2*colsAB];
	  likelihood[j]/=buffer;
	  likelihood[j+colsAB]/=buffer;
	  likelihood[j+2*colsAB]/=buffer;

	  if (nFemales > 0)
	    if (h == 2 && intInSet(&intbuffer, ptr2fIndex, &nFemales) > 0)
	      likelihood[j] = likelihood[j+colsAB] = likelihood[j+2*colsAB] = 0;

	  // IDENTICAL UNTIL HERE
	  ptr2A[elem] = genotypeCall(&likelihood[j], &likelihood[j+colsAB], &likelihood[j+2*colsAB]);
	  buffer = fmax(fmax(likelihood[j], likelihood[j+colsAB]), likelihood[j+2*colsAB]);
	  if (ibSnpLevel1 == 1) buffer *= 0.995;
	  if (ibSnpLevel2 == 1) buffer *= 0.000;
	  ptr2B[elem] = genotypeConfidence(&buffer);
	}
	ibSnpLevel1 = ibSnpLevel2 = 0;
      }
  }
  return(R_NilValue);
}


SEXP krlmmComputeM(SEXP A, SEXP B){

  /*
    ARGUMENTS
    ---------
    A: intensity matrix for allele A
    B: intensity matrix for allele B
    M: log-ratio for a particular SNP (outgoing)

    INTERNAL VARIABLES
    -------- ---------
    rowsAB, colsAB: dimensions of A and B, which is number of SNP and number of sample
  */

  int rowsAB, colsAB;
  rowsAB = INTEGER(getAttrib(A, R_DimSymbol))[0];
  colsAB = INTEGER(getAttrib(A, R_DimSymbol))[1];

  int i, j;

  int *ptr2A, *ptr2B;
  ptr2A = INTEGER_POINTER(AS_INTEGER(A));
  ptr2B = INTEGER_POINTER(AS_INTEGER(B));

  SEXP Rval;
  PROTECT(Rval = allocMatrix(REALSXP, rowsAB, colsAB));

  double *ptr2M;
  long ele;
  ptr2M = NUMERIC_POINTER(Rval);

  for (i = 1; i <= rowsAB; i++){
    for (j = 1; j <= colsAB; j++){
      // elem is an index for A, B and M
      ele = Cmatrix(i, j, rowsAB);
      ptr2M[ele] = (log2(ptr2A[ele])-log2(ptr2B[ele]));
    }
  }
  
  UNPROTECT(1);
  return Rval;
}


SEXP krlmmComputeS(SEXP A, SEXP B){

  /*
    ARGUMENTS
    ---------
    A: intensity matrix for allele A
    B: intensity matrix for allele B
    S: average log-intensity for a particular SNP (outgoing)

    INTERNAL VARIABLES
    -------- ---------
    rowsAB, colsAB: dimensions of A and B, which is number of SNP and number of sample
  */

  int rowsAB, colsAB;
  rowsAB = INTEGER(getAttrib(A, R_DimSymbol))[0];
  colsAB = INTEGER(getAttrib(A, R_DimSymbol))[1];

  int i, j;

  int *ptr2A, *ptr2B;
  ptr2A = INTEGER_POINTER(AS_INTEGER(A));
  ptr2B = INTEGER_POINTER(AS_INTEGER(B));

  SEXP Rval;
  PROTECT(Rval = allocMatrix(REALSXP, rowsAB, colsAB));

  double *ptr2S;
  long ele;
  ptr2S = NUMERIC_POINTER(Rval);

  for (i = 1; i <= rowsAB; i++){
    for (j = 1; j <= colsAB; j++){
      // elem is an index for A, B and S
      ele = Cmatrix(i, j, rowsAB);
      ptr2S[ele] = (log2(ptr2A[ele])+log2(ptr2B[ele]))/2;
    }
  }
  
  UNPROTECT(1);
  return Rval;
}

void calculate_multiple_cluster_scores(int row, double *intensity, double mean_intensity, int *clustering, int num_SNP, int num_sample, int *ptr, double **dist, int *clustercount){
    int p, q, o;
    long vectorPos;
    double a_dist;
    for(p=1; p <= num_sample; p++){
        dist[p][p] = 0;
    }
    for(p=1; p<num_sample; p++){
        for(q=p+1; q<=num_sample; q++){
   	    a_dist = fabs(intensity[Cmatrix(row, p, num_SNP)] - intensity[Cmatrix(row, q, num_SNP)]);
            dist[p][q] = a_dist;
            dist[q][p] = a_dist;
        }
    }

    double sum[4];
    double within_cluster;
    double between_cluster;
    double temp_between_cluster;
    double max;
    for (p=1; p <= num_sample; p++){
        vectorPos = Cmatrix(row, p, num_SNP);
        sum[1] = 0;
        sum[2] = 0;
        sum[3] = 0;
        for (q=1; q<=num_sample; q++){
	  sum[clustering[Cmatrix(row, q, num_SNP)]] += dist[p][q];
        }
       
        within_cluster = sum[clustering[vectorPos]] / (clustercount[clustering[vectorPos]]- 1);
        between_cluster = -1.0;
        for (o=1; o<=3; o++){
            if ((o != clustering[vectorPos]) && (clustercount[o] > 0)){
                temp_between_cluster = sum[o] / clustercount[o];
                if ((between_cluster < 0) || (temp_between_cluster < between_cluster)) {
                    between_cluster = temp_between_cluster;
                }                                           
            }                                               
        }
       
        if (clustercount[clustering[vectorPos]] > 1) {
            if (between_cluster > within_cluster) {
                max = between_cluster;
            } else {
                max = within_cluster;
            }
            ptr[vectorPos] = genotypeConfidence2((between_cluster - within_cluster) / max);                                                   
        }
    }   
}

void calculate_unique_cluster_scores(int row, double *intensity, double mean_intensity, double intensity_range, int num_SNP, int num_sample, int *ptr)				   
{
    int p;
    long vectorPos;

    for(p = 1; p <= num_sample; p++){
        vectorPos = Cmatrix(row, p, num_SNP); 
        ptr[vectorPos] = genotypeConfidence2(1 -  fabs(fabs(intensity[vectorPos] - mean_intensity) / intensity_range));
    }
}


double calculate_SNP_mean(int row, double *intensity, int num_SNP, int num_sample)
{
  double sum_intensity;
  sum_intensity = 0;
  int p;
  for(p = 1; p <= num_sample; p++){
    sum_intensity = sum_intensity + intensity[Cmatrix(row, p, num_SNP)];
  }
  return(sum_intensity / num_sample);
}

double calculate_SNP_range(int row, double *intensity, int num_SNP, int num_sample)
{
  double min_intensity;
  double max_intensity;
  max_intensity = intensity[Cmatrix(row, 1, num_SNP)];
  min_intensity = intensity[Cmatrix(row, 1, num_SNP)];
  int p;  
  for(p = 2; p <= num_sample; p++){
    if (intensity[Cmatrix(row, p, num_SNP)] < min_intensity){
      min_intensity = intensity[Cmatrix(row, p, num_SNP)];  
    }	
    if (intensity[Cmatrix(row, p, num_SNP)] > max_intensity){
      max_intensity = intensity[Cmatrix(row, p, num_SNP)];
    }
  }
  return(max_intensity - min_intensity);
}


SEXP krlmmConfidenceScore(SEXP M, SEXP clustering)
{
    int num_SNP, num_sample;
    num_SNP = INTEGER(getAttrib(M, R_DimSymbol))[0];
    num_sample = INTEGER(getAttrib(M, R_DimSymbol))[1];

    double *ptr2M;
    int *ptr2cluster;
    ptr2M = NUMERIC_POINTER(AS_NUMERIC(M));
    ptr2cluster = INTEGER_POINTER(AS_INTEGER(clustering));

    int i, j;
    int cluster;
    double mean_intensity;
    double intensity_range;
    int k;
    
    SEXP Rval;
    PROTECT(Rval = allocMatrix(INTSXP, num_SNP, num_sample));

    int *ptr2score;
    ptr2score = INTEGER_POINTER(Rval);
    
    double **dist;
    // allocate memory to dist
    dist = (double **)malloc((num_sample + 1)*sizeof(double *));
    for (i = 1; i <= num_sample; i++){
      dist[i] = (double *)malloc((num_sample + 1) * sizeof(double)); 
    }

    int cluster_count[4]; // extra element to cope with zero-base notation
   // int clid[4];

    for (i = 1; i <= num_SNP; i++) {
        cluster_count[1] = 0;
        cluster_count[2] = 0;
        cluster_count[3] = 0;
        for (j=1; j<= num_sample; j++){
	    cluster = ptr2cluster[Cmatrix(i, j, num_SNP)];
            cluster_count[cluster]++;
        }
	mean_intensity = calculate_SNP_mean(i, ptr2M, num_SNP, num_sample);
	intensity_range = calculate_SNP_range(i, ptr2M, num_SNP, num_sample);

        k = 0;
        for (j=1; j<=3; j++){
            if (cluster_count[j] > 0){
                k++;
            //    clid[k] = j;
            }
        }

	if (k==1) {
	  calculate_unique_cluster_scores(i, ptr2M, mean_intensity, intensity_range, num_SNP, num_sample, ptr2score);
	} else {
	  calculate_multiple_cluster_scores(i, ptr2M, mean_intensity, ptr2cluster, num_SNP, num_sample, ptr2score, dist, cluster_count);
	}       
      
    } 

    UNPROTECT(1);
    return Rval;
}

