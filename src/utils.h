int genotypeConfidence(const double *prob);
int genotypeConfidence2(double probability);
int intInSet(const int *x, const int *set, const int *n);
int genotypeCall(const double *pAA, const double *pAB, const double *pBB);
double sdCorrection(const int *n);
int sort_double(const double *a1, const double *a2);
void trimmed_mean(double *datavec, int *classvec, int class, double trim, int cols, int rows, double *m1, double *m2, double *m3, int i_ext);
void trimmed_stats(double *data, double *m1, double *m2, double *m3, int *class, int rows, int cols, double *trim);
long Cmatrix(int row, int col, int totrow);
