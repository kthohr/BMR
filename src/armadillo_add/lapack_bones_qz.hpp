/*
 * Add a QZ Decomposition to the Armadillo C++ Matrix Algebra Library.
 *
 * Keith O'Hara
 * 29/07/15
 */

#ifdef ARMA_USE_LAPACK


#if !defined(ARMA_BLAS_CAPITALS)

  #define arma_sgges  sgges
  #define arma_dgges  dgges

  #define arma_cgges  cgges
  #define arma_zgges  zgges

#else

  #define arma_sgges  SGGES
  #define arma_dgges  DGGES

  #define arma_cgges  CGGES
  #define arma_zgges  ZGGES

#endif



extern "C"
  {
  // QZ decomposition of general real matrices
  void arma_fortran(arma_sgges)(char* jobvsl, char* jobvsr, char* sort, char* selctg, blas_int* n, float* a, blas_int* lda, float* b, blas_int* ldb, blas_int* sdim, float* alphar, float* alphai, float* beta, float* vsl, blas_int* ldvsl, float* vsr, blas_int* ldvsr, float* work, blas_int* lwork, float* bwork, blas_int* info);
  void arma_fortran(arma_dgges)(char* jobvsl, char* jobvsr, char* sort, char* selctg, blas_int* n, double* a, blas_int* lda, double* b, blas_int* ldb, blas_int* sdim, double* alphar, double* alphai, double* beta, double* vsl, blas_int* ldvsl, double* vsr, blas_int* ldvsr, double* work, blas_int* lwork, double* bwork, blas_int* info);
      
      // QZ decomposition of general complex matrices
  void arma_fortran(arma_cgges)(char* jobvsl, char* jobvsr, char* sort, char* selctg, blas_int* n, void* a, blas_int* lda, void* b, blas_int* ldb, blas_int* sdim, void* alpha, void* beta, void* vsl, blas_int* ldvsl, void* vsr, blas_int* ldvsr, void* work, blas_int* lwork,  float* rwork, float* bwork, blas_int* info);
  void arma_fortran(arma_zgges)(char* jobvsl, char* jobvsr, char* sort, char* selctg, blas_int* n, void* a, blas_int* lda, void* b, blas_int* ldb, blas_int* sdim, void* alpha, void* beta, void* vsl, blas_int* ldvsl, void* vsr, blas_int* ldvsr, void* work, blas_int* lwork, double* rwork, double* bwork, blas_int* info);

  }

#endif