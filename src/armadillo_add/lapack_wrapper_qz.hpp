/*
 * Add a QZ Decomposition to the Armadillo C++ Matrix Algebra Library.
 *
 * Keith O'Hara
 * 29/07/15
 */


#ifdef ARMA_USE_LAPACK

//! \namespace lapack namespace for LAPACK functions
namespace lapack
  {
      
  // qz
     
  template<typename eT>
  inline
  void
  gges
    (
    char* jobvsl, char* jobvsr, char* sort, char* selctg, blas_int* n,
    eT* a, blas_int* lda, eT* b, blas_int* ldb, blas_int* sdim,
    eT* alphar, eT* alphai, eT* beta,
    eT* vsl, blas_int* ldvsl, eT* vsr, blas_int* ldvsr,
    eT* work, blas_int* lwork, eT* bwork,
    blas_int* info
    )
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
          
    if(is_float<eT>::value == true)
    {
      typedef float T;
      arma_fortran(arma_sgges)(jobvsl, jobvsr, sort, selctg, n, (T*)a, lda, (T*)b, ldb, sdim, (T*)alphar, (T*)alphai, (T*)beta, (T*)vsl, ldvsl, (T*)vsr, ldvsr, (T*)work, lwork, (T*)bwork, info);
    }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dgges)(jobvsl, jobvsr, sort, selctg, n, (T*)a, lda, (T*)b, ldb, sdim, (T*)alphar, (T*)alphai, (T*)beta, (T*)vsl, ldvsl, (T*)vsr, ldvsr, (T*)work, lwork, (T*)bwork, info);
      }
    }
   

      
  template<typename eT>
  inline
  void
  cx_gges
    (
    char* jobvsl, char* jobvsr, char* sort, char* selctg, blas_int* n,
    eT* a, blas_int* lda, eT* b, blas_int* ldb, blas_int* sdim,
    eT* alpha, eT* beta,
    eT* vsl, blas_int* ldvsl, eT* vsr, blas_int* ldvsr,
    eT* work, blas_int* lwork, typename eT::value_type* rwork,
    typename eT::value_type* bwork,
    blas_int* info
    )
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
          
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef float T;
      typedef typename std::complex<T> cx_T;
      arma_fortran(arma_cgges)(jobvsl, jobvsr, sort, selctg, n, (cx_T*)a, lda, (cx_T*)b, ldb, sdim, (cx_T*)alpha, (cx_T*)beta, (cx_T*)vsl, ldvsl, (cx_T*)vsr, ldvsr, (cx_T*)work, lwork, (T*)rwork, (T*)bwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef double T;
      typedef typename std::complex<T> cx_T;
      arma_fortran(arma_zgges)(jobvsl, jobvsr, sort, selctg, n, (cx_T*)a, lda, (cx_T*)b, ldb, sdim, (cx_T*)alpha, (cx_T*)beta, (cx_T*)vsl, ldvsl, (cx_T*)vsr, ldvsr, (cx_T*)work, lwork, (T*)rwork, (T*)bwork, info);
      }
    }
  
  }

#endif
