/*
 * Add a QZ Decomposition to the Armadillo C++ Matrix Algebra Library.
 *
 * Keith O'Hara
 * 29/07/15
 */


//! \addtogroup auxlib
//! @{

//! QZ Decomposition of general square real matrix pair (A,B).
//! The argument 'side' specifies which of 'Q' (left) and 'Z'
//! (right) should be returned (or both 'b').
template<typename T, typename T1, typename T2>
inline
bool
auxlib2::qz
  (
  Mat<T>&             Smat,
  Mat<T>&             Tmat,
  Mat<T>&              vsl,
  Mat<T>&              vsr,
  const Base<T,T1>&   X,
  const Base<T,T2>&   Y,
  const char          side
  )
  {
  arma_extra_debug_sigprint();
    
  #if defined(ARMA_USE_LAPACK)
    {
    char jobvsl;
    char jobvsr;
        
    switch(side)
      {
      case 'l':  // left
        jobvsl = 'V';
        jobvsr = 'N';
        break;
                
      case 'r':  // right
        jobvsl = 'N';
        jobvsr = 'V';
        break;
                
      case 'b':  // both
        jobvsl = 'V';
        jobvsr = 'V';
        break;
              
      case 'n':  // neither
        jobvsl = 'N';
        jobvsr = 'N';
        break;
                
      default:
        arma_stop("qz(): parameter 'side' is invalid");
        return false;
      }
        
    char eigsort = 'N';
    char kdum = 'N';  // just a dummy
    blas_int sdim  = 0;
        
    Mat<T> A(X.get_ref());
    Mat<T> B(Y.get_ref());
        
    arma_debug_check( ((A.is_square() == false) || (B.is_square() == false)), "qz(): given matrix is not square" );
        
    arma_debug_check( (A.n_rows != B.n_rows), "qz(): given matrices must have the same size" );
        
    if(A.is_empty())
      {
      Smat.reset();
      Tmat.reset();
      vsl.reset();
      vsr.reset();
      return true;
      }
        
    arma_debug_assert_blas_size(A);
        
    const uword A_n_rows = A.n_rows;
        
    podarray<T> alphar(A_n_rows);
    podarray<T> alphai(A_n_rows);
    podarray<T>   beta(A_n_rows);
        
    vsl.set_size( ((jobvsl == 'V') ? A_n_rows : 1), A_n_rows );
    vsr.set_size( ((jobvsr == 'V') ? A_n_rows : 1), A_n_rows );
        
    blas_int N     = blas_int(A_n_rows);
    blas_int lwork = 3 * ((std::max)(blas_int(1),8*N+16));
    blas_int info  = 0;
        
    podarray<T>  work( static_cast<uword>(lwork) );
    podarray<T>  bwork( static_cast<uword>(N) );
        
    arma_extra_debug_print("lapack::gges()");
    lapack::gges
      (
      &jobvsl, &jobvsr, &eigsort, &kdum, &N,
      A.memptr(), &N, B.memptr(), &N, &sdim,
      alphar.memptr(), alphai.memptr(), beta.memptr(),
      vsl.memptr(), &N, vsr.memptr(), &N,
      work.memptr(), &lwork, bwork.memptr(),
      &info
      );
    
    op_strans::apply_mat_inplace(vsl); // transpose Q
        
      if(info == 0)
        {
        Smat = A;
        Tmat = B;
        return true;
        }
      else
        {
        return false;
        }
    }
  #else
    {
    arma_ignore(vsl);
    arma_ignore(vsr);
    arma_ignore(Smat);
    arma_ignore(Tmat);
    arma_ignore(X);
    arma_ignore(Y);
    arma_ignore(side);
    arma_stop("qz(): use of LAPACK needs to be enabled");
    return false;
    }
  #endif
  }

//! QZ Decomposition of general square complex matrix pair (A,B).
//! The argument 'side' specifies which of 'Q' (left) and 'Z'
//! (right) should be returned (or both 'b').
template<typename T, typename T1, typename T2>
inline
bool
auxlib2::qz
  (
  Mat< std::complex<T> >&            Smat,
  Mat< std::complex<T> >&            Tmat,
  Mat< std::complex<T> >&             vsl,
  Mat< std::complex<T> >&             vsr,
  const Base< std::complex<T>, T1 >& X,
  const Base< std::complex<T>, T2 >& Y,
  const char                         side
  )
  {
  arma_extra_debug_sigprint();
    
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename std::complex<T> eT;
        
    char jobvsl;
    char jobvsr;
        
    switch(side)
      {
      case 'l':  // left
        jobvsl = 'V';
        jobvsr = 'N';
        break;
                
      case 'r':  // right
        jobvsl = 'N';
        jobvsr = 'V';
        break;
                
      case 'b':  // both
        jobvsl = 'V';
        jobvsr = 'V';
        break;
                
      case 'n':  // neither
        jobvsl = 'N';
        jobvsr = 'N';
        break;
                
      default:
        arma_stop("qz(): parameter 'side' is invalid");
        return false;
      }
        
    char eigsort = 'N';
    char kdum = 'N';  // just a dummy
    blas_int sdim  = 0;
        
    Mat<eT> A(X.get_ref());
    Mat<eT> B(Y.get_ref());
        
    arma_debug_check( ((A.is_square() == false) || (B.is_square() == false)), "qz(): given matrix is not square" );
        
    arma_debug_check( (A.n_rows != B.n_rows), "qz(): given matrices must have the same size" );
        
    if(A.is_empty())
      {
      Smat.reset();
      Tmat.reset();
      vsl.reset();
      vsr.reset();
      return true;
      }
        
    arma_debug_assert_blas_size(A);
        
    const uword A_n_rows = A.n_rows;
        
    podarray<eT> alpha(A_n_rows);
    podarray<eT>  beta(A_n_rows);
        
    vsl.set_size( ((jobvsl == 'V') ? A_n_rows : 1), A_n_rows );
    vsr.set_size( ((jobvsr == 'V') ? A_n_rows : 1), A_n_rows );
        
    blas_int N     = blas_int(A_n_rows);
    blas_int lwork = 3 * ((std::max)(blas_int(1),2*N));
    blas_int info  = 0;
        
    podarray<eT>  work( static_cast<uword>(lwork) );
    podarray<T>  bwork( static_cast<uword>(N) );
    podarray<T>  rwork( static_cast<uword>(8*N)   );
        
    arma_extra_debug_print("lapack::cx_gges()");
    lapack::cx_gges
      (
      &jobvsl, &jobvsr, &eigsort, &kdum, &N,
      A.memptr(), &N, B.memptr(), &N, &sdim,
      alpha.memptr(), beta.memptr(),
      vsl.memptr(), &N, vsr.memptr(), &N,
      work.memptr(), &lwork, rwork.memptr(), bwork.memptr(),
      &info
      );
    
    op_htrans::apply_mat_inplace(vsl); // transpose Q
        
    if(info == 0)
      {
      Smat = A;
      Tmat = B;
      return true;
      }
    else
      {
      return false;
      }
    }
  #else
    {
    arma_ignore(vsl);
    arma_ignore(vsr);
    arma_ignore(Smat);
    arma_ignore(Tmat);
    arma_ignore(X);
    arma_ignore(Y);
    arma_ignore(side);
    arma_stop("qz(): use of LAPACK needs to be enabled");
    return false;
    }
  #endif
  }


//! @}
