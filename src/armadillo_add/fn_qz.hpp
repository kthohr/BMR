/*
 * Add a QZ Decomposition to the Armadillo C++ Matrix Algebra Library.
 *
 * Keith O'Hara
 * 29/07/15
 */


//! \addtogroup fn_qz
//! @{


//! QZ decomposition for pair of N-by-N general real matrices (A,B)
template<typename T, typename T1, typename T2>
inline
bool
qz
  (
  Mat<T>&             AA,
  Mat<T>&             BB,
  Mat<T>&              Q,
  Mat<T>&              Z,
  const Base<T,T1>&   A,
  const Base<T,T2>&   B,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
    
  const bool status = auxlib2::qz(AA, BB, Q, Z, A, B, 'b');
    
  if(status == false)
    {
    AA.reset();
    BB.reset();
    Q.reset();
    Z.reset();
    arma_bad("qz(): failed to converge", false);
    }
    
    return status;
  }

//! QZ decomposition for pair of N-by-N general complex matrices (A,B)
template<typename T, typename T1, typename T2>
inline
bool
qz
  (
  Mat< std::complex<T> >&            AA,
  Mat< std::complex<T> >&            BB,
  Mat< std::complex<T> >&             Q,
  Mat< std::complex<T> >&             Z,
  const Base< std::complex<T>, T1 >& A,
  const Base< std::complex<T>, T2 >& B,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
    
  const bool status = auxlib2::qz(AA, BB, Q, Z, A, B, 'b');
    
  if(status == false)
    {
    AA.reset();
    BB.reset();
    Q.reset();
    Z.reset();
    arma_bad("qz(): failed to converge", false);
    }
    
  return status;
  }


//! @}
