/*
 * Add a QZ Decomposition to the Armadillo C++ Matrix Algebra Library.
 *
 * It's not easy to extend classes so we'll just
 * define a new class 'auxlib2' for now.
 *
 * Keith O'Hara
 * 29/07/15
 */


//! \addtogroup auxlib
//! @{

class auxlib2
  {
  public:
  
  template<const uword row, const uword col>
  struct pos
    {
    static const uword n2 = row + col*2;
    static const uword n3 = row + col*3;
    static const uword n4 = row + col*4;
    };
  
  //
  // QZ decomposition
      
  template<typename T, typename T1, typename T2>
  inline static bool qz(Mat<T>& Smat, Mat<T>& Tmat, Mat<T>& vsl, Mat<T>& vsr, const Base<T,T1>& X, const Base<T,T2>& Y, const char side);
      
  template<typename T, typename T1, typename T2>
  inline static bool qz(Mat< std::complex<T> >& Smat, Mat< std::complex<T> >& Tmat, Mat< std::complex<T> >& vsl, Mat< std::complex<T> >& vsr, const Base< std::complex<T>, T1 >& X, const Base< std::complex<T>, T2 >& Y, const char side);

  };


//! @}
