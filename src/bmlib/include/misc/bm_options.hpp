/*################################################################################
  ##
  ##   Copyright (C) 2011-2018 Keith O'Hara
  ##
  ##   This file is part of the BM++ C++ library.
  ##
  ##   BM++ is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   BM++ is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with BM++. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

/*
 * BM++ options
 */

#pragma once

#ifdef USE_RCPP_ARMADILLO
    #include <RcppArmadillo.h>
#else
    #ifndef ARMA_DONT_USE_WRAPPER
        #define ARMA_DONT_USE_WRAPPER
    #endif
    #include "armadillo"
#endif

#ifndef STATS_GO_INLINE
    #define STATS_GO_INLINE
#endif

#ifndef STATS_USE_ARMA
    #define STATS_USE_ARMA
#endif

#if defined(_OPENMP) && !defined(BM_NO_OPENMP) && !defined(BM_USE_OPENMP)
    #define BM_USE_OPENMP
#endif

#ifdef BM_USE_FFTW
    #include <fftw3.h>
#endif

#ifdef BM_USE_OPENMP
    #ifndef OPTIM_USE_OMP
        #define OPTIM_USE_OMP
    #endif
    
    #ifndef MCMC_USE_OMP
        #define MCMC_USE_OMP
    #endif

    #ifndef ARMA_USE_OPENMP
        #define ARMA_USE_OPENMP
    #endif
#endif

namespace bm {
    using uint_t = unsigned int;
    static const double inf = std::numeric_limits<double>::infinity();
}
