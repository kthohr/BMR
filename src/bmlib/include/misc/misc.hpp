/*################################################################################
  ##
  ##   Copyright (C) 2011-2018 Keith O'Hara
  ##
  ##   This file is part of the BMLib C++ library.
  ##
  ##   BMLib is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   BMLib is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with BMLib. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

#ifndef _bmlib_misc_HPP
#define _bmlib_misc_HPP

#include "bm_options.hpp"

#include "optim/optim.hpp"
#include "mcmc/mcmc.hpp"
#include "stats/stats.hpp"

namespace bm
{
    #include "bind.hpp"
    #include "byrow.hpp"
    #include "cube_to_mat.hpp"
    #include "embed.hpp"
    #include "inside_trans.hpp"
    #include "log_GPR.hpp"
    #include "lyapunov_dbl.hpp"
    #include "uvec_linspace.hpp"
    #include "zero_rows.hpp"

    #include "tictoc.hpp"
}

#endif
