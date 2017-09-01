/*################################################################################
  ##
  ##   Copyright (C) 2011-2017 Keith O'Hara
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
  ################################################################################*/

#ifndef _bmlib_HPP
#define _bmlib_HPP

#include "misc/bm_options.hpp"

#include "optim/optim.hpp"
#include "mcmc/mcmc.hpp"
#include "stats/stats.hpp"

namespace bm
{
    #include "misc/bind.hpp"
    #include "misc/byrow.hpp"
    #include "misc/cube_to_mat.hpp"
    #include "misc/embed.hpp"
    #include "misc/inside_trans.hpp"
    #include "misc/log_GPR.hpp"
    #include "misc/lyapunov_dbl.hpp"
    #include "misc/uvec_linspace.hpp"
    #include "misc/zero_rows.hpp"

    // VARs
    #include "var/var_sim.hpp"

    #include "var/bvarm.hpp"
    #include "var/bvars.hpp"
    #include "var/bvarw.hpp"
    #include "var/cvar.hpp"

    #include "var/bvartvp.hpp"

    // filters
    #include "filter/chandrasekhar.hpp"
    #include "filter/kalman.hpp"

    #include "filter/dk_filter.hpp"

    // LRES
    #include "lres/gensys_solver.hpp"
    #include "lres/uhlig_solver.hpp"

    // DSGE
    #include "dsge/dsge_simulate.hpp"

    #include "dsge/gensys.hpp"
    #include "dsge/uhlig.hpp"

    #include "dsge/dsge_class.hpp"

    // DSGE-VAR

    #include "dsgevar/dsgevar_class.hpp"
}

#endif
