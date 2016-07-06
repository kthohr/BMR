/*################################################################################
  ##
  ##   MacroLib by Keith O'Hara Copyright (C) 2011-2016
  ##   This file is part of the C++ MacroLib library.
  ##
  ##   MacroLib is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   MacroLib is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ################################################################################*/

/*
 * BVARM, auxiliary functions
 */

// geometric decay
double decay_geo(double x, double HP_4)
{
    double ret = std::pow(HP_4,-x + 1);
    return ret;
}

// harmonic decay
double decay_harm(double x, double HP_4)
{
    double ret = std::pow(x, HP_4);
    return ret;
}
