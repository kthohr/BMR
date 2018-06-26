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
 * simple tictoc functionality
 */

using comptime = std::chrono::time_point<std::chrono::system_clock>;

inline
comptime
tic()
{
    return std::chrono::system_clock::now();
}

inline
void
tictoc(comptime time_inp)
{
    comptime time_now = std::chrono::system_clock::now();

    std::chrono::duration<double> run_time = time_now - time_inp;

    //

    std::time_t end_time = std::chrono::system_clock::to_time_t(time_now);
        
    std::cout << "finished computation at " << std::ctime(&end_time)
              << "total runtime: " << run_time.count() << "s\n";
}
