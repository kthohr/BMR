################################################################################
##
##   Copyright (C) 2011-2018 Keith O'Hara
##
##   This file is part of the R package BMR.
##
##   The R package BMR is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 2 of the License, or
##   (at your option) any later version.
##
##   The R package BMR is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
################################################################################

FEVD <- function(obj, ...) UseMethod("FEVD")
forecast <- function(obj, ...) UseMethod("forecast")
IRF <- function(obj, ...) UseMethod("IRF")

mode_check <- function(obj, ...) UseMethod("mode_check")
states <- function(obj, ...) UseMethod("states")
