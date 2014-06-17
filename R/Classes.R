# 16/06/2014
BVARM <- function(x, ...) UseMethod("BVARM")
BVARS <- function(x, ...) UseMethod("BVARS")
BVARW <- function(x, ...) UseMethod("BVARW")
CVAR <- function(x, ...) UseMethod("CVAR")
BVARTVP <- function(x, ...) UseMethod("BVARTVP")
#
forecast <- function(x, ...) UseMethod("forecast")
IRF <- function(x, ...) UseMethod("IRF")
#
SDSGE <- function(x, ...) UseMethod("SDSGE")
DSGESim <- function(x, ...) UseMethod("DSGESim")
EDSGE <- function(x, ...) UseMethod("EDSGE")
DSGEVAR <- function(x, ...) UseMethod("DSGEVAR")
modecheck <- function(x, ...) UseMethod("modecheck")