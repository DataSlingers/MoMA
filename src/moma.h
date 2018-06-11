// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// Global #includes and #defines
#include "moma_base.h"

// Logging code
#include "moma_logging.h"

// Prox operators
#include "moma_prox.h"

// Prototypes
// moma_logging.cpp
void moma_set_logger_level_cpp(int);
int moma_get_logger_level_cpp();
void moma_log_cpp(int, Rcpp::StringVector);
