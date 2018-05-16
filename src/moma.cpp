// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"


// [[Rcpp::export]]
int doubleMe(const int & x) {

    return x+x;
}
