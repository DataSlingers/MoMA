/* MoMA Base Header
 *
 * This is where we put all #defines and #includes which will be
 * seen by all other headers
 */

#ifndef MOMA_BASE_H
#define MOMA_BASE_H 1

#include <sstream>
#include <iostream>

// We only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// For difficult smoothing matrices, we may encounter artificially small eigenvalues:
// we add a small "nugget" here to regularize the computations
#define MOMA_EIGENVALUE_REGULARIZATION 0.01
static constexpr double MOMA_INFTY = std::numeric_limits<double>::infinity();
static const arma::vec MOMA_EMPTY_GRID_OF_LENGTH1 = -arma::ones<arma::vec>(1);
#define MOMA_FUSEDLASSODP_BUFFERSIZE 5000
#endif
