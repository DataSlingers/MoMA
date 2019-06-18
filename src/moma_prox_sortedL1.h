/*
 * Copyright 2013, M. Bogdan, E. van den Berg, W. Su, and E.J. Candes
 *
 * This file is part of SLOPE Toolbox version 1.0.
 *
 *   The SLOPE Toolbox is free software: you can redistribute it
 *   and/or  modify it under the terms of the GNU General Public License
 *   as published by the Free Software Foundation, either version 3 of
 *   the License, or (at your option) any later version.
 *
 *   The SLOPE Toolbox is distributed in the hope that it will
 *   be useful, but WITHOUT ANY WARRANTY; without even the implied
 *   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *   See the GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with the SLOPE Toolbox. If not, see
 *   <http://www.gnu.org/licenses/>.
 */

#ifndef MOMA_PROX_SORTEDL1
#define MOMA_PROX_SORTEDL1 1

#include "moma_base.h"

int evaluateProx(const arma::vec &y, const arma::vec &lambda, arma::vec &x, int n, const arma::uvec &order);
#endif // MOMA_PROX_SORTEDL1
