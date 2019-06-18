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

#include "moma_prox_sortedL1.h"
#include "moma_logging.h"

// This is a slight modified version of the code provided by
// M. Bogdan, E. van den Berg, W. Su, and E.J. Candes
// http://statweb.stanford.edu/~candes/SortedL1/

int evaluateProx(const arma::vec &y,
                 const arma::vec &lambda,
                 arma::vec &x,
                 int n,
                 const arma::uvec &order)
{
    double d;

    arma::vec s(n);
    arma::vec w(n);
    arma::uvec idx_i(n);
    arma::uvec idx_j(n);

    int i, j, k;

    k = 0;
    for (i = 0; i < n; i++)
    {
        idx_i(k) = i;
        idx_j(k) = i;
        s(k)     = y(i) - lambda(i);
        w(k)     = s(k);

        while ((k > 0) && (w[k - 1] <= w(k)))
        {
            k--;
            idx_j(k) = i;
            s(k) += s[k + 1];
            w(k) = s(k) / (i - idx_i(k) + 1);
        }

        k++;
    }

    for (j = 0; j < k; j++)
    {
        d = w(j);
        if (d < 0)
        {
            d = 0;
        }
        for (i = idx_i(j); i <= idx_j(j); i++)
        {
            x[order(i)] = d;
        }
    }

    return 0;
}
