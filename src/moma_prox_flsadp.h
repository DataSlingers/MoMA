#ifndef MOMA_PROX_FLSADP
#define MOMA_PROX_FLSADP 1

#include "moma_base.h"

// Copyright (c) 2012, Nicholas A. Johnson
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. All advertising materials mentioning features or use of this software
//    must display the following acknowledgement:
//    This product includes software developed by Nicholas A. Johnson.
// 4. Neither the name of Nicholas A. Johnson nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY Nicholas A. Johnson ''AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL Nicholas A. Johnson BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// This file is an adapted version of the package cited above.
// Better modularity and readability. Commments are added
// in order to explain the algorithm. Avoid explicit
// memory management.

// Reference:
// Johnson, Nicholas A.
// "A dynamic programming algorithm for the fused lasso and l 0-segmentation."
// Journal of Computational and Graphical Statistics 22.2 (2013): 246-260.

struct MsgElt
{
    // the location of the knot
    double x_;
    // the sign variable that tells us
    // whether this was a left or right end-point of the
    // segment
    bool sgn_;
    // a delta which can be used to reconstruct the function
    // if we move from the first knot to the last or from
    // the last to the first
    double lin_;
    double quad_;
};

class Msg
{
  public:
    std::vector<MsgElt> buf_;
    arma::vec back_pointers;
    int start_idx_;
    int len_;

    MsgElt init_knot_;
    MsgElt end_knot_;

    void InitMsg(int n, int init_sz, double lin, double quad, double lambda2);
    void UpdMsg(double lambda2, double lin, double quad, int bp_idx);
    void UpdMsgOpt(double lambda2, double lin, double quad, int bp_idx);
    double Argmax(double *max_val);

    // This data structure supports prepend and append;
    // To implement this, first allocate a long array,
    // then start filling data in the middle and extend to both ends. Shift
    // the space filled with data towards the center of the array periodically
    void ShiftMsg(int check_freq);
    arma::vec BackTrace(int seq_len, double last_msg_max);
};

arma::vec myflsadp(const arma::vec &x,
                   double lambda2,
                   int init_buf_sz = MOMA_FUSEDLASSODP_BUFFERSIZE);

#endif  // MOMA_PROX_FLSADP
