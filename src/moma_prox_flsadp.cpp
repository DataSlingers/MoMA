
#include "moma_prox_flsadp.h"
#include "moma_logging.h"

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

// This code is an adapted version of the package cited above.
// Better modularity and readability. Commments are added
// in order to explain the algorithm. Avoid explicit
// memory management.

// Reference:
// Johnson, Nicholas A.
// "A dynamic programming algorithm for the fused lasso and l 0-segmentation."
// Journal of Computational and Graphical Statistics 22.2 (2013): 246-260.

double Msg::Argmax(double *max_val)
{
    // printf("Enter MaxMsg");
    const std::vector<MsgElt> &buf = buf_;

    double lin_left  = 0.0;
    double quad_left = 0.0;

    int last_idx = start_idx_ + len_ - 1;

    // Step from left to right
    for (int knot_idx = start_idx_, m = len_; m; --m, ++knot_idx)
    {
        bool end_knot = (knot_idx == last_idx);
        const MsgElt &k =
            (knot_idx == start_idx_) ? init_knot_ : (end_knot ? end_knot_ : buf[knot_idx]);

        double x1 = (knot_idx == last_idx - 1) ? end_knot_.x_ : buf[knot_idx + 1].x_;

        if (k.sgn_)
        {
            lin_left += k.lin_;
            quad_left += k.quad_;
        }
        else
        {
            lin_left -= k.lin_;
            quad_left -= k.quad_;
        }

        if (quad_left == 0.0)
        {
            continue;
        }

        double hit_x = -lin_left / (2.0 * quad_left);
        if (hit_x < x1)
        {
            if (max_val)
            {
                *max_val = hit_x * lin_left + hit_x * hit_x * quad_left;
            }
            return (hit_x);
        }
    }
    MoMALogger::error("MaxMsg : failed to maximize message");
    return (-1);
}

void Msg::InitMsg(int n, int init_sz, double lin, double quad, double lambda2)
{
    len_       = 2;
    start_idx_ = init_sz / 2;
    buf_       = std::vector<MsgElt>(init_sz);
    back_pointers.resize(n * 2);

    int i = start_idx_;

    buf_[i].x_    = R_NegInf;
    buf_[i].sgn_  = true;
    buf_[i].lin_  = lin;
    buf_[i].quad_ = quad;

    buf_[i + 1].x_    = R_PosInf;
    buf_[i + 1].sgn_  = false;
    buf_[i + 1].lin_  = lin;
    buf_[i + 1].quad_ = quad;

    init_knot_.x_    = R_NegInf;
    init_knot_.sgn_  = true;
    init_knot_.lin_  = lambda2;
    init_knot_.quad_ = 0.0;

    end_knot_.x_    = R_PosInf;
    end_knot_.sgn_  = false;
    end_knot_.lin_  = -lambda2;
    end_knot_.quad_ = 0.0;
    MoMALogger::debug("Finish initialization.");
}

void Msg::ShiftMsg(int check_freq)
{
    if (len_ > buf_.size() - 20 * check_freq)
    {
        std::vector<MsgElt> new_buf(buf_.size() * 3);
        std::vector<MsgElt> *old_buf = &buf_;

        int new_start = (new_buf.size() / 4);
        int old_start = start_idx_;

        for (int k = 0; k < len_; ++k)
        {
            new_buf[k + new_start] = (*old_buf)[k + old_start];
        }

        buf_.swap(new_buf);
        start_idx_ = new_start;
    }

    if (start_idx_ < 5 * check_freq)
    {
        int new_start            = ((buf_.size() - len_) / 2);
        int old_start            = start_idx_;
        std::vector<MsgElt> *buf = &buf_;

        for (int k = len_ - 1; k >= 0; --k)
        {
            (*buf)[k + new_start] = (*buf)[k + old_start];
        }
        start_idx_ = new_start;
    }
    else if (start_idx_ + len_ > buf_.size() - 5 * check_freq)
    {
        int new_start            = ((buf_.size() - len_) / 2);
        int old_start            = start_idx_;
        std::vector<MsgElt> *buf = &buf_;

        for (int k = 0; k < len_; ++k)
        {
            (*buf)[k + new_start] = (*buf)[k + old_start];
        }
        start_idx_ = new_start;
    }
}

void Msg::UpdMsg(double lambda2, double lin, double quad, int bp_idx)
{
    std::vector<MsgElt> &buf = buf_;

    double lin_left    = lin;
    double quad_left   = quad;
    int new_knot_start = -3;

    // Algorithm 2, line 4-14
    for (int knot_idx = start_idx_, m = len_; m; --m, ++knot_idx)
    {
        const MsgElt &k = buf[knot_idx];
        double x1       = buf[knot_idx + 1].x_;
        if (k.sgn_)
        {
            lin_left += k.lin_;
            quad_left += k.quad_;
        }
        else
        {
            lin_left -= k.lin_;
            quad_left -= k.quad_;
        }

        double hit_x = (lambda2 - lin_left) / (2.0 * quad_left);
        if (hit_x < x1)
        {
            // place a knot here
            new_knot_start        = knot_idx - 1;
            back_pointers(bp_idx) = hit_x;
            break;
        }
    }

    double lin_right  = lin;
    double quad_right = quad;
    int new_knot_end  = -2;
    int end_idx       = start_idx_ + len_ - 1;

    double neg_lam2 = -lambda2;

    // Algorihtm 2, line 15-25
    for (int knot_idx = end_idx, m = len_; m; --m, --knot_idx)
    {
        const MsgElt &k = buf[knot_idx];
        double x1       = buf[knot_idx - 1].x_;

        if (k.sgn_)
        {
            lin_right -= k.lin_;
            quad_right -= k.quad_;
        }
        else
        {
            lin_right += k.lin_;
            quad_right += k.quad_;
        }

        double hit_x = (neg_lam2 - lin_right) / (2.0 * quad_right);
        if (hit_x > x1)
        {
            // place a knot here
            new_knot_end              = knot_idx + 1;
            back_pointers(bp_idx + 1) = hit_x;
            break;
        }
    }

    // Prepend and append, Algorithm 2 line 26-31
    MsgElt *k0 = &buf[new_knot_start];
    k0->x_     = R_NegInf;
    k0->sgn_   = true;
    k0->lin_   = lambda2;
    k0->quad_  = 0.0;
    init_knot_ = *k0;

    k0        = &buf[new_knot_start + 1];
    k0->x_    = back_pointers(bp_idx);
    k0->sgn_  = true;
    k0->lin_  = lin_left - lambda2;
    k0->quad_ = quad_left;

    k0        = &buf[new_knot_end];
    k0->x_    = R_PosInf;
    k0->sgn_  = false;
    k0->lin_  = -lambda2;
    k0->quad_ = 0.0;
    end_knot_ = *k0;

    k0        = &buf[new_knot_end - 1];
    k0->x_    = back_pointers(bp_idx + 1);
    k0->sgn_  = false;
    k0->lin_  = lin_right + lambda2;
    k0->quad_ = quad_right;

    start_idx_ = new_knot_start;
    len_       = 1 + new_knot_end - new_knot_start;
}

void Msg::UpdMsgOpt(double lambda2, double lin, double quad, int bp_idx)
{
    // An optimized version

    std::vector<MsgElt> &buf = buf_;

    buf[start_idx_].x_            = R_NegInf;
    buf[start_idx_ + len_ - 1].x_ = R_PosInf;

    double lin_left    = lin;
    double quad_left   = quad;
    int new_knot_start = -3;

    int knot_idx    = start_idx_;
    const MsgElt *k = &(init_knot_);
    double x1       = buf[knot_idx + 1].x_;
    if (k->sgn_)
    {
        lin_left += k->lin_;
        quad_left += k->quad_;
    }
    else
    {
        lin_left -= k->lin_;
        quad_left -= k->quad_;
    }

    double hit_x = (lambda2 - lin_left) / (2.0 * quad_left);
    if (hit_x < x1)
    {
        // place a knot here
        new_knot_start        = knot_idx - 1;
        back_pointers(bp_idx) = hit_x;
    }
    else
    {
        ++knot_idx;
        k = &buf[knot_idx];
        for (;; ++knot_idx, ++k)
        {
            x1 = k[1].x_;
            if (k->sgn_)
            {
                lin_left += k->lin_;
                quad_left += k->quad_;
            }
            else
            {
                lin_left -= k->lin_;
                quad_left -= k->quad_;
            }

            hit_x = (lambda2 - lin_left) / (2.0 * quad_left);
            if (hit_x < x1)
            {
                // place a knot here
                new_knot_start = knot_idx - 1;

                back_pointers(bp_idx) = hit_x;
                break;
            }
        }
    }

    double lin_right  = lin;
    double quad_right = quad;
    int new_knot_end  = -2;
    int end_idx       = start_idx_ + len_ - 1;

    double neg_lam2 = -lambda2;

    knot_idx = end_idx;
    k        = &(end_knot_);
    x1       = buf[knot_idx - 1].x_;

    if (k->sgn_)
    {
        lin_right -= k->lin_;
        quad_right -= k->quad_;
    }
    else
    {
        lin_right += k->lin_;
        quad_right += k->quad_;
    }

    hit_x = (neg_lam2 - lin_right) / (2.0 * quad_right);
    if (hit_x > x1)
    {
        // place a knot here
        new_knot_end              = knot_idx + 1;
        back_pointers(bp_idx + 1) = hit_x;
    }
    else
    {
        --knot_idx;
        k = &buf[knot_idx];
        for (;; --knot_idx, --k)
        {
            x1 = k[-1].x_;

            if (k->sgn_)
            {
                lin_right -= k->lin_;
                quad_right -= k->quad_;
            }
            else
            {
                lin_right += k->lin_;
                quad_right += k->quad_;
            }

            hit_x = (neg_lam2 - lin_right) / (2.0 * quad_right);
            if (hit_x > x1)
            {
                // place a knot here
                new_knot_end              = knot_idx + 1;
                back_pointers(bp_idx + 1) = hit_x;
                break;
            }
        }
    }

    MsgElt *k0 = &buf[new_knot_start + 1];
    k0->x_     = back_pointers(bp_idx);
    k0->sgn_   = true;
    k0->lin_   = lin_left - lambda2;
    k0->quad_  = quad_left;

    k0        = &buf[new_knot_end - 1];
    k0->x_    = back_pointers(bp_idx + 1);
    k0->sgn_  = false;
    k0->lin_  = lin_right + lambda2;
    k0->quad_ = quad_right;

    start_idx_ = new_knot_start;
    len_       = 1 + new_knot_end - new_knot_start;
}

arma::vec Msg::BackTrace(int seq_len, double last_msg_max)
{
    arma::vec x_hat(seq_len);
    double z = x_hat(seq_len - 1) = last_msg_max;
    int i                         = seq_len - 2;

    int bp_idx = (2 * (seq_len - 2));
    for (int idx = seq_len - 1; idx; --idx, bp_idx -= 2, --i)
    {
        if (z < back_pointers(bp_idx))
        {
            z = x_hat[i] = back_pointers(bp_idx);
        }
        else if (z > back_pointers(bp_idx + 1))
        {
            z = x_hat[i] = back_pointers(bp_idx + 1);
        }
        else
        {
            x_hat[i] = z;
        }
    }
    return x_hat;
}

arma::vec myflsadp(const arma::vec &x, double lambda2, int init_buf_sz)
{
    // lambda2 is the penalty level on the difference
    // of adjacent elements.
    int seq_len = x.n_elem;

    if (lambda2 == 0.0)
    {
        return x;
    }

    if (seq_len < 2)
    {
        MoMALogger::error("Fused lasso: input vector has length less than two");
    }

    int check_freq = 40;
    if (init_buf_sz < 30 * check_freq)
    {
        MoMALogger::error("Fused lasso (DP): initial buffer size is too small");
    }

    Msg msg;
    msg.InitMsg(seq_len, init_buf_sz, 0.0, 0.0, lambda2);

    int check_msg = check_freq - 1;

    msg.UpdMsg(lambda2, x(0), -0.5, 0);
    MoMALogger::debug("After first UpdMsg");

    for (int j = 1, bp = 2; j < seq_len; ++j, bp += 2, --check_msg)
    {
        msg.UpdMsgOpt(lambda2, x(j), -0.5, bp);
        // call ShiftMsg periodically
        if (!check_msg)
        {
            check_msg = check_freq - 1;
            msg.ShiftMsg(check_freq);
        }
    }

    double last_msg_max = msg.Argmax(NULL);
    return msg.BackTrace(seq_len, last_msg_max);
}
