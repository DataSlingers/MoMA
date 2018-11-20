#ifndef MOMA_PROX_FLSADP
#define MOMA_PROX_FLSADP 1

#include "moma_base.h"

struct MsgElt {
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

class Msg {
public:
    std::vector<MsgElt> buf_;
    arma::vec back_pointers;
    int start_idx_;
    int len_;

    MsgElt init_knot_;
    MsgElt end_knot_;

    void InitMsg(int n, int init_sz, double lin, double quad, double lambda2);
    void UpdMsg   (double lambda2, double lin, double quad, int bp_idx);
    void UpdMsgOpt(double lambda2, double lin, double quad, int bp_idx);
    double Argmax(double * max_val);

    // This data structure supports prepend and append;
    // To implement this, first allocate a long array, 
    // then start filling data in the middle and extend to both ends. Shift 
    // the space filled with data towards the center of the array periodically
    void ShiftMsg(int check_freq);
    arma::vec BackTrace(int seq_len, double last_msg_max);
    void print();
};

arma::vec myflsadp(const arma::vec& x, double lambda2, int init_buf_sz = 5000);

#endif  // MOMA_PROX_FLSADP
