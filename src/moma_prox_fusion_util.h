#ifndef MOMA_PROX_FUSION_UTIL
#define MOMA_PROX_FUSION_UTIL 1
#include "moma_base.h"
#include "moma_heap.h"

class FusionGroups;
class Group{
    int head;
    int tail;
    int parent;
    double lambda;
    double beta;
    double slope;
    friend class FusionGroups;
public:
    Group(int t=-1,int b=-1,int p=-1,double lambda=-1,double beta = -1,double slope = 0):head(t),tail(b),parent(p),lambda(lambda),beta(beta),slope(slope){};
    void print(){
        Rcpp::Rcout<<"head: " << head 
        << "tail: " << tail 
        << "parent: " << parent
        << "lambda:\t" << lambda
        << "beta:\t" << beta
        << "slope: " << slope
        << "\n";
    }
};

class FusionGroups{
public:
    // The only heap manipulations needed
    friend int heap_change_lambda(std::vector<HeapNode> &, int id,double);
    friend void heap_delete(std::vector<HeapNode> &, int id);

    // Constructor
    FusionGroups(const arma::vec &x);

    // Manipulation on a group
    void print();
    bool is_valid(int this_node);
    int pre_group(int this_group);
    int next_group(int this_group);
    int group_size(int this_group);

    // Print beta
    arma::vec find_beta_at(double target_lam);

    // Calculation concerning lines
    double lines_meet_at(double x1,double x2,double k1,double k2,double y1,double y2);
    double line_value_at(double x,double y,double k,double x_);

    // Merge node dst with the group next to it
    void merge();
    double next_lambda();
    bool all_merged();
    // Some macro
    const int NO_PRE = -2;
    const int NO_NEXT = -3;
    const int INFTY = 2 << 17;

    // A vector stroing all the beta values
    std::vector<Group> g;
    std::vector<HeapNode> pq;
};
#endif
