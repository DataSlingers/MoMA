#ifndef MOMA_PROX_FUSION_UTIL
#define MOMA_PROX_FUSION_UTIL 1
#include "moma_base.h"
#include "moma_heap.h"
class FusedGroups;
class Group
{
  public:
    // range of the group, note they are continuous
    int head;
    int tail;
    // When two groups (say A and B) merge,
    // `parent` of the last node of group B will point to A
    // Note all nodes are initialized with `parent` pointing to itself
    int parent;
    // The following infomation is valid only when `parent` points to itself
    double lambda;
    double beta;
    double slope;
    int map_to_heap;
    friend class FusedGroups;
    Group(int h         = -1,
          int t         = -1,
          int p         = -1,
          double lambda = -1,
          double beta   = -1,
          double slope  = 0)
        : head(h), tail(t), parent(p), lambda(lambda), beta(beta), slope(slope){};
    void print()
    {
        MoMALogger::debug("") << "[" << head << "," << tail << "] map_to_heap: " << map_to_heap
                              << "(lambda:" << lambda << ",beta:" << beta << ",slope: " << slope
                              << ")";
    }
};

class FusedGroups
{
  public:
    // Constructor
    FusedGroups(const arma::vec &x);
    // Merge the next two nodes.
    // Note if multiple pairs of nodes is to be merged at the same lambda, only
    // one pair will be merged
    void merge();
    // Return the next lambda at which merge happens
    double next_lambda();
    // Check if all beta's are merged
    bool all_merged();
    // Evaluate beta by extending the lines
    arma::vec find_beta_at(double target_lam);

    // Manipulation on a group
    void print();
    bool is_valid(int this_node);
    int pre_group(int this_group);
    int next_group(int this_group);
    int group_size(int this_group);

    // Calculation concerning lines
    double line_value_at(double x, double y, double slope, double x_);
    double lines_meet_at(double x1, double x2, double k1, double k2, double y1, double y2);

    // Some constants
    // Used when the group includes beta_1
    const int NO_PRE = -2;
    // Used when the group includes beta_p
    const int NO_NEXT            = -3;
    static const int NOT_IN_HEAP = -4;
    // This constant is used in the function
    // lines_meet_at, where we might bump
    // into situation that the two lines are
    // parallel. In order to deal with this,
    // we return MOMA_INFTY.

    // A vector stroing all the beta values
    std::vector<Group> g;
    Heap heap;
};
#endif
