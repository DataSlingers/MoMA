
// a non-existing child of the node of a heap is far away
#include "moma_prox_fusion_util.h"
#include "moma_heap.h"
int sgn(double val)
{
    return (double(0) < val) - (val < double(0));
}

// Constructor
FusedGroups::FusedGroups(const arma::vec &x) : heap(x.n_elem - 1)
{
    int n = x.n_elem;
    if (n <= 1)
    {
        MoMALogger::error("TODO: deal with scalar");
    }
    g.resize(n);

    // Initialize `head`, `tail`, `parent`, `lambda` and `beta`
    for (int i = 0; i < g.size(); i++)
    {
        g[i] = Group(i, i, i, 0, x(i));
    }
    // slope
    g[0].slope = -sgn(double(x(0) - x(1)));
    for (int i = 1; i < g.size() - 1; i++)
    {
        double s   = -(sgn(x(i) - x(i - 1)) + sgn(x(i) - x(i + 1)));
        g[i].slope = s;
    }
    g[n - 1].slope = -(sgn(x(n - 1) - x(n - 2)));

    // Heap lambda;
    for (int i = 0; i < heap.heap_storage.size(); i++)
    {
        // next merge point of group i and i+1
        double h = 0;
        h        = lines_meet_at(0, 0, g[i + 1].slope, g[i].slope, g[i + 1].beta, g[i].beta);
        heap.heap_storage[i] = HeapNode(i, h);
    }
    heap.heapify();
    for (int i = 0; i < heap.heap_storage.size(); i++)
    {
        int index_in_g            = heap.heap_storage[i].id;
        g[index_in_g].map_to_heap = i;
    }
    g[n - 1].map_to_heap = NOT_IN_HEAP;
    return;
}

void FusedGroups::print()
{
    MoMALogger::debug("") << "Grouping now is";
    for (int i = 0; i < g.size(); i++)
    {
        if (is_valid(i))
        {
            g[i].print();
            if (g[i].map_to_heap != NOT_IN_HEAP && g[i].map_to_heap >= heap.heap_storage.size())
            {
                MoMALogger::error("Exceeds heap limit")
                    << g[i].map_to_heap << "while heap size is " << heap.heap_storage.size();
            }
        }
        else
        {
            MoMALogger::debug("") << "=====";
            g[i].print();
        }
    }
    MoMALogger::debug("");
}

bool FusedGroups::is_valid(int this_node)
{
    return g[this_node].parent == this_node;
}

int FusedGroups::pre_group(int this_group)
{
    if (!is_valid(this_group))
    {
        MoMALogger::error("Only valid groups can be accessed");
    }
    if (this_group == 0)
    {
        return NO_PRE;
    }
    return g[this_group - 1].parent;
}

int FusedGroups::next_group(int this_group)
{
    if (!is_valid(this_group))
    {
        MoMALogger::error("Only valid groups be accessed");
    }
    if (g[this_group].tail == g.size() - 1)
    {
        return NO_NEXT;
    }
    else
    {
        return g[this_group].tail + 1;
    }
}

int FusedGroups::group_size(int this_group)
{
    if (!is_valid(this_group))
    {
        MoMALogger::error("Only valid groups be accessed");
    }
    return g[this_group].tail - g[this_group].head + 1;
}

// line_value_at evaluates the y-value of a line,
// who has slope is k and goes through point (x,y),
// at x_.
double FusedGroups::line_value_at(double x, double y, double slope, double x_)
{
    return y + slope * (x_ - x);
}

arma::vec FusedGroups::find_beta_at(double target_lam)
{
    int n       = (this->g).size();
    arma::vec x = arma::zeros<arma::vec>(n);
    for (int i = 0; i != NO_NEXT;)
    {
        double betaj = line_value_at(g[i].lambda, g[i].beta, g[i].slope, target_lam);
        for (int j = g[i].head; j <= g[i].tail; j++)
        {
            x(j) = betaj;
        }
        i = next_group(i);
    }
    return x;
}

// Find the x value of the intersection of two lines.
double FusedGroups::lines_meet_at(double x1,
                                  double x2,
                                  double slope1,
                                  double slope2,
                                  double y1,
                                  double y2)
{
    if (std::abs(slope1 - slope2) < 1e-10)
    {
        // Note abs(slope1 - slope2) < 1e-10
        // does not work on Linux.
        return MOMA_INFTY;
    }
    return ((y1 - y2) - (slope1 * x1 - slope2 * x2)) / (-slope1 + slope2);
}

void FusedGroups::merge()
{
    HeapNode node = heap.heap_peek_min();
    // Node `dst` will absorb the info of `src` node, `src` will be then marked
    // invalid
    int dst           = node.id;
    double new_lambda = node.lambda;
    int src           = this->next_group(dst);

    if (!is_valid(dst) || src == NO_NEXT)
    {
        MoMALogger::error("Only valid groups can be merged: merge point is not valid");
    }
    if (dst >= src)
    {
        MoMALogger::error("dst_grp should be in front of src_grp");
    }

    // update beta
    g[dst].beta = g[dst].beta + g[dst].slope * (new_lambda - g[dst].lambda);

    // update lambda
    g[dst].lambda = new_lambda;

    // update slope
    int pre_group  = this->pre_group(dst);
    int next_group = this->next_group(src);
    int sgn1       = 0;
    int sgn2       = 0;
    if (next_group != NO_NEXT)
    {
        sgn2 = sgn(g[dst].beta - g[next_group].beta);
    }
    if (pre_group != NO_PRE)
    {
        sgn1 = sgn(g[dst].beta - g[pre_group].beta);
    }
    g[dst].slope = -1 / double(this->group_size(dst) + this->group_size(src)) * (sgn1 + sgn2);

    // set up pointers
    int last_node       = g[src].tail;
    g[dst].tail         = last_node;
    g[src].parent       = dst;
    g[last_node].parent = dst;

    // update heap
    if (pre_group != NO_PRE)
    {
        double lambda_pre = lines_meet_at(g[pre_group].lambda, g[dst].lambda, g[pre_group].slope,
                                          g[dst].slope, g[pre_group].beta, g[dst].beta);
        heap.change_lambda_by_id(g[pre_group].map_to_heap, lambda_pre, this);
    }
    if (next_group != NO_NEXT)
    {
        double lambda_next = lines_meet_at(g[next_group].lambda, g[dst].lambda, g[next_group].slope,
                                           g[dst].slope, g[next_group].beta, g[dst].beta);
        heap.change_lambda_by_id(g[dst].map_to_heap, lambda_next, this);
        heap.remove(g[src].map_to_heap, this);
    }
    else
    {
        heap.remove(g[dst].map_to_heap, this);
    }
}

double FusedGroups::next_lambda()
{
    HeapNode n = heap.heap_peek_min();
    return n.lambda;
};

bool FusedGroups::all_merged()
{
    return heap.is_empty();
};
