#ifndef MOMA_HEAP
#define MOMA_HEAP 1
#include "moma_base.h"
#include "moma_logging.h"

class HeapNode
{
  public:
    HeapNode(int i = -1, double l = -1.0) : id(i), lambda(l){};
    HeapNode &operator=(const HeapNode &source)
    {
        id     = source.id;
        lambda = source.lambda;
        return *this;
    }
    int id;         // value: id-th beta
    double lambda;  // key: for id-th group and its next groupto merge at lambda
    void print() { MoMALogger::debug("") << "lambda: " << lambda << "id: " << id; }
};

// comparision between heap nodes
bool gt(const HeapNode &left, const HeapNode &right);

class FusedGroups;
class Heap
{
  public:
    Heap(int n = 0);
    void heap_print();
    HeapNode heap_peek_min();
    bool is_empty();
    void heapify();
    std::vector<HeapNode> heap_storage;
    void remove(int id, FusedGroups *fg);

    int change_lambda_by_id(int id, double new_lambda, FusedGroups *fg);
    bool is_minheap();

  private:
    void swap(int i, int j, FusedGroups *fg);
    void siftup(int i, FusedGroups *fg);
    void siftdown(int current_node, FusedGroups *fg);

    int min_child(int i);

    // A constant, where the non-existing
    // child is assumed to be located.
    // Used only in function `min_child`
    const int NO_CHILD = -11;
};
#endif
