#ifndef MOMA_HEAP
#define MOMA_HEAP 1
#include "moma_base.h"
#include "moma_logging.h"

// The non-existing child is located at infinity
#define NO_CHILED 2 << 19
class HeapNode{
public:
    HeapNode(int i=-1,double l=-1.0):id(i),lambda(l){};
    HeapNode& operator = (const HeapNode &source){
        id = source.id;
        lambda = source.lambda;
        return *this;
    } 
    int id;     // value: id-th beta
    double lambda;  // key: for id-th group and its next groupto merge at lambda
    void print(){
        Rcpp::Rcout<<"lambda: " << lambda << "id: " << id << "\n";
    }
};
bool gt(const HeapNode &left, const HeapNode &right);
int min_child(std::vector<HeapNode> &h, int i);
void swap(std::vector<HeapNode> &h,int i, int j);
void siftup(std::vector<HeapNode> &heap, int i);
void siftdown(std::vector<HeapNode> &h, int current_node);
bool is_minheap(std::vector<HeapNode> &heap);
bool is_empty(std::vector<HeapNode> &heap);
HeapNode heap_peek_min(std::vector<HeapNode> &heap);
void heap_print(const std::vector<HeapNode> &q);

void heap_delete(std::vector<HeapNode> &heap, int id);
int heap_change_lambda(std::vector<HeapNode> &heap, int id, double new_lambda);
#endif
