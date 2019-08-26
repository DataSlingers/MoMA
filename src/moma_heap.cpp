#include "moma_heap.h"
#include "moma_prox_fusion_util.h"

bool operator>(const HeapNode &left, const HeapNode &right)
{
    return left.lambda > right.lambda;
}

bool gt(const HeapNode &left, const HeapNode &right)
{
    return left > right;
}

Heap::Heap(int n)
{
    heap_storage.resize(n);
};

void Heap::heapify()
{
    std::make_heap(heap_storage.begin(), heap_storage.end(), gt);
}

// Find the smaller child of an element in a heap. Used in siftdown
int Heap::min_child(int i)
{
    int cur_size = heap_storage.size();
    int child    = i * 2 + 1;
    if (child >= cur_size)
    {
        // no children
        return NO_CHILD;
    }
    else if (child + 1 >= cur_size || !(heap_storage[child] > heap_storage[child + 1]))
    {
        // only child or first child is biggest child
        return child;
    }
    else
    {
        // second child exists and is smallest child
        return child + 1;
    }
}

// TODO: extra copy can be avoided in siftdown
void Heap::swap(int i, int j, FusedGroups *fg)
{
    // // DEBUG INFO
    MoMALogger::debug("Swapping ") << heap_storage[i].lambda << "and " << heap_storage[j].lambda;
    (*fg).g[heap_storage[i].id].map_to_heap = j;
    (*fg).g[heap_storage[j].id].map_to_heap = i;
    HeapNode tmp                            = heap_storage[i];
    heap_storage[i]                         = heap_storage[j];
    heap_storage[j]                         = tmp;
}

// In a min-heap, if the key (lambda in our case) decreases, sift it up
void Heap::siftup(int i, FusedGroups *fg)
{
    int parent = (i - 1) / 2;
    while (i != 0 && (heap_storage[parent] > heap_storage[i]))
    {
        Heap::swap(parent, i, fg);
        i      = parent;
        parent = (i - 1) / 2;
    }
}

// In a min-heap, if the key (lambda in our case) increases, sift it down
void Heap::siftdown(int current_node, FusedGroups *fg)
{
    int child = min_child(current_node);
    while (child != NO_CHILD && (heap_storage[current_node] > heap_storage[child]))
    {
        Heap::swap(child, current_node, fg);
        current_node = child;
        child        = min_child(child);
    }
}

// Change the key of any nodes;
int Heap::change_lambda_by_id(int i, double new_lambda, FusedGroups *fg)
{
    if (i < 0 || i >= heap_storage.size())
    {
        MoMALogger::error("Try to change lambda: no such id in current heap: ") << i;
    }
    double old_lambda      = heap_storage[i].lambda;
    heap_storage[i].lambda = new_lambda;
    if (old_lambda < new_lambda)
    {
        // // DEBUG INFO
        // MoMALogger::debug("(") << old_lambda << "," << heap[i].id << ")" << "->"
        // << new_lambda << " siftdown";
        siftdown(i, fg);
    }
    else
    {
        // //  DEBUG INFO
        // MoMALogger::debug("") << old_lambda << "," << heap[i].id << ")" << "->"
        // << new_lambda << " siftup";
        siftup(i, fg);
    }
    return i;
}

//  To delete an element, move it to the tail, pop it out, and then sift down
//  the node that replaces it
void Heap::remove(int i, FusedGroups *fg)
{
    if (i < 0 || i >= heap_storage.size())
    {
        MoMALogger::error("Try to delete: no such id in current heap: ") << i;
    }
    double old_lambda = heap_storage[i].lambda;
    Heap::swap(i, heap_storage.size() - 1, fg);
    (*fg).g[heap_storage[heap_storage.size() - 1].id].map_to_heap = FusedGroups::NOT_IN_HEAP;
    heap_storage.pop_back();
    if (old_lambda < heap_storage[i].lambda)
    {
        siftdown(i, fg);
    }
    else
    {
        siftup(i, fg);
    }
    return;
}

// Check if an array is a min heap
bool Heap::is_minheap()
{
    int i = 0;
    while (2 * i + 1 < heap_storage.size())
    {
        if (heap_storage[i] > heap_storage[2 * i + 1])
        {
            MoMALogger::warning("Not a min-heap")
                << heap_storage[i].lambda << "and" << heap_storage[2 * i + 1].lambda;
            return 0;
        }
        if (2 * i + 2 < heap_storage.size())
        {
            if (heap_storage[i] > heap_storage[2 * i + 2])
            {
                MoMALogger::warning("Not a min-heap")
                    << heap_storage[i].lambda << "and" << heap_storage[2 * i + 2].lambda;
                return 0;
            }
        }
        i++;
    }
    return 1;
}

bool Heap::is_empty()
{
    return heap_storage.size() == 0;
}

// Get the currently minimun value without deleting the node
HeapNode Heap::heap_peek_min()
{
    if (is_empty())
    {
        MoMALogger::error("You are peaking at an empty heap!");
    }
    HeapNode n = heap_storage.front();
    return n;
}

// Print the heap
void Heap::heap_print()
{
    MoMALogger::debug("") << "(lambda, id)\n";
    int cnt  = 0;
    int thre = 1;
    for (auto i : heap_storage)
    {
        Rcpp::Rcout << i.lambda << ", " << i.id + 1 << "\t";
        cnt++;
        if (cnt == thre)
        {
            Rcpp::Rcout << "\n";
            thre *= 2;
            cnt = 0;
        }
    }
    Rcpp::Rcout << "\n";
}
