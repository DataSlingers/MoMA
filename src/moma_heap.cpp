#include "moma_heap.h"

/*
* Create a full order at nodes of the form <lambda,id>
*/
bool gt(const HeapNode &left, const HeapNode &right){
    return left.lambda > right.lambda;
}

/*
* Find the smaller child of an element in a heap. Used in siftdown
*/
int min_child(std::vector<HeapNode> &h, int i) {
    int cur_size = h.size();
	int child = i * 2 + 1;
	if (child >= cur_size) {
		// no children
		return NO_CHILED;
	} else if (child+1 >= cur_size || !gt(h[child],h[child+1])){
		// only child or first child is biggest child
		return child;
	} else {
		// second child exists and is smallest child
		return child+1;
	}
}

/*
* TODO: extra copy can be avoided in siftdown
*/
void swap(std::vector<HeapNode> &h,int i, int j){
    Rcpp::Rcout << h[i].lambda << "and" << h[j].lambda << "swapped\n";
    HeapNode tmp = h[i];
    h[i] = h[j];
    h[j] = tmp;
}

/*
* In a min-heap, if the key (lambda in our case) decreases, sift it up
*/
void siftup(std::vector<HeapNode> &heap, int i) {
    HeapNode tmp = heap[i];
    int parent = (i - 1) / 2;
    while (i != 0 && !gt(heap[i],heap[parent])) {
        heap[i] = heap[parent];
        i = parent;
        parent = (i - 1) / 2;
    }
    heap[i] = tmp;
}

/*
* In a min-heap, if the key (lambda in our case) increases, sift it down
*/
void siftdown(std::vector<HeapNode> &h, int current_node) {
    int cur_size = h.size();
	int child = min_child(h, current_node);
	while (child < cur_size && gt(h[current_node], h[child])){
		swap(h, child, current_node);
        current_node = child;
		child = min_child(h, child);
	}
}

/*
* Change the key of any nodes; TODO: not use trasversal
*/
int heap_change_lambda(std::vector<HeapNode> &heap, int id, double new_lambda){
    for(int i = 0; i < heap.size(); i++){
        if(heap[i].id == id){
            double old_lambda = heap[i].lambda;
            heap[i].lambda = new_lambda;
            if(old_lambda < new_lambda)
                siftdown(heap,i);
            else
                siftup(heap,i);
            return i;
        }
    }
    MoMALogger::error("No such id in current heap: ")<<id;
    return -1;
}

/*
* To delete an element, move it to the tail, pop it out, and then sift down 
* the node that replaces it
*/
void heap_delete(std::vector<HeapNode> &heap, int id){
    for(int i = 0; i < heap.size(); i++){
        if(heap[i].id == id){
            swap(heap,i,heap.size()-1);
            heap.pop_back();
            siftdown(heap,i);
            return;
        }
    }
    //MoMALogger::error("No such id in current heap: ")<<id;
}

/*
* Check if an array is a min heap
*/
bool is_minheap(std::vector<HeapNode> &heap){
    int i = 0;
    while(2 * i + 1 < heap.size()){
        if(gt(heap[i],heap[2 * i + 1])){
            Rcpp::Rcout << heap[i].lambda << "and"<< heap[2*i+1].lambda;
            return 0;
        }
        if(2 * i + 2 < heap.size()){
            if(gt(heap[i],heap[2 * i + 2])){
                Rcpp::Rcout << heap[i].lambda << "and"<< heap[2*i+2].lambda;
                return 0;
            }
        }
        i++;
    }
    return 1;
}

bool is_empty(std::vector<HeapNode> &heap){
    return heap.size() == 0;
}
/*
* Get the currently minimun value without deleting the node
*/
HeapNode heap_peek_min(std::vector<HeapNode> &heap){
    if(heap.size() == 0){
        MoMALogger::error("Empty heap!");
    }
    HeapNode n = heap.front();
    return n;
}

/*
* Print the heap
*/
void heap_print(const std::vector<HeapNode> &q){
    Rcpp::Rcout << "Heap lambda is\n";
    for (auto i : q) Rcpp::Rcout << i.lambda << "\t";
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "Heap id is\n";
    for (auto i : q) Rcpp::Rcout << i.id << "\t";
    Rcpp::Rcout << "\n";
}
