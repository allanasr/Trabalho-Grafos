#ifndef MIN_HEAP_H
#define MIN_HEAP_H

#include "MinHeapNode.h"

class MinHeap
{
    MinHeapNode **heap;
    int capacity;
    int size;
public:
    MinHeap(int capacity);

    ~MinHeap();

    void MinHeapify(int);

    int parent(int i)
    {
        return (i - 1) / 2;
    }

    int left(int i)
    {
        return (2 * i + 1);
    }

    int right(int i)
    {
        return (2 * i + 2);
    }

    MinHeapNode *extractMin();

    void decreaseKey(int i, int new_val);

    MinHeapNode *getMin()
    {
        return heap[0];
    }

    void deleteKey(int i);

    void insertKey(MinHeapNode *k);

    bool isEmpty()
    {
        return size == 0;
    }

    int getIndexOf(int id);

    void clear();
};
#endif
