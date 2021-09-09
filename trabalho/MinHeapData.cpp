#include <climits>
#include <iostream>
#include "MinHeapData.h"
#include "MinHeapNode.h"

using  namespace std;

void swap(MinHeapNode **x, MinHeapNode **y)
{
    MinHeapNode *temp = *x;
    *x = *y;
    *y = temp;
}

MinHeap::MinHeap(int capacity)
{

    size = 0;
    this->capacity = capacity;
    heap = new MinHeapNode *[capacity];
}

MinHeap::~MinHeap()
{
    for (int i = 0; i < size; i++)
    {
        delete heap[i];
    }
    delete[] heap;
}

void MinHeap::insertKey(MinHeapNode *k)
{
    if (size == capacity)
    {
        std::cout << "Heap is full, key not inserted" << std::endl;
        return;
    }

    int i = size;
    heap[i] = k;
    size++;

    while (i != 0 && heap[parent(i)]->getWeight() > heap[i]->getWeight())
    {
        swap(&heap[i], &heap[parent(i)]);
        i = parent(i);
    }
}

void MinHeap::decreaseKey(int i, int new_val)
{
    heap[i]->setWeight(new_val);
    while (i != 0 && heap[parent(i)]->getWeight() > heap[i]->getWeight())
    {
        swap(&heap[i], &heap[parent(i)]);
        i = parent(i);
    }
}

MinHeapNode *MinHeap::extractMin()
{
    if (size <= 0)
        return nullptr;
    if (size == 1)
    {
        size--;
        return heap[0];
    }

    MinHeapNode *root = heap[0];
    heap[0] = heap[size - 1];
    size--;
    MinHeapify(0);

    return root;
}

void MinHeap::deleteKey(int i)
{
    decreaseKey(i, INT_MIN);
    extractMin();
}

void MinHeap::MinHeapify(int i)
{
    int l = left(i);
    int r = right(i);
    int smallest = i;
    if (l < size && heap[l] < heap[i])
        smallest = l;
    if (r < size && heap[r] < heap[smallest])
        smallest = r;
    if (smallest != i)
    {
        swap(&heap[i], &heap[smallest]);
        MinHeapify(smallest);
    }
}

int MinHeap::getIndexOf(int id)
{
    for (int i = 0; i < size; i++)
    {
        if (heap[i]->getId() == id)
            return i;
    }
    return -1;
}

void MinHeap::clear()
{
    for (int i = 0; i < size; i++)
    {
        delete heap[i];
    }
    size = 0;
}
