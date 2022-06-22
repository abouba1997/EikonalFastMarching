#include "Heap.h"

void swap(double* _first, double* _second)
{
    double* temp = _first;
    *_first = *_second;
    *_second = *temp;
}


Heap::Heap(int cap)
{
    size = 0;
    capacity = cap;
    heap.resize(capacity);
}

void Heap::insert(double k) {
    if (size == capacity) {
        std::cout << "MIN HEAP FULL" << std::endl;
        return;
    }

    size++;
    int i = size - 1;
    heap[i] = k;

    while (i != 0 && heap[parent(i)] > heap[i]) {
        std::swap(heap[i], heap[parent(i)]);
        i = parent(i);
    }
}

double Heap::extractMin() {
    if (size == 0) {
        std::cout << "EMPTY HEAP" << std::endl;
        return -1;
    }
    else if (size == 1) {
        size--;
        return heap[0];
    }
    else {
        double root = heap[0];
        heap[0] = heap[size - 1];
        size--;
        heapify(0);
        return root;
    }
}

void Heap::heapify(int i)
{
    int l = left(i);
    int r = right(i);
    int smallest = i;

    if ((l < size) && (heap[l] < heap[smallest])) {
        smallest = l;
    }
    if ((r < size) && (heap[r] < heap[smallest])) {
        smallest = r;
    }

    if (smallest != i) {
        swap(heap[i], heap[smallest]);
        heapify(smallest);
    }
}

void Heap::printHeap()
{
    int power = 0;
    int value = 1;

    for (int i = 0; i < size; i++)
    {
        if (i == value) {
            std::cout << std::endl;
            power++;
            value += (1 << power);
        }
        std::cout << heap[i] << "\t";
    }
    std::cout << std::endl;
}