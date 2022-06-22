#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <string>
#include <ctime>

using namespace std;

class Heap {
private:
    int size;
    int capacity;

    int parent(int i) { return (i - 1) / 2; }
    int left(int i) { return i * 2 + 1; }
    int right(int i) { return i * 2 + 2; }

public:
    Heap(int cap);
    void insert(double k);
    double extractMin();
    void heapify(int i);
    void printHeap();
    vector<double> heap;
};
