#ifndef TOWERSKETCH_H
#define TOWERSKETCH_H

#include "CountMin_CU.h"

class TowerSketch
{
    private:
        int niveles;
        std::vector<CountMin_CU> Tower;
        void tower_create(int niveles);
    public:
        TowerSketch(int niveles);
        void insert(uint64_t elemento);
        double estimar_freq(uint64_t elemento);
        size_t get_memory_size() const;
};

#endif