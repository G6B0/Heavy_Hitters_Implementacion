#ifndef COUNTMIN_CU_H
#define COUNTMIN_CU_H

#include <vector>
#include <cstdint>
#include <climits>

class CountMin_CU
{
    private:
    int d;
    int w;
    std::vector<std::vector<uint64_t>> C;

    public:
    CountMin_CU(int d, int w);

    void insertCMin(uint64_t elemento);
    double estimar_freq(uint64_t elemento) const;
    size_t get_memory_size() const;
};

#endif