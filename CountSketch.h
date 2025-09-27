#ifndef COUNTSKETCH_H
#define COUNTSKETCH_H

#include <vector>
#include <cstdint>

class CountSketch
{
    private:
    int d;
    int w;
    int c = 1;
    std::vector<std::vector<int64_t>> C;
    std::vector<int> g;
    std::vector<int> init_g(int tamaño);

    public:
    CountSketch(int d, int w);

    void insert(uint64_t elemento);
    double estimar_freq(uint64_t elemento);
    double get_size_MB() const;
};

#endif