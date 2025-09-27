#include "murmurhash32.hpp"
#include <random>
#include <limits.h>
#include <algorithm>
#include "CountSketch.h"



CountSketch::CountSketch(int d, int w):
d(d), w(w), c(1), C(d,std::vector<int64_t>(w)), g(init_g(d)){}

double median(std::vector<int64_t> m) {
    int n = m.size();
    int mitad = n/2;
    std::nth_element(m.begin(), m.begin() + mitad, m.end());
    if(n%2 == 1){
        return m[mitad]; //caso impar
    }
    else{
        int64_t elemento_mitad = m[mitad];
        std::nth_element(m.begin(), m.begin() + mitad - 1, m.begin() + mitad);
        int64_t elemento_anterior = m[mitad - 1];
        return (elemento_mitad + elemento_anterior) / 2.0;
    }
}

std::vector<int> CountSketch::init_g(int tamaño){
    std::vector<int> g(tamaño);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0,1);

    for(int i = 0; i < tamaño; i++){
        g[i] = dist(gen);
        if(g[i] == 0){
            g[i] = -1;
        }
        else{
            g[i] = 1;
        }
    }
    return g;
};

void CountSketch::insert(uint64_t elemento){
    for( int i = 0; i < d ; i++){
        uint32_t h = murmurhash(&elemento, i) % w;
        C[i][h] += g[i]*c;
    }
}

double CountSketch::estimar_freq(uint64_t elemento){
    std::vector<int64_t> estimaciones(d);
    for (size_t i = 0; i < d; i++){
        uint32_t h = murmurhash(&elemento, i) % w;
        int64_t kmer = C[i][h];
        estimaciones[i] = g[i] * kmer;
    }
    return median(estimaciones);
}

double CountSketch::get_size_MB() const {
    size_t sizeC = 0;
    for (const auto& row : C)
        sizeC += row.capacity() * sizeof(uint64_t);

    size_t sizeG = g.capacity() * sizeof(int);

    // bytes → MB
    return static_cast<double>(sizeC + sizeG) / (1024.0 * 1024.0);
}

