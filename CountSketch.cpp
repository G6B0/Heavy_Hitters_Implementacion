#include "murmurhash32.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <limits.h>
#include <algorithm>


double median(std::vector<int> m){
    int n = m.size();
    int mitad = n/2;
    std::nth_element(m.begin(), m.begin() + mitad, m.end());
    if(n%2 == 1){
        return m[mitad]; //caso impar
    }
    else{
        int mitad1 = m[mitad];
        int max = *std::max_element(m.begin(), m.begin() + mitad);
        return (mitad1 + max) / 2.0;
    }

}
class CountSketch
{
private:
    int c = 1; //constante
    int d;
    int w;
    std::vector<std::vector<int>> C;
    std::vector<int> g;

    std::vector<int> init_g(int tamaño){
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

public:
    CountSketch(int d, int w):
    d(d), w(w), C(d,std::vector<int>(w)), g(init_g(d)){

    }
    void insert(uint64_t elemento){
        for( int i = 0; i < d ; i++){
            uint32_t h = murmurhash(&elemento, i) % w;
            C[i][h] += g[i]*c;
        }
    }
    double estimar_freq(uint64_t elemento){
        std::vector<int> estimaciones(d);
        for (int i = 0; i < d; i++){
            uint32_t h = murmurhash(&elemento, i) % w;
            estimaciones[i] = g[i]*C[i][h];

        }
        return median(estimaciones);
    }
};

int main(int argc, char const *argv[])
{
    CountSketch CK(10, 300);

    for(int j=0; j<100 ; j++){
        for(int i=0; i<j; i++){
            CK.insert(j);
        }
    }
    int x;
    std::cin >> x;
    std::cout << "la frecuencia del elemento " << x << " es: " << CK.estimar_freq(x) << ".\n";

    return 0;
}
