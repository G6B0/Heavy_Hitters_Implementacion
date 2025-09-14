#include "murmurhash32.hpp"
#include <iostream>
#include <limits.h>
#include <vector>

class CountMin_CU
{
private:
    int d;
    int w;
    std::vector<std::vector<int>> C;
public:
    CountMin_CU(int d,int w):
    d(d), w(w), C(d,std::vector<int>(w)){

    }
    void insertCMin(uint64_t elemento){
        std::vector<uint32_t> hj(d); 
        int min = INT_MAX;
        for(int i = 0; i < d; i++){
            uint32_t h = murmurhash(&elemento, i) % w;
            hj[i] = h;
            if(C[i][h] < min){
                min = C[i][h];
            }
        }
        for(int i = 0; i < d; i++){
            if(C[i][hj[i]] == min){
                C[i][hj[i]]++;
            }
        }
    }

    double estimar_freq(uint64_t elemento){
        int freq_est = INT_MAX;
        for(int i = 0; i < d; i++){
            uint32_t h = murmurhash(&elemento, i) % w;
            if(C[i][h] < freq_est) freq_est = C[i][h];
        }
        return freq_est;
    }
};


int main (int argc, char *argv[]) {
    CountMin_CU CM(10, 300);
    for(int j=0; j < 100; j++){
        for(int i=0; i<j; i++){
            CM.insertCMin(j);
        }
    }
    int x;
    std::cin >> x;
    std::cout << "La frecuencia del elemento " << x << " es: " << CM.estimar_freq(x) <<
    ".\n";
    return 0;
}



