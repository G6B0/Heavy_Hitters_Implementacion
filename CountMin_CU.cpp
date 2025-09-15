#include "CountMin_CU.h"


CountMin_CU::CountMin_CU(int d,int w):
d(d), w(w), C(d,std::vector<int>(w)){}

void CountMin_CU::insertCMin(uint64_t elemento){
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

double CountMin_CU::estimar_freq(uint64_t elemento) const {
    int freq_est = INT_MAX;
    for(int i = 0; i < d; i++){
        uint32_t h = murmurhash(&elemento, i) % w;
        if(C[i][h] < freq_est) freq_est = C[i][h];
    }
    return freq_est;
}



