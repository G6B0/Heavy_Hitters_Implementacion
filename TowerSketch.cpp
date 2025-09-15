#include "TowerSketch.h"

TowerSketch::TowerSketch(int niveles):
    niveles(niveles){
        tower_create(niveles);
}

void TowerSketch::tower_create(int niveles){
    int j = 2;
    Tower.clear();
    Tower.reserve(niveles);
    for(int i = 0; i < niveles; i++){
        Tower.emplace_back(5,j);
        j *= 2; 
    }
}
    
void TowerSketch::insert(uint64_t elemento){
    for(int i = 0; i < niveles; i++){
        Tower[i].insertCMin(elemento);
    }
}
    
double TowerSketch::estimar_freq(uint64_t elemento){
    int freq_est = INT_MAX;
    for(int i = 0; i < niveles; i++){
        int val = Tower[i].estimar_freq(elemento); 
        if(val < freq_est){
            freq_est = val;
        }
    }
    return freq_est;
}