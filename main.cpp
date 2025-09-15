#include <iostream>
#include "TowerSketch.h"

int main (int argc, char *argv[]) {
    TowerSketch TowerSK(5);
    for(int j=0; j < 100; j++){
        for(int i=0; i<j; i++){
            TowerSK.insert(j);
        }
    }
    int x;
    std::cin >> x;
    std::cout << "La frecuencia del elemento " << x << " es: " << TowerSK.estimar_freq(x) <<
    ".\n";
    return 0;
}

//int main (int argc, char *argv[]) {
//    CountMin_CU  CU(5, 300);
//    for(int j=0; j < 100; j++){
//        for(int i=0; i<j; i++){
//            CU.insertCMin(j);
//        }
//    }
//    int x;
//    std::cin >> x;
//    std::cout << "La frecuencia del elemento " << x << " es: " << CU.estimar_freq(x) <<
//    ".\n";
//    return 0;
//}