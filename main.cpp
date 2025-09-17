#include <iostream>
#include "TowerSketch.h"
#include "Extraccionkmer.h"

//int main (int argc, char *argv[]) {
//    TowerSketch TowerSK(5);
//    for(int j=0; j < 100; j++){
//        for(int i=0; i<j; i++){
//            TowerSK.insert(j);
//        }
//    }
//    int x;
//    std::cin >> x;
//    std::cout << "La frecuencia del elemento " << x << " es: " << TowerSK.estimar_freq(x) <<
//    ".\n";
//    return 0;
//}

int main (int argc, char *argv[]) {
    TowerSketch TW(5);
    std::string archivo_csv = "resultados_totales.csv";
    std::vector<uint64_t> kmers = leer_kmers(archivo_csv);
    for(size_t i=0; i < kmers.size(); i++){
        TW.insert(kmers[i]);
        std::cout << "insertando kmer de la posición" << i << "\n";
    }
    std::string x;
    std::cin >> x;
    std::cout << "La frecuencia del elemento " << x << " es: " << TW.estimar_freq(string_to_uint64(x)) <<
    ".\n";
    return 0;
}