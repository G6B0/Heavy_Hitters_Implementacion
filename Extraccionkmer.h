#ifndef EXTRACCIONKMER_H
#define EXTRACCIONKMER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

//Inicializa base_table (Hecho en el act1.cpp)
void initialize_base_table(uint8_t base_to_int[128]) {
    base_to_int['A'] = 0; base_to_int['a'] = 0; // A: 00
    base_to_int['C'] = 1; base_to_int['c'] = 1; // C: 01
    base_to_int['G'] = 2; base_to_int['g'] = 2; // G: 10
    base_to_int['T'] = 3; base_to_int['t'] = 3; // T: 11
}

// Convierte un uint64_t a un k-mer de cadena (para csv) (Hecho en el act1.cpp)
std::string uint64_to_string(uint64_t kmer_val, int k) {
    std::string result(k, ' ');
    const char int_to_base[4] = {'A', 'C', 'G', 'T'};
    for (int i = 0; i < k; ++i) {
        uint64_t base_bits = (kmer_val >> (2 * i)) & 0b11;
        result[k - 1 - i] = int_to_base[base_bits];
    }
    return result;
}

// Convierte un k-mer de cadena a uint64_t (Hecho en el act1.cpp)
uint64_t string_to_uint64(const std::string& kmer) {
    uint8_t base_to_int[128];
    initialize_base_table(base_to_int);
    uint64_t result = 0;
    for (char base : kmer) {
        result = (result << 2) | base_to_int[static_cast<uint8_t>(base)];
    }
    return result;
}

uint64_t encode_kmer_con_k(uint64_t kmer_val, int k){
    if(k == 21){
        return kmer_val << 2;
    }
    else if (k == 31){
        return (kmer_val << 2) | 0b01;
    }
    else {
        throw std::runtime_error ("k invalido, solo se permite 21 o 31");
    }
}

void decode_kmer_con_k(uint64_t encoded, uint64_t &kmer_val, int &k){
    int k_bits = encoded & 0b11;
    kmer_val = encoded >> 2;
    if(k_bits == 0b00){
        k=21;
    }
    else if(k_bits == 0b01){
        k = 31;
    }
    else{
        throw std::runtime_error("k codificado invalido");
    }
}

//Lee los kmers que estan guardados en el .csv (obtenidos por act1.cpp)
std::vector<uint64_t> leer_kmers(const std::string& csv_filename){
    std::vector<uint64_t> kmers;
    std::ifstream file(csv_filename);

    if (!file.is_open()){
        std::cerr << "Error: nose pudo abrir el archivo" << csv_filename << "\n";
        return kmers;
    }

    std::string line;
    bool first_line = true;
    while (std::getline(file, line)){
        if (first_line){
            first_line = false;
            continue;
        }

        std::stringstream ss(line);
        std::string campo;
        int col_x = 0;
        std::string kmer_str;
        int k_value = 0;

        while (std::getline(ss, campo, ',')){
            if (col_x == 1){
                kmer_str = campo;
            }
            if (col_x == 3){
                k_value = std::stoi(campo);
            }
            col_x++;
        }
        
        if(!kmer_str.empty() && (k_value == 21 || k_value == 31)){
            uint64_t kmer_val = string_to_uint64(kmer_str);
            kmers.push_back(encode_kmer_con_k(kmer_val, k_value));
        }
    }
    
    file.close();
    return kmers;
}

#endif