#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <filesystem>
#include <cstdlib>
#include "TowerSketch.h"
#include "Extraccionkmer.h"

// Calcula el reverso complementario de un k-mer en uint64_t
uint64_t reverse_complement_uint64(uint64_t kmer_val, int k) {
    uint64_t rc_val = 0;
    for (int i = 0; i < k; ++i) {
        uint64_t base = (kmer_val >> (2 * i)) & 0b11;
        rc_val = (rc_val << 2) | (base ^ 0b11);
    }
    return rc_val;
}


int main(int argc, char* argv[]) {

    // Tabla de conversión para la codificación y su complemento
    uint8_t base_to_int[128];
    const char int_to_base[4] = {'A', 'C', 'G', 'T'};
    initialize_base_table(base_to_int);
    TowerSketch towersketch(10);
    
    if (argc < 2) {
        std::cerr << "Uso: " << argv[0] << " <archivo_genoma>\n";
        return 1;
    }


    std::string genome_filename = argv[1];
    std::string sequence;
    
    // lee genoma input file
    std::ifstream genomeFile(genome_filename);
    if (!genomeFile.is_open()) {
        std::cerr << "Error: No se pudo abrir el archivo " << genome_filename << "\n";
        return 1;
    }

    int bases_read = 0;
    std::string line;
    while (std::getline(genomeFile, line)) {
        if (line.empty() || line[0] == '>') continue;
        for (char base : line) {
            if (base == 'A' || base == 'C' || base == 'G' || base == 'T') {
                sequence += base;
                bases_read++;
            }
        }
    }
    genomeFile.close();

    double phi = 0.000022; //phi para sacar los HH (bueno: 0.000022 con 100000 bases)
    uint64_t N = 0;
    std::string archivo_csv = "resultados_totales.csv";
    std::vector<std::pair<uint64_t,bool>> kmers = leer_kmers(archivo_csv, genome_filename);

    // extrae k-mers e inserta en los sketchs
    for (int k : {21, 31}) {

        for (size_t i = 0; i <= sequence.length() - k; ++i) {
            uint64_t kmer_val = string_to_uint64(sequence.substr(i, k));
            uint64_t rc_kmer_val = reverse_complement_uint64(kmer_val, k);
            uint64_t canonical_kmer_val = std::min(kmer_val, rc_kmer_val);
            towersketch.insert(canonical_kmer_val);
            N++;
        }
    }

    std::vector<uint64_t> heavy_hitters;
    

    int TP = 0, FP = 0, FN = 0, TN = 0;

for (size_t i = 0; i < kmers.size(); i++) {
    int k;
    uint64_t kmer_canonico;
    decode_kmer_con_k(kmers[i].first, kmer_canonico, k);

    int freq_est = towersketch.estimar_freq(kmer_canonico);
    bool es_HH_estimado = (freq_est >= phi * N);
    bool es_HH_real = kmers[i].second;

    if (es_HH_estimado && es_HH_real) TP++;
    else if (es_HH_estimado && !es_HH_real) FP++;
    else if (!es_HH_estimado && es_HH_real) FN++;
    else TN++;
}

double precision = (TP + FP > 0) ? (double)TP / (TP + FP) : 0.0;
double recall    = (TP + FN > 0) ? (double)TP / (TP + FN) : 0.0;
double f1        = (precision + recall > 0) ? 2 * precision * recall / (precision + recall) : 0.0;

std::cout << "Precision: " << precision << "\n";
std::cout << "Recall: " << recall << "\n";
std::cout << "F1-score: " << f1 << "\n";

    return 0;
}