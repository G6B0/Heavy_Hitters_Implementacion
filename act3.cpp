#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <filesystem>
#include <cstdlib>
#include "TowerSketch.h"
#include "CountSketch.h"
#include "Extraccionkmer.h"
#include <fstream>


// Calcula el reverso complementario de un k-mer en uint64_t
uint64_t reverse_complement_uint64(uint64_t kmer_val, int k) {
    uint64_t rc_val = 0;
    for (int i = 0; i < k; ++i) {
        uint64_t base = (kmer_val >> (2 * i)) & 0b11;
        rc_val = (rc_val << 2) | (base ^ 0b11);
    }
    return rc_val;
}


// Función para calcular memoria usada por heavy hitters
size_t calculate_hh_memory(const std::unordered_map<uint64_t, int>& kmer_frequencies, 
                          long long frequency_threshold, int k) {
    size_t memory_bytes = 0;
    
    for (const auto& pair : kmer_frequencies) {
        if (pair.second >= frequency_threshold) {
            // Memoria por cada HH:
            // - uint64_t key: 8 bytes
            // - int frequency: 4 bytes  
            // - overhead del unordered_map: ~24 bytes por entrada
            // - string del k-mer para almacenamiento: k bytes
            memory_bytes += sizeof(uint64_t) + sizeof(int) + 24 + k;
        }
    }
    return memory_bytes;
}

std::pair<std::vector<std::pair<uint64_t, bool>>, long long> leer_subgenoma(const std::string& genome_filename, int SUBSET_SIZE, long long N_total_genoma) {
    std::ifstream genomeFile(genome_filename);

    std::string sequence;
    size_t bases_read = 0;
    std::string line;

    while (std::getline(genomeFile, line) && bases_read < SUBSET_SIZE) {
        if (line.empty() || line[0] == '>') continue;
        for (char base : line) {
            if (bases_read >= SUBSET_SIZE) break;
            if (base == 'A' || base == 'C' || base == 'G' || base == 'T') {
                sequence += base;
                bases_read++;
            }
        }
    }
    genomeFile.close();

    const size_t MEMORY_BYTES = 1024 * 1024;

    double phi = 0.0000002; //phi para sacar los HH (bueno: 0.000022 con 100000 bases)
    long long frequency_threshold = static_cast<long long>(phi * N_total_genoma);
    
    size_t total_hh_memory = 0;
    int total_hh_count = 0;
    std::vector<std::pair<uint64_t, bool>> subconjunto;

    // extrae k-mers y escribe en output file
    for (int k : {21, 31}) {
        // Usa un mapa de uint64_t para las frecuencias
        std::unordered_map<uint64_t, int> kmer_frequencies;
        
        if (sequence.length() < k) continue;

        for (size_t i = 0; i <= sequence.length() - k; ++i) {
            uint64_t kmer_val = string_to_uint64(sequence.substr(i, k));
            uint64_t rc_kmer_val = reverse_complement_uint64(kmer_val, k);
            uint64_t canonical_kmer_val = std::min(kmer_val, rc_kmer_val);
            kmer_frequencies[canonical_kmer_val]++;
        }

        

        // Calcular memoria de HH para este k
        size_t hh_memory = calculate_hh_memory(kmer_frequencies, frequency_threshold, k);
        total_hh_memory += hh_memory;
        int hh_count = 0;

        for (const auto& pair : kmer_frequencies) {
            bool is_hh = (pair.second >= frequency_threshold);
            subconjunto.push_back({pair.first,is_hh});
            if (is_hh) {
                hh_count++;
                total_hh_count++;
            }
        }
    }

    std::cout << "Total HH: " << total_hh_count << std::endl;
    std::cout << "Memoria total HH: " << total_hh_memory/1024.0 << " KB" << std::endl;
    return {subconjunto, frequency_threshold};
}

int main(int argc, char* argv[]) {

    // Tabla de conversión para la codificación y su complemento
    uint8_t base_to_int[128];
    const char int_to_base[4] = {'A', 'C', 'G', 'T'};
    initialize_base_table(base_to_int);
    TowerSketch towersketch(8);
    CountSketch countsketch(15,500000);
    
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
    std::cout << "Secuencia completa leída: " << sequence.length() << " bases" << std::endl;


    long long N = 0;
    for (int k : {21, 31}) {
        if (sequence.length() >= k) {
            N += sequence.length() - k + 1;
        }
    }
    std::cout << "N total (k-mers en genoma completo): " << N << std::endl;

    std::string archivo_csv = "resultados_totales.csv";
    int subconjunto_length = 2000000;
    std::pair<std::vector<std::pair<uint64_t,bool>>, long long> subconjunto = leer_subgenoma(genome_filename, subconjunto_length, N);
    std::vector<std::pair<uint64_t,bool>> kmers = subconjunto.first;
    long long threshold_archivo = subconjunto.second;


    
    // extrae k-mers e inserta en los sketchs
    for (int k : {21, 31}) {
        for (size_t i = 0; i <= sequence.length() - k; ++i) {

            uint64_t kmer_val = string_to_uint64(sequence.substr(i, k));
            uint64_t rc_kmer_val = reverse_complement_uint64(kmer_val, k);
            uint64_t canonical_kmer_val = std::min(kmer_val, rc_kmer_val);
            countsketch.insert(canonical_kmer_val);
            towersketch.insert(canonical_kmer_val);
            N++;
        }
    }

    // Analisis del ground truth
int real_HH_count = 0;
double phi = 0.0000002;
for (const auto& kmer_pair : kmers) {
   if (kmer_pair.second) real_HH_count++;
}

std::cout << "=== GROUND TRUTH ANALYSIS ===" << std::endl;
std::cout << "Total k-mers en subconjunto: " << kmers.size() << std::endl;
std::cout << "Heavy hitters reales: " << real_HH_count << std::endl;

int TP = 0, FP = 0, FN = 0, TN = 0;
int TP1 = 0, FP1 = 0, FN1 = 0, TN1 = 0;

// CountSketch analysis
for (size_t i = 0; i < kmers.size(); i++) {
    double freq_est = countsketch.estimar_freq(kmers[i].first);
    bool es_HH_estimado = (freq_est >= threshold_archivo);
    bool es_HH_real = kmers[i].second;

    if (es_HH_estimado && es_HH_real) TP++;
    else if (es_HH_estimado && !es_HH_real) FP++;
    else if (!es_HH_estimado && es_HH_real) FN++;
    else TN++;
}

// TowerSketch analysis
for (size_t i = 0; i < kmers.size(); i++) {
    double freq_est = towersketch.estimar_freq(kmers[i].first);
    bool es_HH_estimado = (freq_est >= threshold_archivo);
    bool es_HH_real = kmers[i].second;

    if (es_HH_estimado && es_HH_real) TP1++;
    else if (es_HH_estimado && !es_HH_real) FP1++;
    else if (!es_HH_estimado && es_HH_real) FN1++;
    else TN1++;
}

// Calcular métricas
double precision = (TP + FP > 0) ? (double)TP / (TP + FP) : 0.0;
double recall    = (TP + FN > 0) ? (double)TP / (TP + FN) : 0.0;
double f1        = (precision + recall > 0) ? 2 * precision * recall / (precision + recall) : 0.0;

double precision1 = (TP1 + FP1 > 0) ? (double)TP1 / (TP1 + FP1) : 0.0;
double recall1    = (TP1 + FN1 > 0) ? (double)TP1 / (TP1 + FN1) : 0.0;
double f11        = (precision1 + recall1 > 0) ? 2 * precision1 * recall1 / (precision1 + recall1) : 0.0;

double memory = countsketch.get_size_MB();
double memory1 = towersketch.get_memory_size();

// Imprimir en consola (mantener output actual)
std::cout << "CountSketch:" << std::endl;
std::cout << "Precision: " << precision << std::endl;
std::cout << "Recall: " << recall << std::endl;
std::cout << "F1-score: " << f1 << std::endl;
std::cout << "TP: " << TP << ", FP: " << FP << ", FN: " << FN << ", TN: " << TN << std::endl;
std::cout << "Memory: " << memory << "MB" << std::endl;

std::cout << "\nTowerSketch:" << std::endl;
std::cout << "Precision: " << precision1 << std::endl;
std::cout << "Recall: " << recall1 << std::endl;
std::cout << "F1-score: " << f11 << std::endl;
std::cout << "TP: " << TP1 << ", FP: " << FP1 << ", FN: " << FN1 << ", TN: " << TN1 << std::endl;
std::cout << "Memory: " << memory1 << "MB" << std::endl;

// Escribir en CSV
std::string csv_filename = "resultados_experimento.csv";
std::ofstream csvFile;

// Verificar si el archivo ya existe para decidir si escribir headers
bool file_exists = std::filesystem::exists(csv_filename);
csvFile.open(csv_filename, std::ios::app); // Modo append

// Escribir headers solo si es un archivo nuevo
if (!file_exists) {
    csvFile << "Genoma,Total_kmers,Real_HH,Phi,Threshold,";
    csvFile << "Algoritmo,Precision,Recall,F1_Score,TP,FP,FN,TN,Memory_MB\n";
}

// Extraer nombre del genoma del path
std::string genome_name = genome_filename;
size_t last_slash = genome_name.find_last_of("/\\");
if (last_slash != std::string::npos) {
    genome_name = genome_name.substr(last_slash + 1);
}

// Escribir datos de CountSketch
csvFile << genome_name << "," << kmers.size() << "," << real_HH_count << ",";
csvFile << phi << "," << threshold_archivo << ",";
csvFile << "CountSketch," << precision << "," << recall << "," << f1 << ",";
csvFile << TP << "," << FP << "," << FN << "," << TN << "," << memory << "\n";

// Escribir datos de TowerSketch
csvFile << genome_name << "," << kmers.size() << "," << real_HH_count << ",";
csvFile << phi << "," << threshold_archivo << ",";
csvFile << "TowerSketch," << precision1 << "," << recall1 << "," << f11 << ",";
csvFile << TP1 << "," << FP1 << "," << FN1 << "," << TN1 << "," << memory1 << "\n";

csvFile.close();
std::cout << "\nResultados guardados en: " << csv_filename << std::endl;

return 0;

}