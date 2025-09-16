#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <filesystem>

// Tabla de conversión para la codificación y su complemento
uint8_t base_to_int[128];
const char int_to_base[4] = {'A', 'C', 'G', 'T'};

void initialize_base_table() {
    base_to_int['A'] = 0; base_to_int['a'] = 0; // A: 00
    base_to_int['C'] = 1; base_to_int['c'] = 1; // C: 01
    base_to_int['G'] = 2; base_to_int['g'] = 2; // G: 10
    base_to_int['T'] = 3; base_to_int['t'] = 3; // T: 11
}

// Convierte un k-mer de cadena a uint64_t
uint64_t string_to_uint64(const std::string& kmer) {
    uint64_t result = 0;
    for (char base : kmer) {
        result = (result << 2) | base_to_int[static_cast<uint8_t>(base)];
    }
    return result;
}

// Convierte un uint64_t a un k-mer de cadena (para csv)
std::string uint64_to_string(uint64_t kmer_val, int k) {
    std::string result(k, ' ');
    for (int i = 0; i < k; ++i) {
        uint64_t base_bits = (kmer_val >> (2 * i)) & 0b11;
        result[k - 1 - i] = int_to_base[base_bits];
    }
    return result;
}

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

int main(int argc, char* argv[]) {
    initialize_base_table();
    
    if (argc < 3) {
        std::cerr << "Uso: " << argv[0] << " <archivo_genoma> <archivo_salida_csv>\n";
        return 1;
    }

    std::string genome_filename = argv[1];
    std::string output_filename = argv[2];
    std::string sequence;
    
    // lee genoma input file
    std::ifstream genomeFile(genome_filename);
    if (!genomeFile.is_open()) {
        std::cerr << "Error: No se pudo abrir el archivo " << genome_filename << "\n";
        return 1;
    }

    const int SUBSET_SIZE = 100000; // Porción de 100k bases
    int bases_read = 0;
    std::string line;
    while (std::getline(genomeFile, line) && bases_read < SUBSET_SIZE) {
        if (line.empty() || line[0] == '>') continue;
        for (char base : line) {
            if (bases_read >= SUBSET_SIZE) break; // Detener si se alcanza el límite
            if (base == 'A' || base == 'C' || base == 'G' || base == 'T') {
                sequence += base;
                bases_read++;
            }
        }
    }
    genomeFile.close();

    // abre output file en append
    std::ofstream outputFile(output_filename, std::ios::app);
    if (!outputFile.is_open()) {
        std::cerr << "Error: No se pudo abrir o crear el archivo de salida " << output_filename << "\n";
        return 1;
    }
    
    if (outputFile.tellp() == 0) {
        outputFile << "Genoma,kmer,frecuencia,k_value,is_heavy_hitter\n";
    }


    double phi = 0.000022; //phi para sacar los HH


    size_t total_hh_memory = 0;
    int total_hh_count = 0;

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

        long long total_kmers = 0;
        for (const auto& pair : kmer_frequencies) {
            total_kmers += pair.second;
        }

        long long frequency_threshold = static_cast<long long>(phi * total_kmers);

        // Calcular memoria de HH para este k
        size_t hh_memory = calculate_hh_memory(kmer_frequencies, frequency_threshold, k);
        total_hh_memory += hh_memory;
        int hh_count = 0;

        for (const auto& pair : kmer_frequencies) {
            bool is_hh = (pair.second >= frequency_threshold);
            if (is_hh) {
                hh_count++;
                total_hh_count++;
            }
            outputFile << genome_filename << "," << uint64_to_string(pair.first, k) << "," << pair.second << "," << k << "," << (is_hh ? "true" : "false") << "\n";
        }
    }

    outputFile.close();
    std::cout << "Resultados de " << genome_filename << " guardados en " << output_filename << std::endl;

    std::cout << "Total HH: " << total_hh_count << std::endl;
    std::cout << "Memoria total HH: " << total_hh_memory/1024.0 << " KB" << std::endl;
    return 0;
}