#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <filesystem>

// Función para calcular el reverso complementario de un k-mer
std::string reverseComplement(const std::string& kmer) {
    std::string rc = kmer;
    std::reverse(rc.begin(), rc.end());
    for (char& base : rc) {
        switch (base) {
            case 'A': base = 'T'; break;
            case 'T': base = 'A'; break;
            case 'C': base = 'G'; break;
            case 'G': base = 'C'; break;
        }
    }
    return rc;
}

int main(int argc, char* argv[]) {
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

    std::string line;
    while (std::getline(genomeFile, line)) {
        if (line.empty() || line[0] == '>') continue;
        for (char base : line) {
            if (base == 'A' || base == 'C' || base == 'G' || base == 'T') {
                sequence += base;
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
        outputFile << "Genoma,kmer,frecuencia,k_value\n";
    }

    // extrae k-mers y escribe en output file
    for (int k : {21, 31}) {
        std::unordered_map<std::string, int> kmer_frequencies;
        
        for (size_t i = 0; i <= sequence.length() - k; ++i) {
            std::string kmer = sequence.substr(i, k);
            std::string rc_kmer = reverseComplement(kmer);
            std::string canonical_kmer = (kmer < rc_kmer) ? kmer : rc_kmer;
            kmer_frequencies[canonical_kmer]++;
        }
        for (const auto& pair : kmer_frequencies) {
            outputFile << genome_filename << "," << pair.first << "," << pair.second << "," << k << "\n";
        }
    }

    outputFile.close();
    std::cout << "Resultados de " << genome_filename << " guardados en " << output_filename << std::endl;

    return 0;
}