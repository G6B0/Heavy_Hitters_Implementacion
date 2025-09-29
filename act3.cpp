#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <filesystem>
#include <cstdlib>
#include <cmath>
#include "TowerSketch.h"
#include "CountSketch.h"
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

// Función para calcular memoria usada por heavy hitters
size_t calculate_hh_memory(const std::unordered_map<uint64_t, int>& kmer_frequencies, 
                          long long frequency_threshold, int k) {
    size_t memory_bytes = 0;
    
    for (const auto& pair : kmer_frequencies) {
        if (pair.second >= frequency_threshold) {
            memory_bytes += sizeof(uint64_t) + sizeof(int) + 24 + k;
        }
    }
    return memory_bytes;
}

std::pair<std::vector<std::pair<uint64_t,bool>>, long long>
leer_subgenoma(const std::string& genome_filename, size_t subset_size, double phi) {
    std::ifstream genomeFile(genome_filename);
    if(!genomeFile.is_open()) {
        std::cerr << "Error al abrir " << genome_filename << std::endl;
        return {};
    }

    std::string sequence;
    size_t bases_read = 0;
    std::string line;

    while(std::getline(genomeFile, line) && bases_read < subset_size) {
        if(line.empty() || line[0]=='>') continue;
        for(char base : line) {
            if(bases_read >= subset_size) break;
            if(base=='A'||base=='C'||base=='G'||base=='T') {
                sequence += base;
                bases_read++;
            }
        }
    }
    genomeFile.close();

    // Extraer k-mers canónicos y contar frecuencias
    std::unordered_map<uint64_t,int> kmer_counts;
    for(int k : {21,31}) {
        if(sequence.size()<k) continue;
        for(size_t i=0; i<=sequence.size()-k; i++){
            uint64_t val = string_to_uint64(sequence.substr(i,k));
            uint64_t rc_val = reverse_complement_uint64(val,k);
            uint64_t canon = std::min(val, rc_val);
            kmer_counts[canon]++;
        }
    }

    // Threshold basado en subset
    long long threshold = static_cast<long long>(phi * sequence.size());

    // Construir vector de ground truth
    std::vector<std::pair<uint64_t,bool>> kmers;
    for(const auto& p : kmer_counts) {
        bool is_HH = (p.second >= threshold);
        kmers.push_back({p.first,is_HH});
    }

    return {kmers, threshold};
}

int main(int argc, char* argv[]) {

    // Tabla de conversión para la codificación y su complemento
    uint8_t base_to_int[128];
    const char int_to_base[4] = {'A', 'C', 'G', 'T'};
    initialize_base_table(base_to_int);
    TowerSketch towersketch(8);
    CountSketch countsketch(15,500000);
    double phi = 0.0000044;
    
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

    int subconjunto_length = static_cast<int>(sequence.length() * 0.3);
    std::cout << "subconjunto de tamaño " << subconjunto_length << " bases" << std::endl;
    auto result = leer_subgenoma(genome_filename, subconjunto_length, phi);
    std::vector<std::pair<uint64_t,bool>> kmers = result.first;
    long long threshold_archivo = result.second;

    std::cout << "threshold = " << threshold_archivo << std::endl;

    // PRE-CALCULAR frecuencias reales MIENTRAS insertamos en los sketches
    std::unordered_map<uint64_t, int> frecuencias_reales;
    
    std::cout << "Extrayendo k-mers y calculando frecuencias..." << std::endl;
    for (int k : {21, 31}) {
        for (size_t i = 0; i <= sequence.length() - k; ++i) {
            uint64_t kmer_val = string_to_uint64(sequence.substr(i, k));
            uint64_t rc_kmer_val = reverse_complement_uint64(kmer_val, k);
            uint64_t canonical_kmer_val = std::min(kmer_val, rc_kmer_val);
            
            // Insertar en sketches
            countsketch.insert(canonical_kmer_val);
            towersketch.insert(canonical_kmer_val);
            
            // Contar frecuencia real
            frecuencias_reales[canonical_kmer_val]++;
        }
    }
    std::cout << "Total de k-mers únicos: " << frecuencias_reales.size() << std::endl;

    // Analisis del ground truth
    int real_HH_count = 0;
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

    // *** CÁLCULO DE ERRORES DE FRECUENCIA ***
    
    double total_abs_error_cs = 0.0;
    double total_rel_error_cs = 0.0;
    int hh_detected_cs = 0;

    double total_abs_error_ts = 0.0;
    double total_rel_error_ts = 0.0;
    int hh_detected_ts = 0;

    double mean_abs_error_cs = -1.0;  // -1 indica "no calculado"
    double mean_rel_error_cs = -1.0;
    double mean_abs_error_ts = -1.0;
    double mean_rel_error_ts = -1.0;

    // Determinar si calcular errores
    bool todos_son_hh = (real_HH_count == (int)kmers.size());
    double hh_ratio = (kmers.size() > 0) ? (double)real_HH_count / kmers.size() : 1.0;
    
    // Calcular errores si:
    // 1. Todos son HH (caso especial, eficiente porque iteramos sobre el mapa)
    // 2. O si hay selectividad (menos del 90% son HH)
    if (todos_son_hh) {
        std::cout << "\n=== CALCULANDO ERRORES (CASO: TODOS SON HH) ===" << std::endl;
        std::cout << "Total k-mers = HH: " << real_HH_count << std::endl;
        
        // Iterar sobre el mapa de frecuencias reales (más eficiente)
        for (const auto& pair : frecuencias_reales) {
            uint64_t kmer_id = pair.first;
            double freq_real = static_cast<double>(pair.second);
            
            // CountSketch
            double freq_est_cs = countsketch.estimar_freq(kmer_id);
            double abs_error_cs = std::abs(freq_est_cs - freq_real);
            double rel_error_cs = (freq_real > 0) ? abs_error_cs / freq_real : 0.0;
            
            total_abs_error_cs += abs_error_cs;
            total_rel_error_cs += rel_error_cs;
            hh_detected_cs++;
            
            // TowerSketch
            double freq_est_ts = towersketch.estimar_freq(kmer_id);
            double abs_error_ts = std::abs(freq_est_ts - freq_real);
            double rel_error_ts = (freq_real > 0) ? abs_error_ts / freq_real : 0.0;
            
            total_abs_error_ts += abs_error_ts;
            total_rel_error_ts += rel_error_ts;
            hh_detected_ts++;
        }
        
        // Calcular promedios
        mean_abs_error_cs = (hh_detected_cs > 0) ? total_abs_error_cs / hh_detected_cs : 0.0;
        mean_rel_error_cs = (hh_detected_cs > 0) ? total_rel_error_cs / hh_detected_cs : 0.0;
        mean_abs_error_ts = (hh_detected_ts > 0) ? total_abs_error_ts / hh_detected_ts : 0.0;
        mean_rel_error_ts = (hh_detected_ts > 0) ? total_rel_error_ts / hh_detected_ts : 0.0;
        
    } else if (hh_ratio < 0.9 && real_HH_count > 0) {
        std::cout << "\n=== CALCULANDO ERRORES (CASO SELECTIVO) ===" << std::endl;
        std::cout << "HH ratio: " << (hh_ratio * 100) << "%" << std::endl;
        std::cout << "Procesando " << real_HH_count << " heavy hitters reales..." << std::endl;

        for (size_t i = 0; i < kmers.size(); i++) {
            if (!kmers[i].second) continue; // Solo procesar heavy hitters reales
            
            // Obtener frecuencia real desde el mapa pre-calculado
            uint64_t kmer_id = kmers[i].first;
            auto it = frecuencias_reales.find(kmer_id);
            if (it == frecuencias_reales.end()) {
                std::cerr << "Error: k-mer no encontrado en frecuencias reales" << std::endl;
                continue;
            }
            double freq_real = static_cast<double>(it->second);
            
            // CountSketch
            double freq_est_cs = countsketch.estimar_freq(kmer_id);
            if (freq_est_cs >= threshold_archivo) { // Solo si fue detectado como HH
                double abs_error = std::abs(freq_est_cs - freq_real);
                double rel_error = (freq_real > 0) ? abs_error / freq_real : 0.0;
                
                total_abs_error_cs += abs_error;
                total_rel_error_cs += rel_error;
                hh_detected_cs++;
            }
            
            // TowerSketch
            double freq_est_ts = towersketch.estimar_freq(kmer_id);
            if (freq_est_ts >= threshold_archivo) { // Solo si fue detectado como HH
                double abs_error = std::abs(freq_est_ts - freq_real);
                double rel_error = (freq_real > 0) ? abs_error / freq_real : 0.0;
                
                total_abs_error_ts += abs_error;
                total_rel_error_ts += rel_error;
                hh_detected_ts++;
            }
        }

        // Calcular promedios
        mean_abs_error_cs = (hh_detected_cs > 0) ? total_abs_error_cs / hh_detected_cs : 0.0;
        mean_rel_error_cs = (hh_detected_cs > 0) ? total_rel_error_cs / hh_detected_cs : 0.0;
        mean_abs_error_ts = (hh_detected_ts > 0) ? total_abs_error_ts / hh_detected_ts : 0.0;
        mean_rel_error_ts = (hh_detected_ts > 0) ? total_rel_error_ts / hh_detected_ts : 0.0;
    } else {
        std::cout << "\n=== SALTANDO CÁLCULO DE ERRORES ===" << std::endl;
        std::cout << "Ratio de HH demasiado alto (" << (hh_ratio * 100) << "%) pero no son todos." << std::endl;
        std::cout << "Cálculo sería muy costoso. Errores marcados como -1 (no calculado)." << std::endl;
    }

    // Calcular métricas de clasificación
    double precision = (TP + FP > 0) ? (double)TP / (TP + FP) : 0.0;
    double recall    = (TP + FN > 0) ? (double)TP / (TP + FN) : 0.0;
    double f1        = (precision + recall > 0) ? 2 * precision * recall / (precision + recall) : 0.0;

    double precision1 = (TP1 + FP1 > 0) ? (double)TP1 / (TP1 + FP1) : 0.0;
    double recall1    = (TP1 + FN1 > 0) ? (double)TP1 / (TP1 + FN1) : 0.0;
    double f11        = (precision1 + recall1 > 0) ? 2 * precision1 * recall1 / (precision1 + recall1) : 0.0;

    double memory = countsketch.get_size_MB();
    double memory1 = towersketch.get_memory_size();

    // Imprimir en consola
    std::cout << "\n=== RESULTADOS FINALES ===" << std::endl;
    std::cout << "CountSketch:" << std::endl;
    std::cout << "Precision: " << precision << std::endl;
    std::cout << "Recall: " << recall << std::endl;
    std::cout << "F1-score: " << f1 << std::endl;
    std::cout << "TP: " << TP << ", FP: " << FP << ", FN: " << FN << ", TN: " << TN << std::endl;
    std::cout << "Memory: " << memory << "MB" << std::endl;
    std::cout << "Mean Absolute Error: " << mean_abs_error_cs << std::endl;
    std::cout << "Mean Relative Error: " << mean_rel_error_cs << std::endl;

    std::cout << "\nTowerSketch:" << std::endl;
    std::cout << "Precision: " << precision1 << std::endl;
    std::cout << "Recall: " << recall1 << std::endl;
    std::cout << "F1-score: " << f11 << std::endl;
    std::cout << "TP: " << TP1 << ", FP: " << FP1 << ", FN: " << FN1 << ", TN: " << TN1 << std::endl;
    std::cout << "Memory: " << memory1 << "MB" << std::endl;
    std::cout << "Mean Absolute Error: " << mean_abs_error_ts << std::endl;
    std::cout << "Mean Relative Error: " << mean_rel_error_ts << std::endl;

    // Escribir en CSV
    std::string csv_filename = "resultados_experimento.csv";
    std::ofstream csvFile;

    bool file_exists = std::filesystem::exists(csv_filename);
    csvFile.open(csv_filename, std::ios::app);

    if (!file_exists) {
        csvFile << "Genoma,Total_kmers,Real_HH,Phi,Threshold,";
        csvFile << "Algoritmo,Precision,Recall,F1_Score,TP,FP,FN,TN,Memory_MB,";
        csvFile << "Mean_Abs_Error,Mean_Rel_Error\n";
    }

    // Extraer nombre del genoma
    std::string genome_name = genome_filename;
    size_t last_slash = genome_name.find_last_of("/\\");
    if (last_slash != std::string::npos) {
        genome_name = genome_name.substr(last_slash + 1);
    }

    // Escribir CountSketch
    csvFile << genome_name << "," << kmers.size() << "," << real_HH_count << ",";
    csvFile << phi << "," << threshold_archivo << ",";
    csvFile << "CountSketch," << precision << "," << recall << "," << f1 << ",";
    csvFile << TP << "," << FP << "," << FN << "," << TN << "," << memory << ",";
    csvFile << mean_abs_error_cs << "," << mean_rel_error_cs << "\n";

    // Escribir TowerSketch
    csvFile << genome_name << "," << kmers.size() << "," << real_HH_count << ",";
    csvFile << phi << "," << threshold_archivo << ",";
    csvFile << "TowerSketch," << precision1 << "," << recall1 << "," << f11 << ",";
    csvFile << TP1 << "," << FP1 << "," << FN1 << "," << TN1 << "," << memory1 << ",";
    csvFile << mean_abs_error_ts << "," << mean_rel_error_ts << "\n";

    csvFile.close();
    std::cout << "\nResultados guardados en: " << csv_filename << std::endl;

    return 0;
}