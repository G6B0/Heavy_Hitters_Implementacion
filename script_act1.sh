#!/bin/bash

# Define la carpeta que contiene los genomas
GENOME_FOLDER="Genomas"
# Define el nombre del único archivo de salida CSV
OUTPUT_CSV="resultados_totales.csv"

# Limpia el archivo de salida si ya existe para un inicio limpio
if [ -f "$OUTPUT_CSV" ]; then
    rm "$OUTPUT_CSV"
fi

# Itera sobre cada archivo en la carpeta
for filename in "$GENOME_FOLDER"/*; do
    # Verifica si es un archivo regular
    if [ -f "$filename" ]; then
        echo "Procesando el archivo: $filename"
        # Ejecuta tu programa con el nombre del genoma y el archivo de salida
        ./act1 "$filename" "$OUTPUT_CSV"
    fi
done

echo "Todos los archivos de la carpeta '$GENOME_FOLDER' han sido procesados y los resultados se guardaron en '$OUTPUT_CSV'."