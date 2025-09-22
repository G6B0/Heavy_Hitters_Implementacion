#!/bin/bash

# Define la carpeta que contiene los genomas
GENOME_FOLDER="Genomas"

# Itera sobre cada archivo en la carpeta
for filename in "$GENOME_FOLDER"/*; do
    # Verifica si es un archivo regular
    if [ -f "$filename" ]; then
        echo "Ejecutando act3 con: $filename"
        ./act3 "$filename"
    fi
done

echo "Todos los archivos de la carpeta '$GENOME_FOLDER' han sido procesados y los resultados se guardaron en '$OUTPUT_CSV'."