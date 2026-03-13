#!/bin/bash
# Script para deduplicar FASTQ paired-end por separado con cd-hit-dup
# Entrada: /Volumes/Viroma/cartel/espanoles/*.fastq
# Salida:  /Volumes/Viroma/cartel/espanoles/deduplicados/

input_dir="/Volumes/Viroma/cartel/espanoles"
output_dir="/Volumes/Viroma/cartel/espanoles/deduplicados"

mkdir -p "$output_dir"

for r1 in ${input_dir}/*_1.fastq; do
    sample=$(basename "$r1" _1.fastq)
    r2="${input_dir}/${sample}_2.fastq"

    echo "=== Procesando muestra: $sample ==="

    dedup_r1="${output_dir}/dedup_${sample}_1.fastq"
    dedup_r2="${output_dir}/dedup_${sample}_2.fastq"

    # Deduplicar R1
    cd-hit-dup -i "$r1" -o "$dedup_r1" -u 100 -T 8 -M 0

    # Deduplicar R2
    cd-hit-dup -i "$r2" -o "$dedup_r2" -u 100 -T 8 -M 0

    # Contar número de lecturas (cada 4 líneas = 1 read)
    reads_r1=$(( $(wc -l < "$dedup_r1") / 4 ))
    reads_r2=$(( $(wc -l < "$dedup_r2") / 4 ))

    echo "Reads en ${sample}_1: $reads_r1"
    echo "Reads en ${sample}_2: $reads_r2"

    # Verificar si son iguales
    if [ "$reads_r1" -eq "$reads_r2" ]; then
        echo "✅ Pares balanceados: $sample"
    else
        echo "⚠️  Desbalance detectado en $sample"
    fi
    echo
done
