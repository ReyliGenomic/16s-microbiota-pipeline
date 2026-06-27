# 1. Crear carpetas de salida (si no existen)
mkdir -p clean_reads reports

# 2. Bucle para procesar cada archivo R1
for R1 in *_R1_001.fastq.gz; do
    # Verificar si existen archivos
    [ -e "$R1" ] || continue

    # Construir nombres de archivos
    R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    BASE_NAME="${R1/_R1_001.fastq.gz/}"

    # Verificar que el par R2 existe
    if [ -f "$R2" ]; then
        echo "--> Procesando muestra: $BASE_NAME"
        
        # EJECUTAR FASTP CON CORRECCIÓN:
        # --trim_front1 17: Corta el primer Forward (17 pb) del inicio de R1
        # --trim_front2 21: Corta el primer Reverse (21 pb) del inicio de R2
        # --detect_adapter_for_pe: Sigue buscando y eliminando adaptadores Illumina al final
        
        fastp -i "$R1" -I "$R2" \
              -o "clean_reads/${BASE_NAME}_R1_clean.fastq.gz" \
              -O "clean_reads/${BASE_NAME}_R2_clean.fastq.gz" \
              --trim_front1 17 \
              --trim_front2 21 \
              --detect_adapter_for_pe \
              -h "reports/${BASE_NAME}_fastp_report.html" \
              -j "reports/${BASE_NAME}_fastp_report.json" \
              -w 4
              
        echo "    -> Completado. Guardado en clean_reads/"
    else
        echo "ALERTA: No se encontró pareja R2 para $R1"
    fi
done

echo "--> ¡Proceso finalizado con éxito!"# 1. Crear carpetas de salida (si no existen)
mkdir -p clean_reads reports

# 2. Bucle para procesar cada archivo R1
for R1 in *_R1_001.fastq.gz; do
    # Verificar si existen archivos
    [ -e "$R1" ] || continue

    # Construir nombres de archivos
    R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    BASE_NAME="${R1/_R1_001.fastq.gz/}"

    # Verificar que el par R2 existe
    if [ -f "$R2" ]; then
        echo "--> Procesando muestra: $BASE_NAME"
        
        # EJECUTAR FASTP CON CORRECCIÓN:
        # --trim_front1 17: Corta el primer Forward (17 pb) del inicio de R1
        # --trim_front2 21: Corta el primer Reverse (21 pb) del inicio de R2
        # --detect_adapter_for_pe: Sigue buscando y eliminando adaptadores Illumina al final
        
        fastp -i "$R1" -I "$R2" \
              -o "clean_reads/${BASE_NAME}_R1_clean.fastq.gz" \
              -O "clean_reads/${BASE_NAME}_R2_clean.fastq.gz" \
              --trim_front1 17 \
              --trim_front2 21 \
              --detect_adapter_for_pe \
              -h "reports/${BASE_NAME}_fastp_report.html" \
              -j "reports/${BASE_NAME}_fastp_report.json" \
              -w 4
              
        echo "    -> Completado. Guardado en clean_reads/"
    else
        echo "ALERTA: No se encontró pareja R2 para $R1"
    fi
done

echo "--> ¡Proceso finalizado con éxito!"
