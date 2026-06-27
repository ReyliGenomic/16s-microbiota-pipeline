# -*- coding: utf-8 -*-
import os
import glob
import gzip
from Bio.Seq import Seq

# CONFIGURACION SEGUN JTAX (Paper Kim et al. 2025)
# DJ method: 5'-forward reads-3'-NNNNNNNN-3'-reverse complement of reverse reads-5'
SPACER = "NNNNNNNN"
SPACER_QUAL = "IIIIIIII" 

def run_concatenation():
    # 1. DEFINIR CARPETAS
    # IMPORTANTE: Cambia "clean_reads" por la carpeta donde REALMENTE estan tus archivos E1029_1...
    input_dir = "clean_reads" 
    output_dir = "joined_reads"
    
    # Crear carpeta de salida si no existe
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("--> Iniciando modulo de concatenacion tipo JTax (DJ)...")
    
    # 2. BUSCAR ARCHIVOS FORWARD (Terminacion _1.fastq.gz)
    # <--- CAMBIO: Buscamos especificamente los que terminan en _1.fastq.gz
    search_pattern = os.path.join(input_dir, "*_1.fastq.gz")
    files_r1 = glob.glob(search_pattern)
    
    if not files_r1:
        print(f"ERROR: No se encontraron archivos *_1.fastq.gz en '{input_dir}'.")
        return

    for r1_path in files_r1:
        # 3. GENERAR NOMBRES
        # <--- CAMBIO: Reemplazamos _1 por _2 para hallar el reverso
        r2_path = r1_path.replace("_1.fastq.gz", "_2.fastq.gz")
        
        # <--- CAMBIO: Quitamos _1.fastq.gz para obtener el ID limpio (ej: E1029)
        base_name = os.path.basename(r1_path).replace("_1.fastq.gz", "")
        
        output_path = os.path.join(output_dir, base_name + "_joined.fastq.gz")

        if not os.path.exists(r2_path):
            print("ALERTA: Falta el par R2 para " + base_name)
            continue

        print("--> Procesando: " + base_name)

        try:
            with gzip.open(r1_path, "rt") as f1, \
                 gzip.open(r2_path, "rt") as f2, \
                 gzip.open(output_path, "wt") as fout:
                
                while True:
                    # Leer bloques de 4 lineas
                    h1 = f1.readline().strip()
                    s1 = f1.readline().strip()
                    p1 = f1.readline().strip()
                    q1 = f1.readline().strip()
                    
                    h2 = f2.readline().strip()
                    s2 = f2.readline().strip()
                    p2 = f2.readline().strip()
                    q2 = f2.readline().strip()

                    if not h1: break 

                    # --- LOGICA JTAX DJ ---
                    # 1. Forward Sequence (s1)
                    # 2. Spacer (NNNNNNNN)
                    # 3. Reverse Sequence (s2) -> Reverse Complement
                    
                    # Biopython reverse complement
                    s2_rc = str(Seq(s2).reverse_complement())
                    
                    # Invertir la calidad del R2 manualmente
                    q2_rev = q2[::-1]

                    # Concatenacion
                    joined_seq = s1 + SPACER + s2_rc
                    joined_qual = q1 + SPACER_QUAL + q2_rev

                    # Escribir salida
                    fout.write(h1 + " joined_dj\n")
                    fout.write(joined_seq + "\n")
                    fout.write("+\n")
                    fout.write(joined_qual + "\n")
                    
        except Exception as e:
            print("ERROR procesando " + base_name + ": " + str(e))

    print("--> Exito! Archivos generados en 'joined_reads/'")

if __name__ == "__main__":
    run_concatenation()
