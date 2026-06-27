# -*- coding: utf-8 -*-
import os
import glob
import gzip

class Clear_String_Optimized:
    def __init__(self):
        # 1. Configuración de carpetas
        # Ajusta esto si tus archivos están en otra subcarpeta
        self.input_dir = "joined_reads"           
        self.output_dir = "joined_reads_clean"    
        
        # 2. Cadenas a eliminar (según el paper)
        self.seq_remove = "NNNNNNNN"  # El espaciador de secuencia
        self.qual_remove = "IIIIIIII" # El espaciador de calidad

        # Crear carpeta de salida
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        # 3. Ejecutar limpieza
        self.process_files()

    def process_files(self):
        # Busca archivos comprimidos .fastq.gz
        search_path = os.path.join(self.input_dir, "*_joined.fastq.gz")
        files = glob.glob(search_path)
        
        if not files:
            print(f"ERROR: No encontré archivos .fastq.gz en '{self.input_dir}'")
            return

        print(f"--> Encontrados {len(files)} archivos. Iniciando limpieza DOBLE (Ns e Is)...")

        for filepath in files:
            filename = os.path.basename(filepath)
            output_path = os.path.join(self.output_dir, filename)
            
            try:
                # Lectura y escritura streaming (sin llenar la RAM)
                with gzip.open(filepath, "rt") as fin, \
                     gzip.open(output_path, "wt") as fout:
                    
                    while True:
                        header = fin.readline()
                        if not header: break 
                        
                        seq = fin.readline()
                        plus = fin.readline()
                        qual = fin.readline()

                        # --- ELIMINAR AMBOS ESPACIADORES SIMULTÁNEAMENTE ---
                        # Esto garantiza que el largo de secuencia y calidad sigan coincidiendo
                        clean_seq = seq.strip().replace(self.seq_remove, "")
                        clean_qual = qual.strip().replace(self.qual_remove, "")
                        
                        # Escribir
                        fout.write(header)
                        fout.write(clean_seq + "\n")
                        fout.write(plus)
                        fout.write(clean_qual + "\n")

                print(f"Completado y validado: {filename}")

            except Exception as e:
                print(f"Error procesando {filename}: {e}")

        print(f"--> ¡Listo! Archivos limpios guardados en '{self.output_dir}/'")

if __name__ == "__main__":
    Clear_String_Optimized()
