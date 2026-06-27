import os
import argparse

def create_manifest_single_end(fastq_directory):
    """
    Genera un archivo manifiesto de QIIME 2 para secuencias Single End (o concatenadas).
    Adapta las rutas para que funcionen dentro de Docker (/data/...).
    """
    
    # 1. Verificar directorio
    if not os.path.isdir(fastq_directory):
        print(f"Error: El directorio '{fastq_directory}' no existe.")
        return

    output_file = 'manifest.tsv'
    
    # Obtener la ruta absoluta de la raíz del proyecto (donde corres el script)
    project_root = os.getcwd()
    
    # Obtener lista de archivos y ordenarlos
    files = sorted([f for f in os.listdir(fastq_directory) if f.endswith(".fastq.gz")])
    
    print(f"Procesando {len(files)} archivos en: {fastq_directory}\n")

    with open(output_file, 'w') as f:
        # --- CAMBIO 1: Encabezado para Single End ---
        # Ya no pedimos R1 y R2, solo "absolute-filepath"
        f.write('sample-id\tabsolute-filepath\n')

        for filename in files:
            try:
                # --- CAMBIO 2: Extraer ID ---
                # Tomamos todo antes del primer guion bajo '_'
                sample_id = filename.split('_')[0]
                
                # --- CAMBIO 3: Construir ruta para DOCKER ---
                # 1. Ruta real en el servidor
                host_path = os.path.join(fastq_directory, filename)
                
                # 2. Ruta relativa desde donde estás parado
                # (ej: joined_reads_clean/1001_joined.fastq.gz)
                rel_path = os.path.relpath(host_path, project_root)
                
                # 3. Ruta final para Docker (/data/...)
                docker_path = os.path.join("/data", rel_path)

                # Escribir línea
                f.write(f"{sample_id}\t{docker_path}\n")
                print(f"  + Muestra: {sample_id} -> {docker_path}")

            except IndexError:
                print(f"  - Advertencia: El archivo '{filename}' no tiene formato correcto.")

    print(f"\n¡Listo! Se ha creado '{output_file}' exitosamente.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generar manifest para archivos concatenados (Single End) en Docker.")
    parser.add_argument(
        'directorio_fastq',
        type=str,
        help='Nombre de la carpeta donde están tus archivos concatenados (ej: joined_reads_clean)'
    )
    args = parser.parse_args()
    
    create_manifest_single_end(args.directorio_fastq)
