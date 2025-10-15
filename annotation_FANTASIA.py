#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# nota: previamente a la ejecución de este script tienes que ejecutar los siguientes comandos:
# conda activate gopredsim
# export LD_LIBRARY_PATH=/data/users/sgarjua/00_software/miniconda3/envs/gopredsim/lib/

# 

# importaciones etc
import subprocess
from pathlib import Path
import csv, sys
import os


# configuración ===============================================================
TSV = "/data/users/sgarjua/ann_diamond_test/species.tsv" # archivo con: nombre especie(=nombre carpeta) + ruta al fasta
OUTDIR = Path("/data/users/sgarjua/SofiaFantasia/") # carpeta con las carpetas de las especies con los fastas
GENERATE_GPSM = "/data/users/sgarjua/00_software/FANTASIA/generate_gopredsim_input_files.sh"
LAUNCH_GPSM = "/data/users/sgarjua/00_software/FANTASIA/launch_gopredsim_pipeline.sh"
GPU = "CUDA_VISIBLE_DEVICES=1"

# funciones ===================================================================
# limpiar el fasta
def fasta_cleaner(fasta: Path, clean_fasta: Path):

    cmd = f"sed -r 's/ .+//' {fasta} > {clean_fasta}"
    
    #EJECUTAAAAAR, pero comprobar primero si ya existe


# ejecutar primer comando
def firt_step(clean_fasta: Path, prefix: str, fantasia_run: Path):

    cmd = [ GENERATE_GPSM,
            "--infile", clean_fasta, 
            "--outpath", fantasia_run, 
            "--prott5",
            "--prefix", prefix,
            "--mode GPU"
    ]

    print(cmd)

# ejecutar segundo comando
def second_step(prefix: str, fantasia_run: Path):

    cmd = [ GPU,
            F"screen -L -Logfile {prefix}.log",
            LAUNCH_GPSM,
            "-c", fantasia_run, 
            "-x", prefix, 
            "-m prott5",
            "-o", fantasia_run
    ]

    print(cmd)

# main ========================================================================
def main():
    # compruba que existe el TSV con las especies
    tsv_path = Path(TSV)
    if not tsv_path.exists():
        print(f"[ERROR] No encuentro el TSV: {tsv_path}")
        return

    # abre el TSV y almacena esocies y fasta
    with tsv_path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Espera: especie<TAB>ruta_fasta
            parts = line.split("\t")
            if len(parts) < 2:
                print(f"[WARN] Línea ignorada (formato esperado: especie<TAB>fasta): {line}")
                continue

            species, fasta = parts[0].strip(), parts[1].strip()

            if not species or not fasta:
                print(f"[WARN] Falta especie o fasta en la línea: {line}")
                continue

            if not Path(fasta).exists() or Path(fasta).stat().st_size == 0:
                print(f"[WARN] FASTA no existe o está vacío para {species}: {fasta}")
                continue

            # creamos el prefijo
            parts = species.split("_")
            prefix = parts[0][:2] + parts[1][:3]
            if len(parts) > 2:
                prefix += parts[2][:3]
            
            # se limpia el fasta
            clean_fasta_name = f"{Path(fasta).stem}.clean.faa"
            clean_fasta = OUTDIR / species / clean_fasta
            fasta_cleaner(fasta, species, clean_fasta)

            # se crea la carpeta run_fantasia dentro de la carpeta de la especie
            run_fantasia = OUTDIR / species / run_fantasia
            # if not run_fantasia.exists():
            #     out.mkdir(parents=True, exist_ok=True)
            #     print("[INFO] Carpeta 'run_fantasia' creada correctamente")

            # nos movemos a esa carpeta
            original = Path.cwd()
            os.chdir(run_fantasia)

            # ejecución del primer comando

            # ejecución del segundo comando


            os.chdir(original)

    print("\nTodo terminado. ✔")

if __name__ == "__main__":
    main()