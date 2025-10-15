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


# configuración ===============================================================
TSV = "/data/users/sgarjua/ann_diamond_test/species.tsv" # archivo con: nombre especie(=nombre carpeta) + ruta al fasta
OUTDIR = "/data/users/sgarjua/SofiaFantasia/" # carpeta con las carpetas de las especies con los fastas


# funciones ===================================================================
# crear carpeta

# ejecutar primer comando

# ejecutar segundo comando


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
            prefix = parts[0][:1] + parts[1][:2]
            print(prefix)
            
            # se limpia el fasta

            # se crea la carpeta run_fantasia dentro de la carpeta de la especie

            # nos movemos a esa carpeta

            # ejecución del primer comando

            # ejecución del segundo comando

    print("\nTodo terminado. ✔")

if __name__ == "__main__":
    main()