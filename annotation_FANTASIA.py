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
    
    if clean_fasta.exists():
        print(f"[ALREADY DONE] El archivo fasta ya estaba limpio")
    else:
        print(f"[RUN] Se va a limpiar el fasta")
        try:
            subprocess.run(cmd, shell=True, check=True)
            print(f"[DONE] Limpieza del fasta completada")
        except subprocess.CalledProcessError as e:
            print(f"[FAIL] Revisa parámetros/rutas. Detalle: {e}")

# ejecutar primer comando
def firt_step(species: str, clean_fasta: str, prefix: str, fantasia_run: str):
    out_path = OUTDIR / species / f"{clean_fasta.stem}_cdhit100.pep"
    print(out_path)

    cmd = [ GENERATE_GPSM,
            "--infile", clean_fasta, 
            "--outpath", fantasia_run, 
            "--prott5",
            "--prefix", prefix,
            "--mode GPU"
    ]

    if out_path.exists():
        print(f"[ALREADY DONE] El primer paso para esta especie ya había sido realizado")
    else:
        print(f"[RUN] Se va a ejecutar el primer paso de FANTASIA")
        try:
            subprocess.run(cmd, check=True)
            print(f"[DONE] Primer paso completado; cmd={cmd}")
        except subprocess.CalledProcessError as e:
            print(f"[FAIL] Revisa parámetros/rutas. Detalle: {e}")

# ejecutar segundo comando
def second_step(prefix: str, fantasia_run: str):
    out_path = fantasia_run / f"{prefix}.log"
    screen_name = f"fantasia_{prefix}"

    cmd = F"{GPU} screen -S {screen_name} -L -Logfile {prefix}.log {LAUNCH_GPSM} -c {fantasia_run} -x {prefix} -m prott5 -o {fantasia_run}"

    if out_path.exists():
        print(f"[ALREADY DONE] El segundo paso para esta especie ya había sido realizado")
    else:
        print(f"[RUN] Se va a ejecutar el segundo paso de FANTASIA")
        try:
            subprocess.run(cmd, shell=True, check=True)
            print(f"[WAIT] Esperando a que el screen '{screen_name}' termine...")
            # Esperar a que el screen termine
            while True:
                result = subprocess.run(f"screen -ls | grep -q {screen_name}", shell=True)
                if result.returncode != 0:
                    print(f"[DONE] Segundo paso completado; cmd={cmd}; Log: {out_path}")
                    break
                time.sleep(10)
        except subprocess.CalledProcessError as e:
            print(f"[FAIL] Revisa parámetros/rutas. Detalle: {e}")




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

            print(f"EJECUTANDO FANTASIA PARA LA ESPECIE {species}")

            # se limpia el fasta
            clean_fasta_name = f"{Path(fasta).stem}.clean.faa"
            clean_fasta = OUTDIR / species / clean_fasta_name
            fasta_cleaner(fasta, clean_fasta)

            # se crea la carpeta fantasia_run dentro de la carpeta de la especie
            fantasia_run = OUTDIR / species / "fantasia_run"
            if not fantasia_run.exists():
                fantasia_run.mkdir(parents=True, exist_ok=True)
                print("[INFO] Carpeta 'fantasia_run' creada correctamente")

            # nos movemos a esa carpeta
            original = Path.cwd()
            os.chdir(fantasia_run)

            # ejecución del primer comando
            firt_step(species, clean_fasta, prefix, fantasia_run)

            # ejecución del segundo comando
            second_step(prefix, fantasia_run)

            os.chdir(original)

    print("\nTodo terminado. ✔")

if __name__ == "__main__":
    main()