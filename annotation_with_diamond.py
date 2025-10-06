#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Uso:
1) Edita la sección CONFIG de abajo (TSV, DB1, DB2, OUTDIR, etc.).
2) Ejecuta:  python annotation_with_diamond.py
3) Los resultados se guardan en OUTDIR con nombres: especie.dbname.o6.txt
   (dbname es el nombre del archivo .dmnd sin la ruta)
"""

import subprocess
from pathlib import Path

# =========================
# CONFIG (edita aquí)
# =========================
TSV = "/data/users/sgarjua/ann_diamond_test/species.tsv"   # Formato: especie<TAB>ruta_fasta
DB1 = "/data/shared_dbs/swissprot/uniprot_sprot_r2025_01.dmnd"
DB2 = "/data/shared_dbs/swissprot/uniprot_trembl_r2025_01.dmnd"
OUTDIR = "/data/users/sgarjua/ann_diamond_test/"

MODE = "blastp"            # "blastp" para proteomas; "blastx" para ADN
THREADS = 24
EVALUE = "1e-20"
MAX_TARGET_SEQS = 1
SENSITIVITY = "--sensitive"  # o "--ultra-sensitive"
OUTFMT = 6      # columnas estándar + título del sujeto

# =========================
# NO TOCAR DE AQUÍ ABAJO :)
# =========================

def run_diamond(species: str, fasta: str, db: str, outdir: Path):
    dbname = Path(db).stem  # nombre corto de la DB para el fichero de salida
    out_path = outdir / f"{species}.{dbname}.o6.txt"

    cmd = [
        "diamond", MODE,
        "--query", fasta,
        "--db", db,
        "--outfmt", OUTFMT,
        "--max-target-seqs", str(MAX_TARGET_SEQS),
        "--evalue", EVALUE,
        "--out", str(out_path),
        "--threads", str(THREADS),
        SENSITIVITY,
    ]

    print(f"[RUN ] {species} vs {dbname} → {out_path}")
    try:
        subprocess.run(cmd, check=True)
        print(f"[DONE] {species} vs {dbname}")
    except subprocess.CalledProcessError as e:
        print(f"[FAIL] {species} vs {dbname}. Revisa parámetros/rutas. Detalle: {e}")

def main():
    outdir = Path(OUTDIR)
    outdir.mkdir(parents=True, exist_ok=True)

    tsv_path = Path(TSV)
    if not tsv_path.exists():
        print(f"[ERROR] No encuentro el TSV: {tsv_path}")
        return

    if not Path(DB1).exists() or not Path(DB2).exists():
        print(f"[ERROR] Revisa que existan DB1 y DB2:\n  DB1={DB1}\n  DB2={DB2}")
        return

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

            # Ejecuta contra las dos bases
            run_diamond(species, fasta, DB1, outdir)
            run_diamond(species, fasta, DB2, outdir)

    print("\nTodo terminado. ✔")

if __name__ == "__main__":
    main()
