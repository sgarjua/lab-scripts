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
import csv, sys
from tempfile import NamedTemporaryFile



# ========================= CONFIGURACIÓN =====================================
TSV = "/data/users/sgarjua/ann_diamond_test/species.tsv"   # Formato: especie<TAB>ruta_fasta
DB1 = "/data/shared_dbs/swissprot/uniprot_sprot_r2025_01.dmnd"
DB2 = "/data/shared_dbs/swissprot/uniprot_trembl_r2025_01.dmnd"
OUTDIR = "/data/users/sgarjua/ann_diamond_test/"

MODE = "blastp"            # "blastp" para proteomas; "blastx" para ADN
THREADS = "24"
EVALUE = "1e-20"
MAX_TARGET_SEQS = "1"
SENSITIVITY = "--sensitive"  # o "--ultra-sensitive"
OUTFMT = "6"      # columnas estándar + título del sujeto

# Post-proceso: mapper GO + AHRD (HAY QUE CONFIGURARLO BIEN TODAVÍA !!!!!!!!)
UNIPROT_GO_MAPPER = "/data/users/sgarjua/ann_diamond_test/uniprot_GO_mapper.py"

AHRD_JAR = "/data/software/AHRD/dist/ahrd.jar"
GO_GAF = "/data/shared_dbs/swissprot/goa_uniprot_all.gaf"
UNIPROT_SPROT = "/data/shared_dbs/swissprot/uniprot_sprot_r2025_01.fasta"
UNIPROT_TREMBL = "/data/shared_dbs/swissprot/uniprot_trembl_r2025_01.fasta"
BLACKLIST = "/data/software/AHRD/test/resources/blacklist_descline.txt"
FILTER_SPROT = "/data/software/AHRD/test/resources/filter_descline_sprot.txt"
FILTER_TREMBL = "/data/software/AHRD/test/resources/filter_descline_trembl.txt"
TOKEN_BLACKLIST = "/data/software/AHRD/test/resources/blacklist_token.txt"
JAVA_XMX = "2g"   # sube a "8g" o más si lo necesitas


# =============================================================================

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

    if out_path.exists():
        print(f"[ALREADY DONE] {species} vs {dbname} → {out_path}")
    else:
        print(f"[RUN ] {species} vs {dbname} → {out_path}")
        try:
            subprocess.run(cmd, check=True)
            print(f"[DONE] {species} vs {dbname}")
        except subprocess.CalledProcessError as e:
            print(f"[FAIL] {species} vs {dbname}. Revisa parámetros/rutas. Detalle: {e}")


def write_ahrd_yaml(tmp_yaml_path: Path,
                    proteins_fasta: Path,
                    go_gaf: Path,
                    out_file: Path,
                    sprot_tsv: Path,
                    trembl_tsv: Path,
                    uniprot_sprot_fa: Path,
                    uniprot_trembl_fa: Path,
                    blacklist: Path,
                    filter_sprot: Path,
                    filter_trembl: Path,
                    token_blacklist: Path):
    """Crea YAML temporal de AHRD (mismos pesos que en el bash)."""
    yaml_text = f"""proteins_fasta: {proteins_fasta}
token_score_bit_score_weight: 0.468
token_score_database_score_weight: 0.2098
token_score_overlap_score_weight: 0.3221
gene_ontology_result: {go_gaf}
reference_go_regex: '^UniProtKB\\t(?<shortAccession>[^\\t]+)\\t[^\\t]+\\t(?!NOT\\|)[^\\t]*\\t(?<goTerm>GO:\\\\d{{7}})'
prefer_reference_with_go_annos: true
output: {out_file}
blast_dbs:
  swissprot:
    weight: 653
    description_score_bit_score_weight: 2.717061
    file: {sprot_tsv}
    database: {uniprot_sprot_fa}
    blacklist: {blacklist}
    filter: {filter_sprot}
    token_blacklist: {token_blacklist}

  trembl:
    weight: 904
    description_score_bit_score_weight: 2.590211
    file: {trembl_tsv}
    database: {uniprot_trembl_fa}
    blacklist: {blacklist}
    filter: {filter_trembl}
    token_blacklist: {token_blacklist}
"""
    tmp_yaml_path.write_text(yaml_text)

def run_ahrd(ahrd_jar: Path, yaml_path: Path, xmx: str = "2g"):
    """java -Xmx<xmx> -jar ahrd.jar tmp.yml"""
    if not ahrd_jar.exists():
        raise FileNotFoundError(f"No se encuentra AHRD JAR: {ahrd_jar}")
    if not yaml_path.exists():
        raise FileNotFoundError(f"No se encuentra YAML de AHRD: {yaml_path}")
    cmd = ["java", f"-Xmx{xmx}", "-jar", str(ahrd_jar), str(yaml_path)]
    print(">> Ejecutando AHRD:", " ".join(cmd))
    subprocess.run(cmd, check=True)

# =============================================================================
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

            # Eliminar los puntos de los archivos de proteinas
            with open(fasta, "r") as infile, open(f"{fasta.stem}.nodots.faa", "w") as outfile:
                for line in infile:
                    line = line.replace(".", "")
                    outfile.write(line)
                fasta = outfile

            # Ejecuta contra las dos bases
            run_diamond(species, fasta, DB1, outdir)
            run_diamond(species, fasta, DB2, outdir)

            # === POST-PROCESO: Extraer IDs → Mapper GO → YAML AHRD → AHRD ===

            # Reconstruimos rutas de salidas DIAMOND que ya generaste:
            sprot_tsv = outdir / f"{species}.{Path(DB1).stem}.o6.txt"
            trembl_tsv = outdir / f"{species}.{Path(DB2).stem}.o6.txt"

            if not sprot_tsv.exists() or not trembl_tsv.exists():
                print(f"[WARN] No encuentro TSVs de DIAMOND para {species}. Saltando post-proceso.")
            else:
                # YAML temporal AHRD + ejecución
                ahrd_out = outdir / f"{species}.proteins.funct_ahrd.tsv"
                with NamedTemporaryFile("w", delete=False, suffix=".yml", dir=str(outdir)) as tmp:
                    yaml_path = Path(tmp.name)

                print("[INFO] Generando YAML temporal para AHRD…")
                write_ahrd_yaml(
                    tmp_yaml_path=yaml_path,
                    proteins_fasta=Path(fasta),  # usa el FASTA que ya estabas usando
                    go_gaf=Path(GO_GAF),
                    out_file=ahrd_out,
                    sprot_tsv=sprot_tsv,
                    trembl_tsv=trembl_tsv,
                    uniprot_sprot_fa=Path(UNIPROT_SPROT),
                    uniprot_trembl_fa=Path(UNIPROT_TREMBL),
                    blacklist=Path(BLACKLIST),
                    filter_sprot=Path(FILTER_SPROT),
                    filter_trembl=Path(FILTER_TREMBL),
                    token_blacklist=Path(TOKEN_BLACKLIST),
                )

                print("[INFO] Ejecutando AHRD…")
                run_ahrd(Path(AHRD_JAR), yaml_path, xmx=JAVA_XMX)
                print(f"[DONE] AHRD → {ahrd_out}")

                # Limpieza del YAML temporal (opcional)
                try:
                    yaml_path.unlink(missing_ok=True)
                except Exception:
                    pass


    print("\nTodo terminado. ✔")

if __name__ == "__main__":
    main()
