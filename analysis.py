#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import csv
import sys
import re
from typing import Dict, List, Set, Tuple

# ================================ configuración ==============================
TSV = "/data/users/sgarjua/fof_analisis_resultados.tsv"
OUTFILE = Path("/data/users/sgarjua/comparative_table.tsv")


# ================================ utilidades =================================
GO_REGEX = re.compile(r"GO:\d{7}")

def extract_gos_from_line(parts: List[str]) -> Set[str]:
    """
    Extrae términos GO de una línea ya dividida en columnas.
    Intenta ser robusto: busca tokens GO:####### en cualquier columna.
    Devuelve un set sin duplicados.
    """
    gos: Set[str] = set()
    for field in parts:
        # separar por comas/espacios/punto y coma
        tokens = re.split(r"[,\s;]+", field.strip())
        for tk in tokens:
            if GO_REGEX.fullmatch(tk):
                gos.add(tk)
    return gos


def parse_results(file_path: Path) -> Dict[str, Set[str]]:
    """
    Devuelve un diccionario {protein_id: set(GO)} para un archivo de resultados.
    - Ignora líneas vacías, comentarios y cabeceras típicas.
    """
    mapping: Dict[str, Set[str]] = {}
    with file_path.open(encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("Protein-Accession"):
                continue
            parts = line.split("\t")
            if not parts:
                continue
            prot = parts[0].strip()
            if not prot:
                continue
            gos = extract_gos_from_line(parts)
            mapping[prot] = gos  # set vacío si no hay GO -> también lo guardamos
    return mapping


def basic_stats(mapping: Dict[str, Set[str]]) -> Tuple[int, int, int, int, float]:
    """
    Calcula métricas básicas para un mapeo {prot: set(GO)}:
    - protes (nº IDs)
    - gos_totales (suma de longitudes de set)
    - id_con_go (con al menos 1 GO)
    - id_sin_go (con 0 GO)
    - gos_por_gen (media de GO por proteína; 0.0 si no hay proteínas)
    """
    protes = len(mapping)
    gos_totales = sum(len(s) for s in mapping.values())
    id_con_go = sum(1 for s in mapping.values() if len(s) > 0)
    id_sin_go = protes - id_con_go
    gos_por_gen = (gos_totales / protes) if protes > 0 else 0.0
    return protes, gos_totales, id_con_go, id_sin_go, gos_por_gen


def overlap_stats(h_map: Dict[str, Set[str]], f_map: Dict[str, Set[str]]):
    """
    Calcula métricas de solape entre homología (h_map) y FANTASIA (f_map).
    Devuelve dict con:
      - prots_both: proteínas presentes en ambos mapeos (con independencia de tener GO)
      - prots_solo_h / prots_solo_f
      - go_overlap: suma de |∩| sobre proteínas presentes en al menos uno
      - go_union: suma de |∪|
      - mean_jaccard: media de Jaccard sobre proteínas con |∪|>0
      - go_unicos_h: suma de GO que están solo en homología (sobre todas las proteínas)
      - go_unicos_f: idem solo en FANTASIA
    """
    prots_h = set(h_map.keys())
    prots_f = set(f_map.keys())
    all_prots = prots_h | prots_f
    both = prots_h & prots_f

    go_overlap_sum = 0
    go_union_sum = 0
    jaccs: List[float] = []
    unique_h_sum = 0
    unique_f_sum = 0

    for p in all_prots:
        h_set = h_map.get(p, set())
        f_set = f_map.get(p, set())
        inter = h_set & f_set
        union = h_set | f_set
        go_overlap_sum += len(inter)
        go_union_sum += len(union)
        unique_h_sum += len(h_set - f_set)
        unique_f_sum += len(f_set - h_set)
        if len(union) > 0:
            jaccs.append(len(inter) / len(union))

    mean_jacc = sum(jaccs) / len(jaccs) if jaccs else 0.0

    return {
        "prots_both": len(both),
        "prots_solo_h": len(prots_h - prots_f),
        "prots_solo_f": len(prots_f - prots_h),
        "go_overlap": go_overlap_sum,
        "go_union": go_union_sum,
        "mean_jaccard": mean_jacc,
        "go_unicos_h": unique_h_sum,
        "go_unicos_f": unique_f_sum,
    }


# =================================== main ====================================
def main():
    tsv_path = Path(TSV)
    if not tsv_path.exists():
        print(f"[ERROR] No encuentro el TSV: {tsv_path}", file=sys.stderr)
        sys.exit(1)

    # Preparamos el OUTFILE
    OUTFILE.parent.mkdir(parents=True, exist_ok=True)
    outfh = OUTFILE.open("w", newline="", encoding="utf-8")
    writer = csv.writer(outfh, delimiter="\t")

    # Cabecera de la tabla comparativa
    writer.writerow([
        "species",
        # homología
        "protes_h", "gos_totales_h", "id_con_go_h", "id_sin_go_h", "gos_por_gen_h",
        # FANTASIA
        "protes_f", "gos_totales_f", "id_con_go_f", "id_sin_go_f", "gos_por_gen_f",
        # solape
        "prots_both", "prots_solo_h", "prots_solo_f",
        "go_overlap", "go_union", "mean_jaccard",
        "go_unicos_h", "go_unicos_f"
    ])

    # Leemos el fichero maestro especie \t ruta_h \t ruta_f
    with tsv_path.open(encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                print(f"[WARN] Línea ignorada (se esperan 3 columnas): {line}", file=sys.stderr)
                continue

            species = parts[0].strip()
            homologia = Path(parts[1].strip())
            fantasia = Path(parts[2].strip())

            if not species:
                print(f"[WARN] Especie vacía: {line}", file=sys.stderr)
                continue

            if not homologia.exists() or homologia.stat().st_size == 0:
                print(f"[WARN] HOMOLOGÍA no existe o está vacío para {species}: {homologia}", file=sys.stderr)
                continue

            if not fantasia.exists() or fantasia.stat().st_size == 0:
                print(f"[WARN] FANTASIA no existe o está vacío para {species}: {fantasia}", file=sys.stderr)
                continue

            # ---- Parseo por especie (diccionarios independientes por especie)
            h_map = parse_results(homologia)
            f_map = parse_results(fantasia)

            # ---- Estadísticos básicos
            protes_h, gos_tot_h, con_h, sin_h, mean_h = basic_stats(h_map)
            protes_f, gos_tot_f, con_f, sin_f, mean_f = basic_stats(f_map)

            # ---- Solape
            ov = overlap_stats(h_map, f_map)

            # ---- Escribir fila
            writer.writerow([
                species,
                protes_h, gos_tot_h, con_h, sin_h, f"{mean_h:.4f}",
                protes_f, gos_tot_f, con_f, sin_f, f"{mean_f:.4f}",
                ov["prots_both"], ov["prots_solo_h"], ov["prots_solo_f"],
                ov["go_overlap"], ov["go_union"], f"{ov['mean_jaccard']:.4f}",
                ov["go_unicos_h"], ov["go_unicos_f"],
            ])

            # Log opcional a stdout para seguimiento
            print(f"[OK] {species} -> H({con_h}/{protes_h}), F({con_f}/{protes_f}), "
                  f"Jaccard medio={ov['mean_jaccard']:.4f}")

    outfh.close()
    print(f"[DONE] Tabla comparativa escrita en: {OUTFILE}")


if __name__ == "__main__":
    main()
