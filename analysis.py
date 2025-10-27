#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import csv, sys

from matplotlib_venn import venn2
import matplotlib.pyplot as plt


# configuración ===============================================================
TSV = "/data/users/sgarjua/fof_analisis_resultados.tsv"
OUTFILE = Path("/data/users/sgarjua/comparative_table_final.tsv")

# creamos un diccionario {secuencia: [gos homologia],[gos fantasia]}
resultados = {}

# diccionarios para guardar resultados de los calculos
calculos_h = {}
calculos_f = {}

# funciones ===================================================================

def calc_stats(file, species: str, destino: int):
    """
    Lee un archivo de resultados y guarda los GO en resultados[prot][destino],
    donde destino=0 -> Homología, destino=1 -> FANTASIA.
    Devuelve: (protes, gos_totales, id_con_go, id_sin_go, gos_por_prote)
    """
    id_sin_go = 0
    gos_totales = 0
    protes = 0

    for line in file:
        line = line.strip()
        parts = line.split("\t")
        # saltar cabeceras o vacías
        if not line or line.startswith("#") or line.startswith("Protein-Accession"):
            continue

        prot = parts[0].strip() if parts else ""
        if not prot:
            continue
        protes += 1

        if prot not in resultados:
            resultados[prot] = [[], []]  # [hom, fan]

        gos_field = parts[-1].strip() if parts else ""
        if "GO:" in gos_field:
            gos = [g.strip() for g in gos_field.split(",") if g.strip()]
            resultados[prot][destino] = gos
            gos_totales += len(gos)
        else:
            id_sin_go += 1

    id_con_go = protes - id_sin_go
    gos_por_prote = (gos_totales / protes) if protes else 0.0
    cobertura = (id_con_go / protes) * 100
    return protes, gos_totales, id_con_go, id_sin_go, gos_por_prote, cobertura


def calc_overlap_por_prote(resultados_dict):
    """
    Calcula el solape de GO por proteína.
    Devuelve:
      - overlaps: dict {prote: [GO... que están en Hom y en Fan]}
      - total_solapados: número total de GO en la intersección sumando todas las proteínas
    """
    overlaps = {}
    total_solapados = 0
    for prot, (gos_h, gos_f) in resultados_dict.items():
        sh, sf = set(gos_h), set(gos_f)
        inter = sh & sf
        if inter:
            overlaps[prot] = sorted(inter)
            total_solapados += len(inter)
    return overlaps, total_solapados


def asegurar_cabecera(outfile: Path, header):
    """Escribe cabecera sólo si el fichero no existe aún."""
    modo = "a" if outfile.exists() else "w"
    with outfile.open(modo, newline="", encoding="utf-8") as tsv:
        w = csv.writer(tsv, delimiter="\t")
        if modo == "w":
            w.writerow(header)


def append_fila(outfile: Path, row):
    with outfile.open("a", newline="", encoding="utf-8") as tsv:
        w = csv.writer(tsv, delimiter="\t")
        w.writerow(row)

def calc_total(outfile: Path):
    protes_h, protes_f, id_con_go_h, id_con_go_f, cobertura_h, cobertura_f, gos_por_prote_h, gos_por_prote_f, gos_totales_h, gos_totales_f, total_solapados, solape_h, solape_f = 0
    fila = [
            "\nMEDIA",
            f"{protes_h} | {protes_f}",
            f"{id_con_go_h} | {id_con_go_f}",
            f"{id_sin_go_h} | {id_sin_go_f}",
            f"{cobertura_h:.3f} | {cobertura_f:.3f}",
            f"{gos_por_prote_h:.3f} | {gos_por_prote_f:.3f}",
            f"{gos_totales_h} | {gos_totales_f}",
            total_solapados,
            f"{solape_h:.3f} | {solape_f:.3f}"
            ]
    with outfile.open(encoding="utf-8") as tsv:
        for line in tsv:
            line = line.strip()
            if line.startswith("Especie"):
                continue
            parts = line.split("\t")

            protes = parts[1].split("|")
            protes_h += protes[0]
            protes_f += protes[1]

            id_con_go = parts[2].split("|")
            id_con_go_h += id_con_go[0] / 2
            id_con_go_f += id_con_go[1] / 2

            id_sin_go = parts[3].split("|")
            id_sin_go_h += id_sin_go[0] / 2
            id_sin_go_f += id_sin_go[1] / 2

            cobertura = parts[4].split("|")
            cobertura_h += cobertura[0] / 2
            cobertura_f += cobertura[1] / 2

            gos_por_prote = parts[5].split("|")
            gos_por_prote_h += gos_por_prote[0] / 2
            gos_por_prote_f += gos_por_prote[1] / 2

            gos_totales = parts[6].split("|")
            gos_totales_h += gos_totales[0]
            gos_por_prote_f += gos_totales[1]

            total_solapados += parts[7] / 2

            solape = parts[8].split("|")
            solape_h += solape[0] / 2
            solape_f += solape[1] / 2

    return fila

def diagrama_venn():
    A = 50
    B = 40
    AB = 12

    solo_A = A - AB
    solo_B = B - AB

    venn2(subsets=(solo_A, solo_B, AB), set_labels=("GO-Homología", "GO-Fantasia"))
    plt.title("Venn de los téminos GO anotados con FANTASIA vs anotación por homología, para la lista de especies analizadas")
    plt.show()

    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.savefig("venn.png", dpi=300, bbox_inches="tight")

# main ========================================================================
def main():
    tsv_path = Path(TSV)
    if not tsv_path.exists():
        print(f"[ERROR] No encuentro el TSV: {tsv_path}")
        return

    cabecera = [
        "Especie",
        "Secuencias (H|F)",
        "Con GO (H|F)",
        "Sin GO (H|F)",
        "Cobertura% (H|F)",
        "Media GO/sec (H|F)",
        "GOs totales (H|F)",
        "GOs solapados (total)",
        "% (H|F)"
    ]
    asegurar_cabecera(OUTFILE, cabecera)

    with tsv_path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # Formato esperado: especie<TAB>ruta_homologia<TAB>ruta_fantasia
            parts = line.split("\t")
            if len(parts) < 3:
                print(f"[WARN] Línea ignorada (esperado: especie<TAB>ruta_homologia<TAB>ruta_fantasia): {line}")
                continue

            species = parts[0].strip()
            homologia = Path(parts[1].strip())
            fantasia = Path(parts[2].strip())

            if not species or not homologia or not fantasia:
                print(f"[WARN] Falta especie o rutas en la línea: {line}")
                continue

            if not homologia.exists() or homologia.stat().st_size == 0:
                print(f"[WARN] RESULTADOS HOMOLOGÍA no existe o está vacío para {species}")
                continue

            if not fantasia.exists() or fantasia.stat().st_size == 0:
                print(f"[WARN] RESULTADOS FANTASIA no existe o está vacío para {species}")
                continue

            # IMPORTANTE: limpiar resultados por especie
            resultados.clear()

            # homología
            with homologia.open(encoding="utf-8") as hom:
                protes_h, gos_totales_h, id_con_go_h, id_sin_go_h, gos_por_prote_h, cobertura_h = calc_stats(hom, species, destino=0)
                calculos_h[species] = [protes_h, gos_totales_h, id_con_go_h, id_sin_go_h, gos_por_prote_h, cobertura_h]

            # FANTASIA
            with fantasia.open(encoding="utf-8") as fan:
                protes_f, gos_totales_f, id_con_go_f, id_sin_go_f, gos_por_prote_f, cobertura_f = calc_stats(fan, species, destino=1)
                calculos_f[species] = [protes_f, gos_totales_f, id_con_go_f, id_sin_go_f, gos_por_prote_f, cobertura_f]

            # calcular el solape de GO por proteína
            overlaps, total_solapados = calc_overlap_por_prote(resultados)
            solape_h = (total_solapados/gos_totales_h)*100
            solape_f = (total_solapados/gos_totales_f)*100

            # construir fila
            fila = [
                species,
                f"{protes_h} | {protes_f}",
                f"{id_con_go_h} | {id_con_go_f}",
                f"{id_sin_go_h} | {id_sin_go_f}",
                f"{cobertura_h:.3f} | {cobertura_f:.3f}",
                f"{gos_por_prote_h:.3f} | {gos_por_prote_f:.3f}",
                f"{gos_totales_h} | {gos_totales_f}",
                total_solapados,
                f"{solape_h:.3f} | {solape_f:.3f}"
            ]
            append_fila(OUTFILE, fila)

        calc_total(OUTFILE)        
        append_fila(OUTFILE, total)


if __name__ == "__main__":
    main()
