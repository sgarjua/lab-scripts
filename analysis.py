#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# importaciones etc
from pathlib import Path
import csv, sys

# configuración ===============================================================
TSV = "/data/users/sgarjua/fof_analisis_resultados.tsv"
OUTFILE = Path("/data/users/sgarjua/comparative_table.tsv")
# creamos un diccionario {secuencia: [gos homologia],[gos fantasia]}
resultados = {}
calculos_h = {}
calculos_f = {}

# funciones ===================================================================

def calc_stats(file: Path, species: str):
    # contadores
    id_sin_go = 0
    gos_totales = 0
    protes = 0

    # para cada linea del archivo
    for line in file:
        line = line.strip()
        parts = line.split("\t")
        # tiene que saltarse la primera linea y la cabecera
        if not line or line.startswith("#") or line.startswith("Protein-Accession"):
            continue
        prot = parts[0].strip()
        protes += 1
        # cogemos el id y la metemos en el diccionario si no estaba todavía
        if prot not in resultados:
            resultados[prot] = [[], []]      # [hom, fan]
        # asociamos los términos go a la lista de [gos homologia que están asociado a ese id]
        # En AHRD los GO están en la última columna, coma-separados
        gos_field = parts[-1].strip() if parts else ""
        if "GO:" in gos_field:
            gos = [g.strip() for g in gos_field.split(",") if g.strip()]
            resultados[prot][0] = gos
            for g in gos:
                gos_totales += 1
        else:
            id_sin_go += 1

    # calculos pertinentes:
    # gos totales
    # ids con al menos 1 go
    id_con_go = protes - id_sin_go
    # media de gos/gen
    gos_por_gen = gos_totales / protes

    return protes, gos_totales, id_con_go, id_sin_go, gos_por_gen


# main ========================================================================
def main():
    # compruba que existe el TSV con las especies
    tsv_path = Path(TSV)
    if not tsv_path.exists():
        print(f"[ERROR] No encuentro el TSV: {tsv_path}")
        return

    # abre el TSV y almacena las rutas de los archivos de resultados
    with tsv_path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Espera: especie<TAB>ruta_homologia<TAB>ruta_fantasia
            parts = line.split("\t")
            if len(parts) < 3:
                print(f"[WARN] Línea ignorada (formato esperado: especie<TAB>ruta_homologia<TAB>ruta_fantasia): {line}")
                continue

            species, homologia, fantasia = parts[0].strip(), Path(parts[1].strip()), Path(parts[2].strip())

            if not species or not homologia or not fantasia:
                print(f"[WARN] Falta especie o rutas en la línea: {line}")
                continue

            if not Path(homologia).exists() or Path(homologia).stat().st_size == 0:
                print(f"[WARN] RESULTADOS HOMOLOGÍA no existe o está vacío para {species}")
                continue

            if not Path(fantasia).exists() or Path(fantasia).stat().st_size == 0:
                print(f"[WARN] RESULTADOS FANTASIA no existe o está vacío para {species}")
                continue

            # abrimos los resultados de homología
            with homologia.open(encoding="utf-8") as hom:
                protes_h, gos_totales_h, id_con_go_h, id_sin_go_h, gos_por_gen_h = calc_stats(hom, species)
                calculos_h[species] = [protes_h, gos_totales_h, id_con_go_h, id_sin_go_h, gos_por_gen_h]
                print(calculos_h)

            # abrimos los resultados de fantasia
            with fantasia.open(encoding="utf-8") as fan:
                protes_f, gos_totales_f, id_con_go_f, id_sin_go_f, gos_por_gen_f = calc_stats(fan, species)
                calculos_f[species] = [protes_f, gos_totales_f, id_con_go_f, id_sin_go_f, gos_por_gen_f]
                print(calculos_f)

            cabecera = ["Especie", "Secuencias (H|F)", "Con GO (H|F)", "Sin GO (H|F)", "media GO/sec (H|F)", "GOs solapados"]
            fila = [species, f"{protes_h} | {protes_f}", f"{id_con_go_h} | {id_con_go_f}", f"{id_sin_go_h} | {id_sin_go_h}", f"{gos_por_gen_h} | {gos_por_gen_f}", ]
            with open(OUTFILE, "x") as t:
                t.write(cabecera)
                t.write(fila)
            # calcular el solape

            # hacer el diagrama de venn

    print(resultados)
if __name__ == "__main__":
    main()