#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# importaciones etc
from pathlib import Path
import csv, sys

# configuración ===============================================================
TSV = "/data/users/sgarjua/fof_analisis_resultados.tsv"
# creamos un diccionario {secuencia: [gos homologia],[gos fantasia]}
resultados = {}

# funciones ===================================================================



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

            species, homologia, fantasia = parts[0].strip(), [1].strip()parts, [2].strip()parts

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
                # contadores
                id_sin_go = 0

                # para cada linea del archivo
                for line in hom:
                    line = line.strip()
                    parts = species.split("\t")
                    # tiene que saltarse la primera linea y la cabecera
                    if not line or line.startswith("#") or line.startswith("Protein-Accession"):
                        continue
                    prot = parts[0]
                    # cogemos el id y la metemos en el diccionario si no estaba todavía
                    resultados[prot] = []
                    print(resultados)
                    # asociamos los términos go a la lista de [gos homologia que está asociado a ese id]
                
                # calculos pertinentes:
                # gos totales
                # ids con al menos 1 go
                # media de gos/gen


            # abrimos los resultados de fantasia

                # para cada linea del archivo
                    # cogemos el id y la metemos en el diccionario si no estaba todavía
                    # asociamos los términos go a la lista de [gos fantasia que está asociado a ese id]

                # calculos pertinentes:
                # gos totales
                # ids con al menos 1 go
                # media de gos/gen

                # calcular el solape

                # hacer el diagrama de venn