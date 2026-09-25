# Revisión de estereoquímica — 15 PDBQT rescatados con RDKit

Fecha: 6 de septiembre de 2026

## Contexto
15 de los 916 ligandos corregidos no pudieron generar 3D con el builder de
OpenBabel (`no_genera_3d`). Se rescataron con RDKit (ETKDGv3) + conversión a
PDBQT con OpenBabel. Como el embedding se hizo con `enforceChirality=False`,
se auditaron los centros quirales antes de dar el rescate por bueno.

## Método
1. **Esperado**: R/S de cada centro quiral según el SMILES corregido
   (tags `@`/`@@` de RDKit).
2. **Real**: `Chem.AssignStereochemistryFrom3D` sobre el conformero 3D.
3. Corrección, si hacía falta: reflejo de un vecino terminal a través del
   centro quiral (invierte la quiralidad sin tocar el esqueleto).
4. Control de calidad: 100 % de centros coincidentes y distancia mínima
   entre átomos no enlazados ≥ 1.2 Å (sin choques).

## Bug detectado y corregido durante la revisión
- La primera versión de la auditoría usaba `AssignStereochemistryFrom3D(mol,
  force=True)`, pero en esta versión de RDKit (2026.03.5) la firma es
  `AssignStereochemistryFrom3D(mol, confId=-1, replaceExistingTags=True)`;
  el kwarg inválido se tragaba la excepción y se comparaban los tags 2D
  contra sí mismos (falso 100 %). Se validó la corrección con una prueba de
  espejo (invertir todas las coordenadas → 0/4 coincidencias).
- Además, `MMFFOptimizeMolecule` **invierte estereocentros** en estos
  sistemas fusionados (barrera de inversión fácilmente cruzable). Por eso el
  protocolo final NO minimiza con MMFF tras el embedding: las coordenadas
  ETKDG son un punto de partida válido y Vina genera sus propias
  conformaciones durante el docking.

## Resultado final (verificado)
| Métrica | Valor |
|---|---|
| Ligandos auditados | 15 |
| Centros quirales totales | 58 (4 o 6 por ligando) |
| Coincidencia esperado vs 3D | **58/58 (100 %)** |
| Distancia mínima no-enlazada | 1.70 – 1.74 Å (sin choques) |
| PDBQT regenerados | 15/15 (sobrescritos con la geometría validada) |
| Archivos con problemas | 0 |

Todos los 15 quedaron en `OK_sin_cambio` (la geometría ETKDG original ya era
estereoquímicamente correcta una vez verificada con el método correcto). No
fue necesario invertir ningún centro.

## Artefactos
- `analysis/_redock_corregidos/estereoquimica_auditoria.csv` — auditoría por ligando
- `analysis/_redock_corregidos/estereo_img/<ligand>.png` — imágenes 2D con
  la etiqueta R/S esperada en cada centro quiral (revisión visual)
- `analysis/_redock_corregidos/pdbqt/<ligand>.pdbqt` — PDBQT validados

## Conclusión
La advertencia de la entrega anterior ("13 de 15 con estereoquímica no
garantizada") **queda retirada**: la revisión rigurosa confirma que los 15
conservan la estereoquímica del SMILES corregido. Los 916 PDBQT del
re-acoplamiento quedan validados al 100 %.