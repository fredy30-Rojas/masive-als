# Notas del fichero `controles_sod1_v3.csv`

Verdad de referencia de SOD1, versión 3 (20 de septiembre de 2026): los
controles positivos separados por quimiotipo, con la evidencia de cada uno y su
control de redocking. Detalle del trabajo previo en
`AUDITORIA_POSITIVOS_2026-09-20.md` y en
`redocking_trp32/INFORME_CONTROLES_CORREGIDOS_2026-09-20.md`.

## Correcciones del 20 de septiembre de 2026

1. **SMILES de las tres catecolaminas mal escritos.** En la primera versión
   figuraban como `c1ccc(O)c1O`, que es un anillo aromático de cinco miembros:
   RDKit **no puede parsearlos** («Can't kekulize mol»). Corregidos a partir del
   ligando tal y como está en su estructura:

   | ligando | SMILES erróneo | SMILES corregido |
   |---|---|---|
   | isoproterenol | `CC(C)NCC(O)c1ccc(O)c1O` | `CC(C)NCC(O)c1ccc(O)c(O)c1` |
   | dopamina | `NCCc1ccc(O)c1O` | `NCCc1ccc(O)c(O)c1` |
   | adrenalina | `CNC[C@H](O)c1ccc(O)c1O` | `CNCC(O)c1ccc(O)c(O)c1` |

   Se escriben sin estereoquímica asignada porque es como aparecen en las
   estructuras depositadas (el mutante I113T se cristalizó con racematos).

2. **RMSD de redocking recalculados.** Los de la primera versión estaban inflados
   por un error de emparejamiento de átomos en el medidor propio. Ahora:
   isoproterenol 0,48–1,28 Å (pasa), adrenalina 0,69–0,72 Å (pasa), dopamina
   3,13–3,15 Å (no pasa; falla la puntuación, no la búsqueda) y 5-fluorouridina
   11,07–11,29 Å (no pasa, y su pose depositada no es un mínimo del potencial:
   tiene un contacto F···N de 2,07 Å en una estructura de 1,06 Å).

3. **La 5-fluorouridina se queda marcada** como no utilizable como patrón de
   RMSD, por lo del punto 2, aunque su evidencia estructural (co-cristalización
   en 4A7S) sigue siendo válida como indicio de unión.

## Criterio de entrada de un control positivo

Sin cita y sin tipo de ensayo, el compuesto no entra. Las series congénicas se
colapsan a un representante: contar 18 veces la misma química no son 18
positivos (la serie pirazolona aporta uno, CHEMBL2165613).
