# Por que ranking_afinidades.csv no coincide con resultados_libreria_total.csv

**24 de septiembre de 2026.** Este ranking se limpio: se retiraron **41 ligandos**
(120 filas) que tenian la afinidad FALSA.

## El porque

La auditoria de tipos de atomo del 23 sep 2026 (`_auditoria_tipos_cribado.log`)
encontro en la libreria ligandos con un elemento critico (boro, silicio, selenio,
telurio...) **nombrado** en el PDBQT pero **tipado** como otro elemento (casi
siempre C). Vina tipa la malla por el tipo, no por el nombre: esos ligandos se
acoplaron y puntuaron fingiendo ser su analogo de carbono, asi que su numero del
ranking no mide su quimica. Es el mismo defecto que dejo a `DEC_CHEMBL4543460`
sin pose en la validacion de TDP-43 (boron no tipable), pero en su variante mas
traicionera: aqui el ligando SI se acopla, y puntua mal sin avisar.

## Que se hizo

- Los 41 ligandos se retiraron de `ranking_afinidades.csv` con
  `_limpiar_ranking_tipos.py`, que aplica la MISMA logica de
  `auditar_tipos_cribado.py` (importada, no copiada).
- Las filas retiradas NO se borraron: estan en
  `ranking_afinidades_descartes_tipos.csv` por si hay que auditar.
- No se tocó `resultados_libreria_total.csv` ni el top100 ni las listas de
  candidatos: se verifico que NINGUNO de estos ligandos llego a lista de
  candidatos (cero contaminacion; el ranking bruto era el unico sitio sucio).

## Impacto por diana

- FUS: 38 filas retiradas
- SOD1: 41 filas retiradas
- TDP43_v2: 41 filas retiradas

## Lo que NO dice esta nota

Que los ligandos retirados sean malos compuestos: varios son quimicamente
interesantes (boro, silicio). Solo que SU AFINIDAD VINA aqui es un artefacto.
Para puntuarlos de verdad habria que re-prepararlos y re-acoplarlos con un
motor que entienda esos elementos, y eso es otra conversacion.
