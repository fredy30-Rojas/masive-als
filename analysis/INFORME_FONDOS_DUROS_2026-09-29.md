# Las dos dianas, medidas con el mismo fondo duro

**29 de septiembre de 2026.** Las dos dianas del proyecto (TDP-43 y SOD1) quedan
puntuadas contra los mismos tres fondos, con el mismo codigo y el mismo criterio.
Esto no existia: el fondo duro de R-BIND 2.0 estaba montado y apagado en SOD1, y en
TDP-43 ya estaba activo pero sin re-verificarse.

## Los numeros

| diana | fondo | n | AUC crudo | AUC/atomo | residual | EF5% | quimias | veredicto |
|---|---|---|---|---|---|---|---|---|
| **TDP-43** | señuelos emparejados | 122 | 0,301 | 0,686 | 0,306 | 0,00 | 2/5 | PASA |
| **TDP-43** | unidores de ARN (duro) | 152 | 0,357 | 0,734 | 0,391 | 0,00 | 2/5 | PASA |
| **TDP-43** | los dos juntos | 274 | 0,332 | 0,713 | 0,356 | 0,00 | 2/5 | PASA |
| **SOD1** | señuelos emparejados | 482 | 0,440 | 0,564 | 0,478 | 3,64 | 5/8 | PASA |
| **SOD1** | unidores de ARN (duro) | 152 | 0,313 | 0,838 | 0,461 | 0,00 | 5/8 | PASA |
| **SOD1** | los dos juntos | 634 | 0,409 | 0,629 | 0,486 | 3,64 | 5/8 | PASA |

## Lo que se ve al ponerlas lado a lado

**1. Las dos se comportan igual, y el motivo es el mismo.** En ambas el **AUC crudo
esta por debajo del azar** (0,357 y 0,313) y el **EF5 es cero en las dos**. El AUC por
atomo sube al pasar del fondo blando al duro (0,686 -> 0,734 en TDP-43; 0,564 -> 0,838
en SOD1). Es decir: contra un unidor de ARN generico, el embudo separa el tamano, y
el tamano es justo lo que un pegamento de ARN tambien tiene.

**2. El PASA de las dos viene entero del emparejado por tamano.** En ninguna de las
dos dianas pasa el corte de cabeza del 5 %: **cero quimias en los dos casos**. TDP-43
solo recupera una quimia si el corte se alarga al 20 %, y alli sigue sin alcanzar las
dos que exige el criterio. El unico criterio que aguanta es el que quita el
confusor de tamano por completo.

**3. TDP-43 es la mas debil de las dos, y el motivo es concreto.** Recupera 2 de 5
quimias frente a 5 de 8 de SOD1. Y el detalle que mas duele: **los tres fragmentos
de RRM2**, que son la unica evidencia estructural directa del sitio de union, caen
en los puestos **152, 155 y 157 de 159**: practicamente al fondo. Son fragmentos
pequenos y el embudo los hunde por su tamano, no por su quimica.

## Que NO se ha hecho aqui

- **TDP-43 no tiene control de redocking posible.** No existe estructura de TDP-43
  con un farmaco pequeno dentro: su bolsillo (4BS2) es una interfaz proteina-proteina
  con un fragmento de Phe74 de otra cadena metido. No hay pose cristalografica que
  redockear, asi que el AUC de TDP-43 no se puede contrastar contra ninguna pose
  conocida. Unico apoyo: que reproduce las energias publicadas (PE859 -7,60 frente
  a -8,49), que NO es lo mismo que ordenar.
- **El tercer fondo (medidos que no unen) sigue sin poder activarse.** Hay 5
  compuestos acoplados en la caja de SOD1, pero estan `pendiente_medicion`: hasta que
  un ensayo los declare negativos no se pueden contar como tales, y eso no se
  resuelve con computo.

## La conclusion, medida y no supuesta

Con las dos dianas puntuadas contra los mismos tres fondos:

**Ordenar una libreria por energia de acoplamiento suelta no sirve para elegir
candidatos.** La cabeza de la lista es grande, no buena. En las dos dianas el AUC
crudo cae por debajo del azar en cuanto el fondo incluye uniones genericas de ARN,
que es el fondo que se parece a lo que un farmaco de verdad tiene delante.

Lo que si sobrevive es el emparejado por tamano. Y eso, medido con el mismo criterio
en las dos dianas, es lo unico que se puede reportar.

## Reproducibilidad

La reejecucion de TDP-43 da un CSV **byte a byte identico** al ya versionado (281
filas, 0 diferencias). Igual que la de SOD1. Los dos resultados son reproducibles,
no heredados.
