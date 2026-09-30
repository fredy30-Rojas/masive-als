# Control de redocking de hCA2: 6T4P

La naftalensulfonamida del cristal, acoplada de nuevo en la misma caja, con el zinc dentro del receptor.

## Resultado

| variante | RMSD de bolsillo | piso | alineado | N-Zn pose | N-Zn cristal | veredicto |
|---|---|---|---|---|---|---|
| caja 20 A, depth 20 | **0.27 A** | 0.27 A | 0.06 A | 1.92 A | 2.00 A | PASA |
| caja 20 A, depth 32 | **0.26 A** | 0.26 A | 0.06 A | 1.92 A | 2.00 A | PASA |
| caja 22 A, depth 20 | **0.23 A** | 0.23 A | 0.03 A | 1.96 A | 2.00 A | PASA |
| caja 22 A, depth 32 | **0.23 A** | 0.23 A | 0.06 A | 1.92 A | 2.00 A | PASA |
| caja 24 A, depth 20 | **0.27 A** | 0.27 A | 0.06 A | 1.92 A | 2.00 A | PASA |
| caja 24 A, depth 32 | **0.26 A** | 0.26 A | 0.06 A | 1.92 A | 2.00 A | PASA |

Liston de bolsillo: 2.0 A. Ademas, para hCA2 hace falta que la pose acoplada se coordine al zinc a menos de 2.6 A: sin eso, un buen RMSD puede ser una pose en el sitio correcto pero sin tocar el metal, que es justo lo que el AUC de hCA2 no debe premiar.

## Que significa

El control pasa. El receptor con su zinc y la caja de 24 A reproducen la pose del cristal, asi que el AUC que venga midera el embudo y no una falta de preparacion. Este es el primer control de redocking del proyecto que sirve para algo, despues del de `andr`.

## Ficheros

- `analysis/preparar_receptor_hca2.py` — receptor con el zinc.
- `analysis/calibracion_hca2/control_redocking_hca2.py` — este redocking.
- `analysis/calibracion_hca2/barrido_redocking_hca2.txt` — todas las variantes.
