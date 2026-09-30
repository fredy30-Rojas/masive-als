# El control de redocking de hCA2 — 30 de septiembre de 2026

Uno de los 29 cristales que pasaron la pericia, elegido con tres medidas que no son opinión: la resolución del cristal, el factor R-free, y si el compuesto está entre los 13.051 activos con pchembl.

## El control

| | |
|---|---|
| cristal | **6t4p** |
| farmaco | MHK (naphthalene-1-sulfonamide) |
| peso | 14 atomos pesados, 1 enlaces rotatorios |
| resolucion | 1.75 A |
| R-free | 0.199 |
| union al zinc | N a 2.00 A del Zn |
| agua puente | 1 a menos de 3,5 A |
| ocupacion | 1 copia, 0 conformaciones alternativas, B medio 13 |
| anclaje | ningun polar sin pareja en la proteina |
| encaje | radio 3.6 A en una caja de 24 A |
| potencia | Ki p6.88 |
| SMILES | `NS(=O)(=O)c1cccc2ccccc12` |

## Por qué este y no otro

La pericia dejó 29. Sobre esos 29 se cae la mayoría por resolución: por debajo de 2,5 Å la posición del nitrógeno que toca el zinc tiene medio angstrom de error, y el redocking compara justo ahí. Se cae la otra mayoría porque el compuesto no está entre los activos medidos: el control tiene que hablar de la misma molécula que luego entra en el AUC, o el AUC y el control hablarían de cosas distintas.

De los que pasan, se elige el más rígido y el más pequeño. El más pequeño porque la caja no lo recortará, el más rígido porque entonces la pose del cristal es un mínimo local obvio del paisaje energético, que es lo que un control tiene que ser para que el redocking signifique algo.

## Plan B

Si el redocking de este falla, el motivo importará y estos son los siguientes con el mismo criterio:

| cristal | farmaco | nombre |
|---|---|---|
| 6u4q | PX7 | 4-methyl-1lambda~6~,2,4-benzothiadiazine-1,1 |
| 6ugr | Q71 | 7-fluoro-1lambda~6~,2,4-benzothiadiazine-1,1 |
| 6ugz | Q71 | 7-fluoro-1lambda~6~,2,4-benzothiadiazine-1,1 |
| 6t4o | MJ2 | 3,5-dimethylbenzenesulfonamide |

## Lo que este control NO demuestra

Que el receptor y la caja estén bien. Si el redocking reproduce la pose del cristal por debajo de 2 Å, eso valida la preparación. Si no, el problema está en el receptor (el zinc, el protonado, la caja) y se dice antes de generar un solo señuelo, igual que se hizo con CDK2.

## Ficheros

- `analysis/elegir_control_hca2.py` — criba en dos fases.
- `analysis/fijar_control_hca2.py` — este paso.
- `analysis/calibracion_hca2/pericia_hca2.csv` — los 217 medidos.
- `analysis/calibracion_hca2/control_fijado.json` — el control, para que el redocking no tenga que volver a decidir.
