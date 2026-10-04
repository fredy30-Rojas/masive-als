INFORME DE AMPLIACION DE LOS POSITIVOS DE hCA2 — 2026-10-04

AMPLIACION DE LOS POSITIVOS DE hCA2 — 2026-10-04 13:20
   objetivo: 300 positivos

SERIE MEDIDA DE hCA2
   ya descargada: 1200 activos (serie_hca2_auditoria.csv)
   ampliada en disco: 7200 activos
   serie de trabajo: 7200 activos medidos

POSITIVOS
   se conservan los 60 de siempre (lo ya acoplado sigue siendo comparable)
   activos con mol legible: 7200
   familias de Murcko de 5+ activos: 318
   familias distintas en total: 2013
   familias ya cubiertas por los 60 positivos: 59

PREPARANDO LOS NUEVOS (receta canonica del proyecto)
   preparados 229 de 240 en 1 s
   DESCARTADOS 11 (siguen dando mas de un fragmento): CHEMBL4164295, CHEMBL4535036, CHEMBL4163381, CHEMBL4535469, CHEMBL4587627, CHEMBL4093916, CHEMBL4448793, CHEMBL4174088, CHEMBL4174548, CHEMBL4167276
   Un positivo que se cae sin avisar baja el AUC sin que se note, asi que se repone con el siguiente de la lista.
   preparados 11 de 11 en 0 s
   anadidos: 240 (objetivo total 300)
   pchembl de los anadidos: max 11.00 | min 8.66

SULFONAMIDA (esto es lo que se puede decir, y lo que no)
   sulfonamidas entre los nuevos: 223 de 240
   ficheros con algun hidrogeno polar (HD/HS): 223 de 223
   OJO CON ESE NUMERO, y por eso no se lee como un aviso:
     `re_hd` mira si el PDBQT trae ALGUN hidrogeno polar, no si el
     hidrogeno esta EN el nitrogeno del sulfonamida. Con cualquier
     otro donador en la molecula ya da que si, asi que no decide
     nada. El banco de 60 lo leyo como '60 de 60 mal' y ese aviso
     era un artefacto de la comprobacion, no una medida.
     Mapear el N del sulfonamida por indice tampoco vale: meeko
     reordena los atomos al escribir el PDBQT (fue el hallazgo del
     control de CDK2 el 29 de septiembre).
   Queda pendiente una comprobacion de verdad, con correspondencia de
   atomos deducida como en `diagnostico_orden_poses.py`. No se
   inventa aqui.

positivos escritos: activos_hca2_ampliado.csv (300)

AUDITORIA DEL BANCO AMPLIADO
   esqueletos de los 300 positivos: 297
   senuelos actuales: 2403
   senuelos que comparten esqueleto con ALGUN positivo: 10
   senuelos limpios (los que se usarian al medir): 2393 -> decoys_hca2_ampliado.csv
   excluidos: 10 -> decoys_hca2_excluidos.csv
   ATENCION: esos 10 inflarian el AUC por construccion (son casi
   indistinguibles de un positivo). Se dejan fuera al medir; el banco original NO se toca.
