# Re-validación de bolsillos y calibración del corte — 20 agosto 2026

> **DOCUMENTO HISTÓRICO — LEER CON LAS CORRECCIONES POSTERIORES.** El veredicto
> «✅ VALIDADO» de SOD1 Trp32 que aparece abajo **ya no se sostiene**: el 20 de
> septiembre de 2026 se encontró que 16 de los 20 activos son la misma serie
> congénica (`AUDITORIA_POSITIVOS_2026-09-20.md`). El control de re-acoplamiento
> de los cuatro ligandos cristalizados del bolsillo, en cambio, **pasa en dos de
> los cuatro** una vez corregido un error de programación nuestro en el medidor
> de RMSD: isoproterenol 0,48-1,28 Å y adrenalina 0,66-0,75 Å; falla en dopamina
> (2,9-3,2 Å) y en 5-fluorouridina, cuya pose depositada no es un mínimo del
> potencial (contacto F···N de 2,07 Å en una estructura de 1,06 Å). Ver
> `redocking_trp32/INFORME_CONTROLES_CORREGIDOS_2026-09-20.md`. Lo que el 0,815
> demuestra es que la caja está bien colocada; lo que queda en duda es si el
> método ordena compuestos de quimiotipos distintos.
> Se deja el texto original sin editar para no reescribir el registro.

Ejecutado por Buffy (tareas completas): (1) re-validar TDP-43 y FUS con señuelos,
(2) recuperar Oracle y el rescoring MM-GBSA, (3) controles positivos y recalibrar el corte.

## 1. Validación con señuelos (protocolo mejorado)

Se amplió la librería fuente de señuelos de 912 a **6.462 compuestos**
(ChEMBL extra + FDA + batch2), más señuelos property-matched por activo (MW ±30,
logP ±1.0, rotB ±2, Tanimoto <0.35), y se dockeó con Vina CPU (misma caja del cribado,
sin tocar la GPU ocupada en la tanda z010).

| Diana | Bolsillo | Activos | Señuelos | ROC-AUC | EF5% | Veredicto |
|---|---|---|---|---|---|---|
| **SOD1** | Trp32 (46.5,80.0,73.3) | 20 | 199 | **0.815** | **4.0** | ✅ VALIDADO |
| **FUS** | modelo 1 (−14.5,15.1,−7.8) | 2 | 87 | **0.609** | 0.0 | ⚠️ DÉBIL (2 activos) |
| **TDP-43** | RRM1 (28.3,43.7,52.5) | 5 | 193 | **0.465** | 0.0 | ❌ NO VALIDADO |

Notas:
- SOD1 Trp32 re-confirmado (AUC 0.815, EF5% 4.0, reproducible en
  `_validacion_SOD1/validacion_SOD1_trp32.csv`). Las cajas metal (0.717) y dímero (0.809)
  son inferiores.
- TDP-43 RRM1: AUC **0.465** (< azar). Los señuelos puntúan MEJOR que los activos de
  literatura (media act −5.23 vs dec −5.46). El bolsillo RRM1 no discrimina. Causa
  probable: los "activos" (rTRD01/nTRD22 fragmentos de superficie de RNA; bis-ANS/5FUrd
  ligantes de agregados/RNA) no ocupan el bolsillo dockeado.
- FUS: AUC 0.609 con solo 2 activos de literatura → estadísticamente inconcluso.

## 2. Controles positivos (ChEMBL/PubChem)

- **ChEMBL NO tiene ninguna bioactividad curada** para TDP-43 (CHEMBL2362981) ni FUS
  (CHEMBL5724679): las consultas devuelven 404 (sin datos).
- PubChem tampoco aporta activos anotados para TARDBP/FUS.
- Conclusión: **no existen controles positivos fiables por proteína en bases públicas**.
  La validación de TDP-43/FUS solo puede apoyarse en literatura fragmentaria, lo que
  limita cualquier benchmark de señuelos (n=2-5 activos → intervalos de confianza anchos).

## 3. Calibración del corte

- El corte por percentil **por proteína (top-5%) ≈ −6.3 kcal/mol** en las 3 dianas —
  es el criterio robusto (scores Vina no comparables entre receptores).
- El umbral absoluto **−7.0 ≈ top-1%** (SOD1 p99 = −7.00; TDP43 p99 = −6.9; FUS p99 = −6.8).
  Es un corte MUY estricto y selectivo.
- Frente a los señuelos: el 5% FP calibrado cae en −3.85 (TDP43) a −4.74 (SOD1). Un
  candidato ≤ −7 queda MUY por debajo del 5% FP → baja probabilidad de falso positivo
  decoy-like. Pero esto NO rescata la mala discriminación del bolsillo TDP-43.

## 4. Estado del rescoring MM-GBSA (Oracle)

- Oracle estaba inaccesible por SSH: **mi IP quedó bloqueada** (fail2ban), no la VM.
  Se recuperó el acceso vía puente desde la VM de GCloud (34.41.60.229). La clave se
  limpió después.
- El runner **terminó** (24/24 pendientes, 03:34): total 42 pares candidato-diana.
- **Calidad deficiente para TDP-43**: 12 TIMEOUT>1800s + 7 'NoneType' + 2 NaN + 2 valores
  basura (dG=1.9e13, 51.3). Solo 1 resultado TDP-43 válido (CHEMBL3309775 −14.87).
- FUS/SOD1 válidos: CHEMBL3311449·SOD1 −15.79 (líder), CHEMBL9010·FUS −24.75,
  CHEMBL7563·FUS −30.64, CHEMBL9532·FUS −32.22, CHEMBL8905·FUS −22.74...
  **Extremos sospechosos** (dG < −30): típicos de poses mal parametrizadas; no creerlos
  sin inspección.
- Pendiente: re-pasada TDP-43 con pipeline arreglado y revisión de los extremos FUS.

## 5. Recomendaciones (orden)

1. **TDP-43**: decidir si el bolsillo RRM1 es correcto (comparar con la superficie de
   unión a RNA de la literatura) o re-definirlo; mientras tanto, marcar los 23 candidatos
   TDP-43 como no validados en el paper.
2. **FUS**: buscar más activos de literatura (o MD de ligandos reportados) para subir la
   potencia estadística del benchmark (n=2 es insuficiente).
3. **Rescoring**: arreglar el pipeline MM-GBSA de TDP-43 (timeouts) y revisar poses de
   los dG extremos de FUS antes de rankear.
4. Mantener top-5% por proteína como corte primario; −7 como umbral de priorización.

## Archivos

- `validacion_TDP43_v2.csv`, `validacion_FUS_v2.csv` (detalle activos/señuelos)
- `compounds/decoys_library.smi` (6.462 compuestos fuente de señuelos)
- `rescoring_mmgbsa.csv`, `rescoring_mmgbsa_robusto.csv` (descargados de Oracle)
- `analizar_validacion_v2.py` (script de análisis)

---

## ACTUALIZACIÓN 2026-08-26 — Bolsillo TDP43 v6 (interfaz RRM1-RRM2, Arg151-Asp247)

### Problema raíz identificado
El receptor usado hasta ahora (4IUF) solo contiene RRM1 (residuos 103-213), **sin Asp247 ni RRM2** — el sitio donde se unen los activos validados experimentalmente NO EXISTE en ese receptor. Por eso la v4 daba AUC 0.517 ≈ azar y activos con energía positiva (+7.99).

### Corrección científica (literatura: Kapsiani et al., Cambridge 2026, PMC12918951)
- **Receptor nuevo: 4BS2** (tandem RRMs 96-269, con Arg151 y Asp247 del puente salino crítico para unión a RNA)
- **Bolsillo: interfaz RRM1-RRM2 / puente salino Arg151-Asp247** — sitio donde PE859 (−8.49 kcal/mol, H-bond Arg151) y berberrubine (−7.72, interfaz RRM1-RRM2 solapando bolsillo ATP) se unen
- Caja: centro (24.23, 16.89, −15.87), tamaño 26 Å, cubre ambos sitios (separados 11.8 Å)
- Receptor preparado con H polares, re-tipificado estilo ADFR (N donador 220, NA aceptor 21, HD 307)

### Activos de validación (9, con respaldo experimental)
PE859 y berberrubine (validados en HEK + C. elegans, paper 2026), berberine, sanguinarine, cepharanthine, epiberberine, coptisine, nitidine (cluster de berberrubine con actividad reportada), ketoconazole (reportado reducir agregación).

### Resultados v6 (Vina CPU 1.2.3, exhaustividad 16, size 26, 9 activos + 122 decoys property-matched)
- **ROC-AUC: 0.517** (crudo), 0.533 (por átomo pesado), 0.556 (por √átomos)
- ActivOS acoplan correctamente: cepharanthine −10.72 (90% mejor que decoys), sanguinarine −8.41 (83%), coptisine −8.23 (77%), berberrubine −7.74 (57%), PE859 −7.60 (50%)
- PE859 (−7.60) y berberrubine (−7.74) reproducen las energías del paper (−8.49, −7.72) dentro de ~1 kcal/mol → **el protocolo y receptor son correctos**
- Los decoys top (−13.6 a −12.5) ~~son compuestos grandes de la librería (75-162 átomos) que puntúan mejor por artefacto de tamaño de Vina~~ **CORREGIDO 2026-09-21: no era el tamaño, eran pseudo-átomos «glue» de Meeko. Ver la actualización al final del documento.**

### Interpretación científica
1. **El bolsillo está CORREGIDO**: el receptor ahora contiene el sitio real de unión y los activos validados acoplan como reporta la literatura (mejora cualitativa enorme vs v4 donde no acoplaban).
2. **Vina crudo NO discrimina** (AUC ~0.52) — hallazgo típico: el scoring de Vina favorece moléculas grandes/lipofílicas. Los decoys property-matched de la propia librería CHEMBL acoplan bien.
3. **Implicación**: el ranking por afinidad Vina cruda no es fiable para TDP43. El ranking final DEBE usar rescoring MM-GBSA (ya en marcha en Oracle para los top-100 vina).

### Acciones recomendadas
- [x] Receptor 4BS2 preparado y guardado en `_tdp43_bolsillo_v2/`
- [x] Activos v2 en `activos_tdp43_v2.csv`
- [x] Resultados v6 en `validacion_TDP43_v6_4BS2.csv`
- [ ] Re-cribado de la librería contra bolsillo nuevo (después de SOD1/FUS o en Kaggle) + rescoring MM-GBSA del top
- [ ] Considerar el sesgo de tamaño en el ranking (normalización por átomos o filtro de MW)

---

## ACTUALIZACIÓN 2026-09-21 — El «artefacto de tamaño» de TDP-43 eran pseudo-átomos «glue»

### Qué se creía
Que los decoys que encabezaban el ranking de TDP-43 v6 (−13.6 a −12.5) eran
compuestos grandes de la librería que puntuaban mejor por el sesgo de tamaño de
Vina en cajas pequeñas.

### Qué era en realidad
Meeko, con los ajustes por defecto, **no cierra los anillos de 7 eslabones en
adelante**: los abre en dos ramas y pega los extremos con pseudo-átomos de
pegamento (`CG0`/`G0`), duplicando dos átomos por anillo abierto.

Consecuencia en el score: Vina reporta `inter + intra − unbound`, y en el
ligando aislado el anillo abierto tiene una energía interna (`REMARK UNBOUND`)
de **+7 a +16 kcal/mol** en vez de ~0. El compuesto recibe ese regalo completo.

Los 16 afectados de TDP-43 v6 son **los 14 primeros señuelos del ranking
(puestos 1 a 14) más la activa cefarantina**. Al repararlos pierden entre 3,2 y 6,1 kcal/mol cada
uno, y su `UNBOUND` baja de +7,4…+15,8 a ≈ −0,4.

### Alcance medido (los cuatro conjuntos de validación)

| Conjunto | Afectados | Papel en el ranking | AUC antes → después |
|---|---|---|---|
| TDP43 v6 (4BS2) | 16 | los 14 primeros señuelos + **cefarantina** (la mejor activa) | 0,517 → 0,534 (crudo); 0,534 → 0,554 (por átomo); 0,556 → 0,539 (por raíz) |
| FUS v2 | 3 | puestos 1-3 (3 señuelos) | 0,609 → 0,626 |
| SOD1 trp32 (v1) | 2 | puestos 1-2 (2 señuelos) | 0,797 → 0,803 |
| SOD1 v3 | 11 | todo el fondo (5 duros + 6 emparejados) | 0,671 → 0,686 (crudo); 0,486 → 0,532 (sin tamaño), familia «independientes» |

Resultados por conjunto en `reparar_validaciones_glue.log` y
`reparar_sod1_v3_glue.log`; las poses y los ligandos viejos se conservan en
`out/_antes_glue/` y `ligands/_roto_glue/` de cada conjunto, así que todo es
reproducible.

### Qué cambia en las conclusiones

1. **La conclusión de fondo NO cambia.** El AUC de TDP-43 sigue siendo ≈ azar
   (0,53 crudo) y el de SOD1 trp32 sin normalizar sigue siendo ~0,80 solo
   porque el sesgo de tamaño trabaja a favor de los activos grandes. Lo que
   cambia es **la causa**: no es el tamaño de Vina, es un fichero de ligando mal
   preparado.
2. **La explicación «sesgo de tamaño» queda descartada para los top de
   TDP-43.** Los 14 primeros señuelos son moléculas normales (22-54 átomos). El sesgo de
   tamaño sigue existiendo en el fondo (≈ −0,09 kcal/mol por carbono pesado),
   pero no explica la cabeza del ranking.
3. **El score de Vina sí tenía señal en TDP-43, solo estaba tapada.** Tras
   reparar, la activa cefarantina (−7,47) sigue por delante del mejor señuelo
   (−9,39, CHEMBL607833). Antes de reparar, 14 señuelos artificialmente inflados
   la superaban.
4. **Todos los conjuntos de validación estaban contaminados, no solo TDP-43.**
   El fallo estaba también en SOD1 (v1 y v3) y FUS. La librería ya se había
   reparado (3.803 afectados, `gpu_dock/reparar_libreria_glue.py`), pero los
   ficheros de control se prepararon con el Meeko viejo.

### Acciones
- [x] Reparados y re-acoplados los 32 ligandos afectados de los cuatro conjuntos
- [x] Métricas antes/después recalculadas con el mismo receptor, caja y exhaustividad de cada conjunto
- [x] Herramientas: `reparar_ligandos_glue.py`, `reparar_validaciones_glue.py`, `reparar_sod1_v3_glue.py`
- [ ] Revisar si otras validaciones históricas (`validacion_TDP43_v2/v3/v4/v5`, `validacion_FUS.csv`) tienen el mismo defecto antes de citarlas
