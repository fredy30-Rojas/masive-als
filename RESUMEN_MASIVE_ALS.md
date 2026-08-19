# MASIVE-ALS — Documento completo del proyecto (18 agosto 2026)

## QUÉ ES
MASIVE-ALS = Massive Virtual Screening for Amyotrophic Lateral Sclerosis Drug Candidates.
Cribado virtual masivo de compuestos contra 3 proteínas clave de la ELA: TDP-43, SOD1, FUS.
Herramientas: AutoDock Vina (Vina-GPU 2.1), AlphaFold-Multimer, GROMACS.
Meta: 3-5 fármacos candidatos para validación experimental.
Acceso abierto (CC-BY 4.0, Zenodo).
Lema: "No busco publicar un artículo. Busco vivir."
Proyecto de Fredy Rojas (diagnóstico ELA). Investigación sin financiación, con infraestructura gratuita.

## CARPETA
C:\Users\Fredy\masive-als\

## ESTADO ACTUAL (18 agosto 2026)
1. Tanda z001 en GPU local (RTX 4080, Vina-GPU 2.1 compilado para Windows): COMPLETA — 4.984 compuestos x 3 proteínas (4.978 acoplados por proteína).
2. TDP-43: COMPLETA — 4.978/4.984 acoplados. Mejor afinidad -8.6 kcal/mol (CHEMBL6356). Resultados: gpu_dock/tanda_z001/resultados_z001.csv.
3. SOD1: COMPLETA — 4.978/4.984 acoplados. Mejor afinidad -9.5 kcal/mol (CHEMBL4559945).
4. FUS: COMPLETA — 4.978/4.984 acoplados. Mejor afinidad -8.2 kcal/mol (CHEMBL6207).
5. Embudo de candidatos RECALCULADO (revisión científica 18/08): POR PROTEÍNA (top-5% propio de cada diana, no global — corrige el sesgo 138/17/16) + filtro CNS (TPSA<=90, MW<=450). z001 completa (14.952): 133 tras PAINS+Lipinski/Veber (TDP43 46, SOD1 45, FUS 42) → 42 con CNS en analysis/candidatos_filtrados.csv. Total histórico (32.715): 244 → 76 CNS, con 6 FÁRMACOS APROBADOS FDA: clozapina (SOD1 -7.5, CNS-MPO 4.8), fluoresceína (SOD1 -7.1), rucaparib (SOD1 -7.1/FUS -6.4), ketazolam y perampanel (FUS -6.4) — en analysis/candidatos_total_cns.csv. Scripts nuevos: analysis/cns_filtro.py. SMILES de 27 fármacos corregidos (el primer hit de ChEMBL era sal/N-óxido, e.g. clozapina; validación de 98 nombres: 69 OK vs librería FDA, 27 corregidos). Poses PDBQT mapeadas a cada candidato (columna pose_pdbqt). Enriquecimiento de SMILES vía API ChEMBL en background (566 IDs del top-5% total). RECONCILIACIÓN 18/08: el CSV de producción es analysis/candidatos_filtrados.csv (42, con poses PDBQT), metodología de la revisión (top-5% por proteína + CNS). Una re-implementación independiente con umbrales más estrictos y sin filtro CNS (38 candidatos, _candidatos_v4.csv) es una variante exploratoria — documentada y superada en el paper. Nota CNS-MPO: el score es aproximación de 4 componentes (no los 6 de Wager con pKa/logD); solo 10/42 llegan a >=4 (mejor CHEMBL9347 4.42); el filtro CNS principal es TPSA<=90/MW<=450.
6. Total regenerado 18/08 01:20: 32.715 pares únicos (resultados_vinagpu_total.csv) = histórico (27.721) + z001 CORRECTA INTEGRADA (4.984×3 = 14.952 pares, 100%; SOD1 -9.5 ahora presente). Nota: es HISTÓRICO, mezcla bolsillos previos; la tanda vigente de referencia es z001. Publicado en GitLab (masive-als-data) 18/08 01:20. Push anteriores 01:07 con 32.699.
6b. Los 6 compuestos faltantes de z001 (B boro / Si silicio: tipos no soportados por Vina) YA RECUPERADOS: parche de tipos en _patch_tipos_pdbqt2.py (añadido caso 'Si' minúscula) + re-dock 18/08 01:14 → 18/18 pares (CHEMBL4553125: TDP-43 -7.4, SOD1 -7.2). Informe: tandas_z001/faltantes_z001.txt (marcado RESUELTO).
6c. Validación de señuelos EXTENDIDA (18/08 ~03:25, revisión científica punto 1): además de SOD1 (Trp32 AUC 0.815 ✅), se validaron TDP-43 y FUS con activos de la literatura (TDP-43: rTRD01, nTRD22 [NMR], bis-ANS, Congo Red, 5-FUrd, isoproterenol → AUC 0.515 ⚠️; FUS: dehidroximetilflazina, cleroindicina C [ML+MD] → AUC 0.545 ⚠️). Ambos cerca del azar: las cajas de TDP-43 (RRM1) y FUS (modelo 1) NO discriminan binders conocidos de señuelos — marcadas como provisionales hasta re-benchmark del bolsillo (igual que se corrigió SOD1). Detalle: analysis/validacion_TDP43.csv y validacion_FUS.csv. Nota: set de activos pequeño (especialmente FUS, poca literatura de ligandos directos).
7. Kaggle v15: cuota semanal GPU (30h) agotada; auto-lanzador en Oracle reintenta el push hasta el sábado 22/08 00:00 UTC.
8. Oracle ARM (79.72.57.253): pipeline cribando ChEMBL 24/7, watchdog de tandas encadena y fusiona resultados.
9. Conversión librería ZINC 2M: COMPLETA 400/400 tandas (z001-z400) = 1.934.817 PDBQT 3D en C:\Users\Fredy\zinc2m\tandas_tgz (18/08 05:03). Una tanda (z070) quedó truncada por interrupción del empaquetado y se re-empaquetó desde sus 4.834 PDBQT (validado gzip, sin huecos).

## RECEPTORES
- TDP-43: dominio RRM1 (TDP43_RRM1_clean.pdb → gpu_dock/TDP43.pdbqt). Caja de la literatura: centro (28.3, 43.7, 52.5), tamaño 25. NOTA paper: el receptor real NO es el PDB completo (6b1n/4IUF), es el dominio RRM1 limpio.
- SOD1: PDB 1HL5 (tetrámero → gpu_dock/SOD1.pdbqt). Caja: centro (46.5, 80.0, 73.3), tamaño 22 — cubre región de agregación C-terminal (residuos 1-33, 96-112, 131-153), decisión deliberada anti-agregación, NO el sitio catalítico Cu/Zn.
- FUS: PDB 6G99 — ATENCIÓN: estructura NMR con 20 modelos, Vina la rechaza; usar modelo 1 (reparado con Open Babel). Caja: centro (-14.5, 15.1, -7.8), tamaño 25.
- Sitios de docking corregidos a los de la literatura (agosto 2026): energías corregidas -8.2 a -9.5 kcal/mol

## COMPUESTOS / LIBRERÍAS
- 87 fármacos FDA aprobados + 50 bioactivos CHEMBL de ELA (primeras tandas)
- 3.311 fármacos aprobados ChEMBL (2.978 PDBQT)
- 827 bioactivos ChEMBL batch2 (813 PDBQT)
- 5 lotes ChEMBL adicionales: fase III (1.073), fase II (6.631), fase I (923), naturales (7.629), batch3 pChEMBL>=6
- 107 tandas (~47.000 ligandos) en GitHub masive-als-data
- Librería ZINC 2M en conversión (400 tandas)
- Total de pares cribados a la fecha: 32.715 pares únicos (histórico 27.721 + z001 completa 14.952, deduplicados e integrados 18/08)

## INFRAESTRUCTURA
- Oracle Cloud ARM (79.72.57.253): pipeline 24/7, watchdog tandas_als_watchdog.py, memoria de resultados /mnt/extra/masive-als/resultados
- Kaggle GPU (P100): 30h/semana gratis, cuota se renueva sábado 00:00 UTC
- Modal: 30 usd/mes gratis (token yograbotodo), crédito agotado este mes
- Google Colab T4: bloqueado temporalmente (renueva ~24h)
- PC local RTX 4080 Laptop (12GB): Vina-GPU 2.1 Windows, tanda z001 ~10-12h resumible
- GitHub repo: fredy30-Rojas/masive-als-data (ligandos y resultados)

## SUPERCOMPUTADORAS SOLICITADAS (seguimiento)
- MareNostrum 5 (BSC-CNS, Barcelona): 200.000h GPU + 50.000h CPU, programa RES (deadline 15 sept, requiere correo institucional)
- Frontier (OLCF, Oak Ridge): Director's Discretionary, 80.000 GPU-node-hours AMD MI250X
- Leonardo (EuroHPC, Cineca, Italia)
- LUMI (EuroHPC, Finlandia)
- Colaboración en contacto: Ana Martínez y Carmen Gil (CIB-CSIC) — "lo revisamos todo en septiembre"; Enric (BSC)
- Dra. Povedano (Hospital Bellvitge): informes de progreso enviados

## CRONOGRAMA
- Sep 2026: fase 0, preparación entorno
- Oct-Nov 2026: cribado fase 1 (5M compuestos)
- Dic 2026: cribado fase 2 (5M compuestos)
- Ene 2027: validación dinámica molecular (top 1.000)
- Feb 2027: análisis y publicación

## ESTRUCTURA DE CARPETAS
- paper/ — paper_masive_als.md (artículo completo)
- solicitudes/ — Frontier OLCF, Leonardo EuroHPC, LUMI + correos
- cartas/ — Bellvitge, IRB Barcelona, VHIR Vall Hebron
- doctora/ — documentos, correo seguimiento, generador de audio
- patente/ — borrador_patente.md
- presentacion/ — presentacion_masive_als.html
- web/ — portal.html (tema oscuro), index.html, presentacion.html
- gpu_dock/ — tandas de docking GPU (tanda_z001 con resultados, configs, logs)
- proteins/ — receptores
- compounds/ — ligandos
- resultados/, results/, md_results/ — resultados de cribados
- analysis/ — merge_results.py, prepare_md.py, analyze_md.py
- scripts/, src/, tools/ — scripts de conversión, fábrica de tandas (fabrica_100_tandas.py), descargadores ChEMBL
- deploy/ — deploy.sh, nginx-sites.conf
- dominio/ — masive-als.json
- colaboradores/ — contactos
- prep/ — preparación de datos

## DOCUMENTOS CLAVE
- README.md — descripción orientada a MareNostrum 5 (BSC-CNS)
- paper/paper_masive_als.md — artículo completo
- 6_FASES_DEL_PROYECTO.txt — fases del proyecto
- guia_registro_res.md — guía registro RES
- informe_progreso_707.md — informe de progreso (707 docks)
- reflexion_buffy_9agosto2026.md — reflexión sobre error del firewall

## NOTAS
- Fredy tiene ELA, escribe con control ocular Tobii 4C + OptiKey Pro.
- Regla: no enviar correos sin aprobación previa de Fredy.
- Pipeline corregido en agosto 2026: receptor FUS 6G99 (NMR) rechazado por Vina; sitios de docking corregidos a la literatura; energías reales.
- Tratamiento: datos abiertos, reproducibilidad total.
