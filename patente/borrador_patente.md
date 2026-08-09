# BORRADOR DE SOLICITUD DE PATENTE
# Título: Método de cribado molecular masivo para la identificación de compuestos terapéuticos contra la Esclerosis Lateral Amiotrófica
# Inventor: Fredy Rojas Gutiérrez
# Fecha: 9 de agosto de 2026

---

## 1. TÍTULO DE LA INVENCIÓN

Método computacional de cribado molecular masivo en infraestructura de supercomputación pública para la identificación de compuestos terapéuticos contra la agregación de TDP-43, la desestabilización de SOD1 y la transición de fase aberrante de FUS en la Esclerosis Lateral Amiotrófica.

---

## 2. CAMPO TÉCNICO

La presente invención se refiere al campo del descubrimiento computacional de fármacos, específicamente a un método de cribado virtual masivo que combina predicción estructural mediante inteligencia artificial (AlphaFold-Multimer), docking molecular acelerado por GPU (AutoDock-GPU) y validación por dinámica molecular (GROMACS) ejecutado en infraestructura pública de supercomputación.

---

## 3. ANTECEDENTES

Actualmente no existe un tratamiento curativo para la ELA. Los métodos tradicionales de descubrimiento de fármacos requieren entre 10 y 15 años y una inversión superior a 1,000 millones de euros por compuesto. La presente invención reduce el tiempo de la fase de descubrimiento a 6 meses mediante el uso de supercomputación pública.

---

## 4. DESCRIPCIÓN DETALLADA

El método comprende las siguientes etapas:

a) Generación de 100,000 conformaciones proteicas de TDP-43, SOD1 y FUS mediante AlphaFold-Multimer con 12 ciclos de reciclaje.

b) Descarga y conversión de 10 millones de compuestos de las librerías ZINC20, DrugBank y Enamine REAL al formato PDBQT.

c) Cribado virtual masivo mediante AutoDock-GPU en una configuración de 200 GPUs simultáneas (NVIDIA H100), con un rendimiento de 400,000 docks por segundo.

d) Selección de los 1,000 compuestos con mejor energía de unión para cada diana proteica.

e) Validación mediante dinámica molecular de 1 microsegundo con GROMACS, incluyendo cálculo MM-GBSA de energía libre de unión, RMSD y radio de giro.

f) Selección final de 3-5 compuestos cabeza de serie con energía de unión < -10 kcal/mol y RMSD < 2 Å.

---

## 5. REIVINDICACIONES

1. Método de cribado molecular masivo caracterizado por la combinación secuencial de AlphaFold-Multimer, AutoDock-GPU y GROMACS sobre infraestructura de supercomputación pública con al menos 200 GPUs simultáneas.

2. Método según reivindicación 1, donde las proteínas diana son TDP-43, SOD1 y FUS.

3. Método según reivindicación 1, donde el número de compuestos cribados es de al menos 10 millones.

4. Compuestos identificados mediante el método de las reivindicaciones anteriores para su uso en el tratamiento de la ELA.

---

## 6. RESUMEN

Método computacional que reduce el tiempo de descubrimiento de fármacos contra la ELA de 14 años a 6 meses mediante cribado masivo de 10 millones de compuestos contra 3 proteínas diana utilizando 200 GPUs en un superordenador público.

---

*Documento preparatorio. Se recomienda consultar con un agente de patentes antes de la presentación formal ante la OEPM (Oficina Española de Patentes y Marcas).*
