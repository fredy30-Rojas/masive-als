# Propuesta de colaboracion
## MASIVE-ALS + Grupo de Ana Martinez y Carmen Gil (CIB-CSIC)

---

## Resumen

Ustedes son el grupo mas avanzado en TDP-43 del mundo. Descubrieron AP-2, el primer farmaco contra TDP-43 que llega a ensayo clinico. Tienen Molefy Pharma, designacion de medicamento huerfano, y la Fase 1 en marcha.

Nosotros hemos solicitado acceso a 4 supercomputadoras, tenemos un pipeline con AlphaFold-Multimer (Nobel 2024) listo para ejecutar, y la capacidad de cribar 10 millones de compuestos. Pero necesitamos su respaldo cientifico para que nos concedan el acceso.

Proponemos unir fuerzas para expandir su trabajo a escala masiva y a nuevas proteinas.

---

## Lo que ofrecemos

### 1. Expandir su busqueda a SOD1 y FUS

Ustedes atacaron TDP-43 con exito. Pero la ELA tiene otras dos proteinas clave que nadie esta cribando a esta escala:

- SOD1: causa el 20% de la ELA familiar. Estabilizar el dimero nativo.
- FUS: transicion liquido-solido patologica. Bloquear la fase aberrante.

Podemos buscar farmacos para las tres proteinas en paralelo.

### 2. Buscar la segunda generacion de AP-2

AP-2 es un gran farmaco. Pero puede tener analogos mejores.

Podemos cribar 10 millones de compuestos contra TDP-43 buscando:
- Moleculas con mayor afinidad que AP-2
- Mejor selectividad (menos efectos secundarios)
- Mejor perfil ADME (absorcion, distribucion, metabolismo, excrecion)
- Compuestos que superen la barrera hematoencefalica

### 3. AlphaFold-Multimer (Premio Nobel de Quimica 2024)

Podemos generar 100,000 conformaciones por proteina, explorando todo el paisaje conformacional. Esto permite encontrar sitios de union que con estructuras estaticas serian invisibles. Nuestro pipeline usa la version mas reciente de AlphaFold-Multimer, galardonada con el Nobel en 2024.

### 4. Escala masiva con supercomputacion

| Recurso | Metodos tradicionales | MASIVE-ALS |
|---|---|---|
| Compuestos cribados | Miles | 10 millones |
| Docks por segundo | Cientos | 500,000 |
| Proteinas | 1 (TDP-43) | 3 (TDP-43, SOD1, FUS) |
| Conformaciones por proteina | 1-10 estaticas | 100,000 (AlphaFold) |
| MD por candidato | 10-50 ns | 1 microsegundo |
| Tiempo total | Meses/anos en cluster | 6 meses en supercomputadora |

### 5. Supercomputadoras solicitadas (pendientes de aprobacion)

- MareNostrum 5 (BSC, Barcelona): 200,000 h GPU
- Frontier (OLCF, EE.UU.): 80,000 GPU-node-hours AMD MI250X
- Leonardo (CINECA, Italia): 80,000 h GPU NVIDIA A100
- LUMI (CSC, Finlandia): EuroHPC Development Access

### 6. Pipeline completo ya construido

- AlphaFold-Multimer: prediccion de estructuras 3D
- AutoDock-GPU: docking masivo optimizado para multi-GPU
- GROMACS + CHARMM36: dinamica molecular de alta precision
- MM-GBSA: calculo de energia libre de union

### 7. Acceso abierto y sin coste

Todo el codigo en GitHub (CC-BY 4.0). Todos los resultados en Zenodo.

No pedimos financiacion. No pedimos equipos. Solo respaldo cientifico.

---

## Lo que pedimos

1. Carta de respaldo del CIB-CSIC para las solicitudes de supercomputacion. Su aval cientifico puede ser decisive para que nos concedan el acceso.
2. Orientacion cientifica en la validacion de candidatos
3. Coautoria en publicaciones
4. Explorar colaboracion con Molefy Pharma si surgen candidatos prometedores

---

## Por que tiene sentido esta colaboracion

- Ustedes ya demostraron que el cribado virtual funciona (AP-2). Nosotros lo llevamos a escala masiva.
- Ustedes abrieron el camino en TDP-43. Nosotros lo extendemos a SOD1 y FUS.
- Su infraestructura (Molefy Pharma) esta lista para llevar candidatos a ensayo clinico.
- No competimos: ustedes aportan la validacion experimental, nosotros la computacion.
- AlphaFold-Multimer (Nobel 2024) no existia cuando empezaron. Ahora esta disponible.
- Su respaldo puede ser la diferencia entre que nos aprueben o nos rechacen las supercomputadoras.

---

## Cronograma propuesto

- Septiembre 2026: confirmacion de acceso a supercomputadoras
- Octubre-Noviembre 2026: cribado Fase 1 (5M compuestos)
- Diciembre 2026: cribado Fase 2 (5M compuestos)
- Enero 2027: validacion MD de los 1,000 mejores
- Febrero 2027: publicacion de resultados

---

## Contacto

Fredy Rojas Gutierrez
fredy_30@hotmail.com | +34 675 31 58 41
Rubi, Barcelona

Paciente de la Unidad de ELA del Hospital de Bellvitge
Dra. Monica Povedano (mpovedano@bellvitgehospital.cat)

---

"No busco publicar un articulo. Busco vivir."
