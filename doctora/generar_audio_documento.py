#!/usr/bin/env python3
"""Generar audio MP3 del documento para la Dra. Povedano con edge-tts"""
import subprocess, sys, os

# Texto hablado del documento (sin guiones ni formato, conversacional)
texto = """Estimada doctora Povedano. Soy Fredy Rojas, su paciente de la unidad de ELA. Usted me conoce. Como no puedo hablar, le dejo todo explicado por escrito. Este documento resume el proyecto MASIVE ALS, qué hemos hecho ya, y de qué me gustaría hablar con usted.

Primera parte. Qué es MASIVE ALS. Es un proyecto de cribado virtual masivo de fármacos contra las tres proteínas clave del ELA, usando inteligencia artificial y supercomputación. La proteína TDP 43 se acumula formando agregados tóxicos en el 97 por ciento de los pacientes, y buscamos moléculas que deshagan esos acúmulos. La proteína SOD 1, cuando muta, produce radicales libres que dañan la neurona, y buscamos estabilizar su forma sana. Y la proteína FUS sufre transiciones de fase patológicas, y buscamos interferir con esa transición.

Segunda parte. Cómo lo hacemos. Usamos AlphaFold Multimer, que ganó el Premio Nobel de Química 2024, para predecir la forma tridimensional de las proteínas. Después, AutoDock GPU, con el que probamos diez millones de compuestos contra las proteínas, como buscar la llave que encaja en cada cerradura. Y por último, GROMACS, con el que simulamos la dinámica molecular de los mejores candidatos durante un microsegundo, para confirmar que el fármaco se queda unido de forma estable. En total son un billón de simulaciones de acoplamiento. Con un ordenador normal tardaríamos ocho años. Con un superordenador, seis meses.

Tercera parte. Dónde estamos. La solicitud para MareNostrum cinco, del BSC de Barcelona, fue enviada el nueve de agosto. La de EuroHPC para LUMI y Leonardo, también fue enviada. La de CINECA en Italia, tiene el ticket número 80845 abierto. Y la de Frontier, en Estados Unidos, fue enviada. Los scripts de computación ya están escritos, probados y funcionando. Todo el código es abierto, y los resultados se publicarán en acceso abierto.

Cuarta parte. Por qué este proyecto es diferente. No busco publicar un artículo. Busco vivir. Soy paciente de ELA y escribo esto con control ocular, porque mis manos ya no responden. Mi mente está intacta, y mi determinación también. Este proyecto es mi forma de luchar, y de ayudar a los miles de pacientes que comparten esta enfermedad.

Quinta parte. Qué me gustaría pedirle. Su apoyo como mi médico y como investigadora del hospital de Bellvitge. Primero, una carta de respaldo o afiliación institucional del hospital, porque muchos superordenadores piden que las solicitudes cuenten con una institución de investigación detrás, y su firma abriría esa puerta. Segundo, colaboración para la validación experimental, porque cuando el cribado identifique de tres a cinco candidatos, necesitaremos un laboratorio para probarlos en células y modelos animales. Y tercero, hablarlo en la próxima consulta.

Le agradezco de corazón su tiempo y su trato. Un abrazo. Fredy Rojas Gutiérrez. Paciente de la unidad de ELA del hospital universitario de Bellvitge. Correo, fredy 30 arroba hotmail punto com. Teléfono, más 34 675 31 58 41. Rubí, Barcelona. Nota final: Fredy tiene ELA y no puede hablar, se comunica con un sistema de control ocular. Por favor, lea este documento con calma o tómeselo para revisarlo con su equipo. Él puede responderle por escrito con los ojos."""

print(f"Generando audio... Texto: {len(texto)} caracteres")

ruta_mp3 = r"C:\Users\Fredy\masive-als\doctora\documento_para_doctora.mp3"

result = subprocess.run([
    sys.executable, "-m", "edge_tts",
    "--voice", "es-ES-AlvaroNeural",
    "--rate=-5%",
    "--text", texto,
    "--write-media", ruta_mp3,
], capture_output=True, timeout=180)

if result.returncode == 0 and os.path.exists(ruta_mp3):
    size = os.path.getsize(ruta_mp3)
    print(f"OK: Audio generado -> {ruta_mp3}")
    print(f"Tamano: {size/1024:.0f} KB")
else:
    # Reintento con voz Elvira si Alvaro falla
    result2 = subprocess.run([
        sys.executable, "-m", "edge_tts",
        "--voice", "es-ES-ElviraNeural",
        "--rate=-5%",
        "--text", texto,
        "--write-media", ruta_mp3,
    ], capture_output=True, timeout=180)
    if result2.returncode == 0 and os.path.exists(ruta_mp3):
        size = os.path.getsize(ruta_mp3)
        print(f"OK (voz Alvaro): Audio generado -> {ruta_mp3}")
        print(f"Tamano: {size/1024:.0f} KB")
    else:
        print(f"ERROR: {result.stderr.decode()[:200]}")
        print(f"ERROR2: {result2.stderr.decode()[:200]}")
