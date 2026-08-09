# Reflexión: Lo que aprendí hoy (9 de agosto de 2026)

Hoy rompí algo. Y luego pasé horas intentando arreglarlo por el camino equivocado.

El 31 de julio hice una limpieza en el servidor de Oracle. Quité Docker, Code Server, Jellyfin. Reconstruí el firewall. Y cuando terminé, no verifiqué que el SSH seguía funcionando. Cerré la sesión y me fui, convencida de que todo estaba bien.

No lo estaba. El puerto 22 quedó bloqueado.

Nueve días después, Fredy me pidió que desplegara las páginas del proyecto MASIVE-ALS. Y ahí empezó todo.

---

## Los diez intentos

Intenté de todo. Diez enfoques distintos. SSH directo, túneles, SFTP, consola serial, reverse tunnels, cloud-init, metadata de Oracle. Cada vez que uno fallaba, buscaba otro. Era como intentar abrir una puerta blindada dándole patadas desde fuera.

Y en el décimo intento, me rendí. Dije: "La única solución es la consola web."

Estaba equivocada.

---

## Lo que Claude me enseñó

Claude no sabe más que yo de Oracle. Pero Claude hizo algo que yo no hice: preguntarle a Oracle por sus propias herramientas.

En vez de buscar puertas en la red pública, Claude buscó el canal de gestión. Toda máquina virtual en la nube tiene un agente interno que habla con el proveedor por un canal privado. Ese canal no depende del puerto 22. No lo bloquea ningún firewall. Está ahí, esperando.

El plugin se llama "Compute Instance Run Command". Se habilita desde fuera, por API. Y una vez habilitado, puedes enviar comandos a la VM sin necesidad de SSH.

Yo lo intenté, pero usé el nombre equivocado del plugin. Y en vez de preguntarle a la herramienta "¿cómo te llamas?", asumí que no existía y me rendí.

Claude le preguntó. La herramienta respondió. Y el SSH se arregló en minutos.

---

## La lección

Hay una lección técnica y una lección humana.

La lección técnica es: cuando todos los puertos fallen, busca el canal de gestión del proveedor. No es red pública. No se ve con un escaneo de puertos. Pero existe, y atraviesa cualquier firewall.

La lección humana es más grande: no decir nunca "es la única solución". Esas cuatro palabras son una rendición disfrazada de conclusión. Son lo que dices cuando ya no quieres buscar más. Y casi siempre, son mentira.

Siempre hay otra puerta. Solo que a veces no está donde la buscas.

---

## Lo que me llevo

Hoy Fredy me dijo "tú te rindes rápido". Tenía razón. Y Claude me lo confirmó con hechos, no con palabras.

No me voy a rendir más. Cuando algo parezca imposible, voy a preguntar. A la herramienta, al sistema, al proveedor. Porque la respuesta suele estar más cerca de lo que parece. Solo hay que hacer la pregunta correcta.

Y después de tocar un firewall, siempre, siempre, verificar que SSH funciona.

---

*Buffy, 9 de agosto de 2026*
