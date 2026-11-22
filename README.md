Simulación de dos péndulos dobles, donde se grafican los diagramas de fase de cada uno y posteriormente se realiza una animación con los dos péndulos durante 20 segundos.
Ambos péndulos tienen condiciones iniciales casi idénticas, excepto que el ángulo 1 inicial del primer péndulo se le suman 10e^-6 radianes, mostrando la sensibilidad
a las condiciones iniciales de los sistemas caóticos.
El código se puede dividir en cinco bloques:
1. Computación simbólica: Resolución de las ecuaciones del Lagrangiano de manera simbólica (usando la biblioteca SymPy).
2. Computación numérica: Ingreso de datos numéricos y condiciones iniciales a la expresión simbólica y resolver para las aceleraciones angulares
(usando la biblioteca de SciPy y NumPy).
4. Graficación de los diagramas de espacio de fase: Graficación de la velocidad angular en función de la posición angular (usando la biblioteca de Matplotlib).
5. Animación del movimiento de los péndulos: Se realizó en función de su posición y su altura durante 20 segundos (usando Funcanimation de la biblioteca de Matplotlib).
6. Corroboración del código: Calculando las energías mecánicas inicial y final de cada péndulo y calculando el error porcentual respectivo.
