import numpy as np
import sympy as smp  # Cálculo Simbolico
from scipy.integrate import solve_ivp

'''
Definimos las variables en juego del sistema de manera simbólica en SymPy,
esto para obtener las ecuaciones de Lagrange.

Las longitudes de las cuerdas (L1, L2; que no tienen peso),
las masas (m1, m2), la gravedad, el tiempo t.
Dejamos por fuera el parámetro de rozamiento porque analizamos
un sistema ideal.
'''

#Definimos las variables

#Tiempo, masas, longitudes y gravedad
t, g = smp.symbols('t, g')
m1, m2 = smp.symbols('m1, m2')
L1, L2 = smp.symbols('L1, L2')

#Definimos los ángulos como funciones del tiempo
theta1, theta2 = smp.symbols('theta1, theta2', cls=smp.Function)

theta1 = theta1(t)
theta2 = theta2(t)

#Hallamos las primeras y segundas derivadas de theta

theta1_d = smp.diff(theta1, t)
theta2_d = smp.diff(theta2, t)

theta1_dd = smp.diff(theta1_d, t)
theta2_dd = smp.diff(theta2_d, t)

#Definimos las coordenadas en el plano cartesiano

x1 = L1*smp.sin(theta1)
y1 = -L1*smp.cos(theta1)

x2 = L1*smp.sin(theta1) + L2*smp.sin(theta2)
y2 = -L1*smp.cos(theta1) - L2*smp.cos(theta2)

#ENERGIA CINÉTICA

E1 = 1/2 * m1 * (smp.diff(x1, t)**2 + smp.diff(y1, t)**2)
E2 = 1/2 * m2 * (smp.diff(x2, t)**2 + smp.diff(y2, t)**2)

E = E1 + E2

#ENERGIA POTENCIAL

U1 = m1*g*y1
U2 = m2*g*y2

U = U1 + U2

#LAGRANGIANO

L = E - U

#Ecuaciones de Euler-Lagrange

#Para theta1
EL1 = smp.diff(smp.diff(L, theta1_d), t) - smp.diff(L, theta1)

#Para theta2
EL2 = smp.diff(smp.diff(L, theta2_d), t) - smp.diff(L, theta2)

#Simplificamos las ecuaciones
EL1 = smp.simplify(EL1)
EL2 = smp.simplify(EL2)

#Resolvemos las ecuaciones de movimiento (acá usamos un diccionario)
soluciones = smp.solve(
    [EL1, EL2],
    [theta1_dd, theta2_dd],
    dict=True
)

#Almacenamos las dos soluciones
theta1_dd_sol = smp.simplify(soluciones[0][theta1_dd])
theta2_dd_sol = smp.simplify(soluciones[0][theta2_dd])


#Lambdify transforma expresiones de SymPy en funciones lambda para calcular valores numéricos
#(simbolico a numerico) funciones ecuaciones de lagrange
#Volvemos a resolver las ecuaciones de E-L pero para los parámetros
soluciones = smp.solve([EL1, EL2], theta1_dd, theta2_dd)
LEF1 = smp.lambdify((theta1, theta2, theta1_d, theta2_d, t, m1, m2, L1, L2, g), soluciones[theta1_dd])
LEF2 = smp.lambdify((theta1, theta2, theta1_d, theta2_d, t, m1, m2, L1, L2, g), soluciones[theta2_dd])

#Parámetros y condiciones iniciales para los dos péndulos
cond_inicial_1 = [2.0+10e-6, 2.0, -1.0, 3.0]  #theta1, theta1_d, theta2, theta2_d, 
cond_inicial_2 = [2.0, 2.0, -1.0, 3.0]
m1_val = 2
m2_val = 2
L1_val = 1
L2_val = 1
g_val = 9.8


'''
solve_ivp de scipy no trabaja directamente con ecuaciones diferenciales de segundo orden.
Tenemos que transformar las dos ecuaciones de segundo orden en un sistema de cuatro ecuaciones de primer orden.
'''

#Parámetros de tiempo y un vector z que contendrá las 4 ecuaciones
def sistema_ec(t, z):

    theta1 = z[0]
    theta1_d = z[1]

    theta2 = z[2]
    theta2_d = z[3]
#Calculamos las aceleraciones angulares con los parámetros ingresados
    theta1_dd = LEF1(
        theta1, theta2,
        theta1_d, theta2_d,
        t, m1_val, m2_val, L1_val, L2_val, g_val
    )

    theta2_dd = LEF2(
        theta1, theta2,
        theta1_d, theta2_d,
        t, m1_val, m2_val, L1_val, L2_val, g_val
    )
#Devuelve una lista con cuatro elementos
    return [
        theta1_d,
        theta1_dd,
        theta2_d,
        theta2_dd
    ]

#Definimos durante cuánto tiempo queremos simular el movimiento del péndulo doble y en cuántos puntos queremos observarlo.
t_inicial = 0
t_final = 20


#Acá usando la estructura de linspace le decimos desde que instante comenzamos hasta qué segundos 
#y cuantos puntos va a tomar en ese tiempo

tiempo = np.linspace(t_inicial, t_final, 1000)

#Resolución numérica del primer péndulo usando solve.ivp() con el método Runge-Kutta
sol1 = solve_ivp(
    sistema_ec,
    (t_inicial, t_final),
    cond_inicial_1,
    t_eval=tiempo,
    method='RK45',
    rtol=1e-9,
    atol=1e-9
)
sol2 = solve_ivp(
    sistema_ec,
    (t_inicial, t_final),
    cond_inicial_2,
    t_eval=tiempo,
    method='RK45',
    rtol=1e-9,
    atol=1e-9
)

#Queremos extraer y separar los resultados de las dos simulaciones numéricas, 
#sol1.y y sol2.y: contiene los valores de las variables del sistema

#Péndulo 1

theta1_sol_1 = sol1.y[0]
theta2_sol_1 = sol1.y[2]

theta1_d_sol_1 = sol1.y[1]
theta2_d_sol_1 = sol1.y[3]


#Péndulo 2

theta1_sol_2 = sol2.y[0]
theta2_sol_2 = sol2.y[2]

theta1_d_sol_2 = sol2.y[1]
theta2_d_sol_2 = sol2.y[3]

#Importamos las librerías necesarias para graficar y realizar la animación

import matplotlib.pyplot as plt
#Viene de Matplotlib, para permite crear animaciones
from matplotlib.animation import FuncAnimation

#Creamos una figura que contenga dos graficos dentro de una: 1 fila y 2 columnas y el ancho 15 y el alto 6
fig, axs = plt.subplots(1, 2, figsize=(15, 6))

#Gráfica del péndulo 1
axs[0].plot(theta1_sol_1, theta1_d_sol_1,
            color='blue', label='Masa 1')

axs[0].plot(theta2_sol_1, theta2_d_sol_1,
            color='red', label='Masa 2')

axs[0].set_title('Diagrama de fase: Péndulo 1')
axs[0].set_xlabel('Posición angular (rad)')
axs[0].set_ylabel('Velocidad angular (rad/s)')
axs[0].legend()
axs[0].grid()


#Gráfica del péndulo 2
axs[1].plot(theta1_sol_2, theta1_d_sol_2,
            color='blue', label='Masa 1')

axs[1].plot(theta2_sol_2, theta2_d_sol_2,
            color='red', label='Masa 2')

axs[1].set_title('Diagrama de fase: Péndulo 2')
axs[1].set_xlabel('Posición angular (rad)')
axs[1].set_ylabel('Velocidad angular (rad/s)')
axs[1].legend()
axs[1].grid()

#Mostramos las gráficas
plt.tight_layout()
plt.show()

'''
Para esta sección realizamos la animación:
Acá en principio, pasamos de los ángulos a coordenadas cartesianas para animar. Cabe destacar que 
theta1_sol_1 (y 2) contiene todos los ángulos que obtuvo la primera simulación (y la segunda)
'''

#Péndulo 1

x1_sol_1 = L1_val * np.sin(theta1_sol_1)
y1_sol_1 = -L1_val * np.cos(theta1_sol_1)

x2_sol_1 = x1_sol_1 + L2_val * np.sin(theta2_sol_1)
y2_sol_1 = y1_sol_1 - L2_val * np.cos(theta2_sol_1)

#Péndulo 2

x1_sol_2 = L1_val * np.sin(theta1_sol_2)
y1_sol_2 = -L1_val * np.cos(theta1_sol_2)

x2_sol_2 = x1_sol_2 + L2_val * np.sin(theta2_sol_2)
y2_sol_2 = y1_sol_2 - L2_val * np.cos(theta2_sol_2)


#Acá creamos una figura y un sistema de ejes donde posteriormente vamos a hacer la animación
#Creamos una figura cuadrada de 8x8
fig, ax = plt.subplots(figsize=(8, 8))

#Configuramos la apariencia y que queden a la misma escala
ax.set_xlim(-2.5, 2.5)
ax.set_ylim(-2.5, 2.5)
ax.set_aspect('equal')
ax.grid()

ax.set_title("Comparación de dos pendulos dobles con condiciones iniciales 'cercanas'")

#Primer pendulo: azul
linea1, = ax.plot([], [], 'o-', color='blue', lw=3, markersize=8)

#Segundo pendulo: rojo
linea2, = ax.plot([], [], 'o-', color='red', lw=3, markersize=8)

#Trayectorias
trayectoria1, = ax.plot([], [], color='blue', alpha=0.3)
trayectoria2, = ax.plot([], [], color='red', alpha=0.3)

#Acá definimos para comenzar la animación, comienza vacío
def iniciar():

    linea1.set_data([], [])
    linea2.set_data([], [])

    trayectoria1.set_data([], [])
    trayectoria2.set_data([], [])

    return linea1, linea2, trayectoria1, trayectoria2

#Acá va a ir actualizando las posiciones
def actualizar(i):

    #Primer pendulo
    linea1.set_data(
        [0, x1_sol_1[i], x2_sol_1[i]],
        [0, y1_sol_1[i], y2_sol_1[i]]
    )

    #Segundo pendulo
    linea2.set_data(
        [0, x1_sol_2[i], x2_sol_2[i]],
        [0, y1_sol_2[i], y2_sol_2[i]]
    )

    #Trayectorias de las masas finales
    trayectoria1.set_data(
        x2_sol_1[:i+1], y2_sol_1[:i+1]
    )

    trayectoria2.set_data(
        x2_sol_2[:i+1], y2_sol_2[:i+1]
    )

    return linea1, linea2, trayectoria1, trayectoria2


animacion = FuncAnimation(
    fig,
    actualizar,
    frames=len(sol1.t),
    init_func=iniciar,
    interval=20,
    blit=True
)

plt.show()

'''
Acá verificamos si el codigo está bien, calculamos la energía mecánica inicial y final de cada péndulo y calculamos
el error porcentual entre ambas, como el sistema es conservativo, la energía debería conservarse así que debería ser razonablemente chico
'''
#Definimos una función que reciba ciertos parámetros
def energia_mecanica(th1, th2, th1_d, th2_d):

    Ec = (
        0.5 * m1_val * (L1_val * th1_d)**2
        + 0.5 * m2_val * (
            (L1_val * th1_d)**2
            + (L2_val * th2_d)**2
            + 2 * L1_val * L2_val * th1_d * th2_d
            * np.cos(th1 - th2)
        )
    )

    Ep = -g_val * (
        (m1_val + m2_val) * L1_val * np.cos(th1)
        + m2_val * L2_val * np.cos(th2)
    )

    return Ec + Ep

#Creamos variables que almacene la energía mecánica del primer péndulo
Em1 = energia_mecanica(
    theta1_sol_1, theta2_sol_1,
    theta1_d_sol_1, theta2_d_sol_1
)
#Creamos variables que almacene la energía mecánica del segundo péndulo
Em2 = energia_mecanica(
    theta1_sol_2, theta2_sol_2,
    theta1_d_sol_2, theta2_d_sol_2
)

#Imprimimos los valores
print("Pendulo 1")
print("Energia inicial:", Em1[0])
print("Energia final:", Em1[-1])
#Imprimimos el error porcentual entre ambas energías
print("Error (%):", abs((Em1[-1]-Em1[0])/Em1[0])*100)

print("Péndulo 2")
print("Energia inicial:", Em2[0])
print("Energia final:", Em2[-1])
#Imprimimos el error porcentual entre ambas energías
print("Error (%):", abs((Em2[-1]-Em2[0])/Em2[0])*100)
