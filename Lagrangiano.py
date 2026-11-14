#FORMULA DEL LAGRANGIANO

import sympy as smp #Calculo Simbolico

'''
Definimos las variables en juego del sistema de manrera simbolica en sympy,
esto para obtener las ecuaciones de Lagrange.
Las longitudes de las cuerdas (L1, L2; que no tienen peso), las masas (m1, m2),  la gravedad, 
el tiempo t y dejamos por fuera el parametro de rozamiento porque analizamos un sistema ideal
'''


#Definimos una funcion para calcular el lagrangiano

def lagrangiano():

    #DEFINIR VARIABLES
    #con el fin de obtener las ecuaciones en términos del nombre de las variables (i.e t**2)

    #tiempo masas largos y gravedad
    t , g = smp.symbols('t, g') #Tiempo y gravedad respectivamente
    m1, m2 = smp.symbols('m1, m2') #Masas 1 y 2
    L1, L2 = smp.symbols('L1, L2') #Longitudes de las cuerdas 1 y 2

    # smp.pprint(t**2) pretty printer

    #Falta definir los angulos que son funciones de tiempo 
    theta1, theta2 = smp.symbols('theta1, theta2', cls=smp.Function)

    #Las escribo como funciones de tiempo
    theta1= theta1(t)
    theta2 = theta2(t)

    #smp.pprint(theta2) pretty printer

    #Hallo las primeras y segundas derivadas de theta
    theta1_d = smp.diff(theta1, t)
    theta2_d = smp.diff(theta2, t)
    theta1_dd = smp.diff(theta1_d, t)
    theta2_dd = smp.diff(theta2_d, t)

    #Definimos las coordenadas en el plano cartesiano usando newton
    x1 = L1*smp.sin(theta1)
    y1 = -L1*smp.cos(theta1)
    x2 = L1*smp.sin(theta1)+L2*smp.sin(theta2)
    y2 = -L1*smp.cos(theta1)-L2*smp.cos(theta2)


    #Energia cinetica
    E1 = 1/2 * m1 * (smp.diff(x1, t)**2 + smp.diff(y1, t)**2)
    E2 = 1/2 * m2 * (smp.diff(x2, t)**2 + smp.diff(y2, t)**2)
    E = E1+ E2

    #Energia potencial
    U1 = m1*g*y1
    U2 = m2*g*y2
    U = U1+U2

    #Lagrangiano
    L = E - U

    #ESCRIBIMOS LAS ECUACIONES DE EULER-LAGRANGE

    #smp.diff(L, thetan_d) calcula la derivada parcial del lagrangiano respecto a la vel angular y la otra respecto a t
    # Para theta 1
    EL1 = smp.diff(smp.diff(L, theta1_d), t) - smp.diff(L, theta1)
    # Para theta 2
    EL2 = smp.diff(smp.diff(L, theta2_d), t) - smp.diff(L, theta2)

    # Sympy automaticamente simplifica las ecuaciones
    EL1 = smp.simplify(EL1)
    EL2 = smp.simplify(EL2)

    # RESOLVEMOS LAS ECUACIONES DE MOVIMIENTO
    soluciones = smp.solve(
        [EL1, EL2],
        [theta1_dd, theta2_dd],
        dict=True) #dict=True sympy retorna las soluciones diccionarios la primera la resuelve respecto a theta 1 y así

    #Almacenamos las dos soluciones como variables separadas

    theta1_dd_sol = soluciones[0][theta1_dd]
    theta2_dd_sol = soluciones[0][theta2_dd]

    #accedo al primer dicceionario de la lista [0] y da el valor asociado a la variable [theta]
    
    return (
        theta1_dd,
        theta2_dd,
        theta1,
        theta2,
        theta1_d,
        theta2_d,
        m1, m2, L1, L2, g,
        soluciones
    )

