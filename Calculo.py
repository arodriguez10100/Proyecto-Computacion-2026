#PARAMETROS PARA CALCULAR EL LAGRANGIANO

import sympy as smp
#Trabaja con vectores y arreglos de numeros
import numpy as np
#resuelve ecuaiones diferenciales con aproximaciones numericas
from scipy.integrate import solve_ivp

from Lagrangiano import lagrangiano

'''
Definimos una función que dados los parametros fisicos y condiciones iniciales, use la funcion Lagrangiano para resolverlo 
con aproximaciones numericas el sistema
'''


#Acá 'recuperamos' las expresiones simbólicas de la funcion Lagrangiano
    #Almacenamos en esta nueva funcion los simbolos de la funcion lagrangiano


#Lambdify transforma expresiones de SymPy en funciones lambda para calcular valores numéricos
#(simbolico a numerico)
#acc reciben números y calculan las aceleraciones angulares

#Tenemos la aceleracion angular de la primera y segunda masa

'''
solve_ivp no trabaja directamente con ecuaciones diferenciales de segundo orden.
Tenemos que transformar las dos ecuaciones de segundo orden en u nsistema de cuatro ecuaciones de primer orden.
'''




    

