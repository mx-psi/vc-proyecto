#! /usr/bin/python3
# Authors: Pablo Baeyens y Antonio Checa
# iterativo.py
# Fichero para métodos de cálculo iterativos para resolución de sistemas de ecuaciones no lineales

import argparse
import math
import numpy as np

def jacobiana(f, x, delta = None):
  """Aproxima la traspuesta del jacobiano de f en x mediante el cálculo de
  diferencias finitas de primer orden en cada componente.

  Argumentos posicionales:
  - f: Función (callable) a la que calcular la derivada. Se asume que es C²
  - x: Punto en el que calcular el jacobiano

  Devuelve:
  - Aproximación de la traspuesta de la matriz jacobiana de f en x"""

  # Calcula tamaño de delta
  # en función de fórmula dada en Multiple View in Geometry (Apéndice A6.2)
  if delta is None:
    delta = np.maximum(1e-6, 1e-4*np.abs(x))
  N = x.size
  fx = f(x)

  ders = []
  for i in range(N): # Añade las derivadas parciales respecto de xi
    v = np.zeros(N)
    v[i] = delta[i] # Dirección de la derivada
    ders.append((f(x + v) - fx)/delta[i])
  return np.vstack(ders)


def lm(f, inicial, objetivo, umbral = 1e-4, max_iter = 100):
  """Implementa el algoritmo de Levenberg-Marquadt.
   Dada f, calcula x tal que |f(x)-objetivo| < umbral en un entorno de inicial o
   devuelve la mejor aproximación encontrada si no ha encontrado tal x en max_iter iteraciones.

  Argumentos posicionales:
  - f: Función
  - inicial: Estimación inicial de x
  - objetivo: valor en el codominio de f al que acercarse

  Argumentos opcionales:
  - umbral: Umbral de error aceptable
  - max_iter: Máximo número de iteraciones

  Devuelve:
  - x: Valor que minimiza |f(x)-objetivo|
  - err: |f(x)-objetivo|"""


  I = np.eye(inicial.size) # identidad

  augment = 0.0001 # valor lambda
  x = inicial # Mejor solución hasta el momento
  iters = 0 # Número de iteraciones
  epsilon = f(x)-objetivo
  norm = np.linalg.norm(epsilon) # Valor a minimizar

  while norm > umbral and iters < max_iter:
    if np.isinf(augment):
      # No se ha conseguido convergencia. Abortamos
      break

    JT = jacobiana(f, x) # Traspuesta de la jacobiana

    try:
      delta = np.linalg.solve(JT.dot(JT.T) + augment*I, -JT.dot(epsilon)).reshape((9,))
    except np.linalg.linalg.LinAlgError:
      # Matriz singular. Aumentamos el parámetro augment
      augment *= 10
      continue

    candidate = x + delta
    cand_norm = np.linalg.norm(f(candidate) - objetivo)
    if cand_norm < norm:
      x = candidate
      epsilon = f(x)-objetivo
      norm = cand_norm
      augment /= 10
    else:
      augment *= 10

    iters += 1

  return x, norm


if __name__ == "__main__":
  print("Este es un módulo auxiliar que no puede ejecutarse por sí solo.")
  print("Por favor, ejecuta 'python3 main.py'")
