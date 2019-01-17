#! /usr/bin/python3
# Authors: Pablo Baeyens y Antonio Checa
# iterativo.py
# Fichero para métodos de cálculo iterativos para resolución de sistemas de ecuaciones no lineales

import argparse
import math
import numpy as np

def e(N, i):
  """Devuelve el i-ésimo vector de la base usual de R^N.
  Argumentos posicionales:
  - N: Dimensión del vector
  - i: Posición en la base usual

  Devuelve:
  Vector de dimensión N con componentes todas nulas salvo en la posición i"""
  # Más eficiente que np.eye(1,N,i) de acuerdo con
  v = np.zeros(N)
  v[i] = 1
  return v


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
    ders.append((f(x + delta*e(N,i)) - fx)/delta[i])

  return np.vstack(ders)


def lm(f, inicial, objetivo, umbral = 1e-4, max_iter = 1000):
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


  I = np.eye(inicial.size, inicial.size) # identidad

  augment = 0.0001 # valor lambda
  x = inicial
  iters = 0
  epsilon = f(x)-objetivo
  norm = np.linalg.norm(epsilon)

  while norm > umbral and iters < max_iter:
    J = jacobiana(f, x)

    delta = np.linalg.solve(J.dot(J.T) + augment*I, -J.dot(epsilon))
    delta = [delta[0][i] for i in range(len(delta[0]))]
    candidate = x + delta
    cand_norm = np.linalg.norm(f(candidate) - objetivo)
    if cand_norm < norm:
      x = candidate
      epsilon = f(x)-objetivo
      norm = cand_norm
      augment /= 10
    else:
      augment *= 10
      if np.isinf(augment):
        augment = 0.001
    iters += 1

  return x, norm
