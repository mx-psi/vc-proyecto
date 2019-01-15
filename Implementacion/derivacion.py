#! /usr/bin/python3
# Authors: Pablo Baeyens y Antonio Checa
# derivacion.py
# Fichero para cálculo aproximado de derivadas mediante diferencias finitas

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
  v = np.zeros(N)
  v[i] = 1
  return v

def derivada(f, x):
  """Aproxima el jacobiano de f en x mediante el cálculo de
  diferencias finitas de primer orden en cada componente.
  Argumentos posicionales:
  - f: Función (callable) a la que calcular la derivada. Se asume que es C²
  - x: Punto en el que calcular el jacobiano
  Devuelve:
  - Aproximación de la matriz jacobiana de f en x"""

  # Calcula tamaño de delta
  # en función de fórmula dada en Multiple View in Geometry (Apéndice A6.2)
  delta = np.maximum(1e-6, 1e-4*x)
  N = x.size
  fx = f(x)

  ders = []
  for i in range(N): # Añade las derivadas parciales respecto de xi
    ders.append((f(x + delta*e(N,i)) - fx)/delta[i])

  return np.vstack(ders).T
