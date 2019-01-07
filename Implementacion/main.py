#! /usr/bin/python3
# Authors: Pablo Baeyens y Antonio Checa
# main.py
# Fichero principal para estimación de homografías

import argparse
import math
import numpy as np

def normaliza(puntos, inv = False):
  """Normaliza un conjunto de puntos.
  Argumentos posicionales:
  - puntos: Conjunto de vectores de dimensión n.
  Argumentos opcionales:
  - inv: Flag que indica si se devuelve la transformación directa o inversa
  Devuelve:
  - Conjunto de puntos con centroide cero y distancia media al centroide $\sqrt{2}$
  - Matriz que lleva el conjunto original al conjunto devuelto (o inversa)"""

  centroid = np.sum(puntos)/len(puntos)
  dist = np.linalg.norm(puntos - centroid) # TODO: Arreglar

  pass

def inicialHom(corr):
  """Obtiene una estimación inicial de la homografía
  Argumentos posicionales:
  - corr: Lista de pares de correspondencias
  Devuelve:
  - Estimación inicial de homografía que ajusta estas correspondencias"""

  orig, T_orig = normaliza(np.fromiter(x for x,y in corr,float))
  dest, T_dest = normaliza(np.fromiter(y for x,y in corr, float), inv = True)

  # TODO: Coger de las prácticas

  return T_dest*H*T_orig

def iteracion(H, err, corr):
  """Realiza un paso del algoritmo iterativo"""
  pass

def getHom(corr, umbral = 1e-4, N = 1000):
  """Obtiene una homografía a partir de una lista de correspondencias.
  Argumentos posicionales:
  - corr: Lista de pares de puntos en coordenadas inhomogéneas
  Argumentos opcionales:
  - umbral: Umbral por debajo del cuál el error es aceptable
  - N: Número máximo de iteraciones
  Devuelve:
  - Matriz 3x3 que representa la homografía
  - Error residual de la homografía
  """

  H = inicialHom(corr)
  err = math.inf
  n = 0

  while err > umbral or n < N:
    H, err = iteracion(H,err,corr)
    n += 1

  return H, err


def showHom(im1, im2):
  """Muestra un mosaico de dos imágenes."""
  corr   = getCorrespondences(im1, im2)
  H, err = getHom(corr)
  # TODO: Mostrar usando funciones de OpenCV


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  subparsers = parser.add_subparsers(help="Acceso al modo manual", dest="modo")
  manual = subparsers.add_parser('manual')
  manual.add_argument('archivo1', type=str, help="Archivo de imagen 1")
  manual.add_argument('archivo2', type=str, help="Archivo de imagen 2")
  args = parser.parse_args()

  if args.modo == "manual":
    print("Modo manual")
    # TODO: Cargar imágenes
  else:
    print("Modo de ejemplo")
    # TODO: Elegir imágenes por defecto
