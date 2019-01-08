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
  - puntos: Array np de vectores de dimensión n.
  Argumentos opcionales:
  - inv: Flag que indica si se devuelve la transformación directa o inversa
  Devuelve:
  - Conjunto de puntos con centroide cero y distancia media al centroide $\sqrt{2}$
  - Matriz que lleva el conjunto original al conjunto devuelto (o inversa)"""

  # Centroide
  centroide = np.mean(puntos, axis = 0)

  # Distancia media al centroide
  distancia_media = np.mean(np.linalg.norm(puntos - centroide, axis = 1))

  # Factor de escalado para normalizar distacia media
  escalado = math.sqrt(2)/distancia_media

  if inv:
    S = np.diag([1/escalado,1/escalado,1])
    tx,ty = centroide
    T = np.array([[1,0,tx], [0,1,ty], [0,0,1]])
    M = T@S
  else:
    tx,ty = - centroide
    T = np.array([[1,0,tx], [0,1,ty], [0,0,1]])
    S = np.diag([escalado,escalado,1])
    M = S@T

  # Calcula valores a devolver
  normalizados = escalado*(puntos - centroide)
  return normalizados, M


def inicialHom(origs, dests):
  """Obtiene una estimación inicial de la homografía
  Argumentos posicionales:
  - origs: Lista de orígenes
  - dests: Lista de coordenadas de destinos
  Devuelve:
  - Estimación inicial de homografía que ajusta estas correspondencias"""

  orig_n, T_orig = normaliza(origs)
  dest_n, T_dest = normaliza(dests, inv = True)

  v = np.zeros(3)
  A = None

  # Para cada correspondencia concatena la matriz Ai
  for src, dst in zip(orig_n, dest_n):
    # TODO: Asumo que las correspondencias vienen dadas en coordenadas inhomogéneas
    src_h = np.append(src, 1)
    f1 = np.concatenate((v, - src_h, dst[1]*src_h))
    f2 = np.concatenate((src_h, v, -dst[0]*src_h))
    if A is None:
      A = np.vstack((f1, f2))
    else:
      A = np.vstack((A, f1, f2))

  # Halla SVD
  *_, V = np.linalg.svd(A)
  H = V[-1,:].reshape((3,3))

  return T_dest@H@T_orig

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
    puntos = np.array([[0,3],[4,3],[-2,3],[-2,-2]])

    norm, T = normaliza(puntos, inv = True)
    centroide = np.mean(norm, axis = 0)
    distancia_media = np.mean(np.linalg.norm(norm, axis = 1))
    print(centroide)
    print(distancia_media)

    corr = np.array([[p,p] for p in puntos], float)
    print(inicialHom(puntos, puntos))
    # TODO: Elegir imágenes por defecto
