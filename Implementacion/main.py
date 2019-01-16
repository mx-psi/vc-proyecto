#! /usr/bin/python3
# Authors: Pablo Baeyens y Antonio Checa
# main.py
# Fichero principal para estimación de homografías

import argparse
import math
import numpy as np
import iterativo


def C_Hx(orig, dest, h):
  """Calcula C_H(X) donde X = (orig, dest) y H es la homografía asociada al vector h.
  Argumentos posicionales:
  - orig: Punto de origen en coordenadas TODO
  - dest: Punto de destino en coordenadas TODO
  - h: Matriz de homografía redimensionada como vector de R⁹
  Devuelve:
  - Evaluación de C_H(X) con X = (orig, dest)"""

  x_i, y_i, w_i = dest
  first_row = np.concatenate(([0,0,0], -w_i*orig, y_i*orig), axis=None)
  second_row = np.concatenate((w_i*orig, [0,0,0], -x_i*orig), axis=None)
  m = np.vstack((first_row, second_row))

  h_t = np.vstack(h)

  return m.dot(h_t)

def JJT(orig, dest, h):
  """Función auxiliar para el cálculo del error de Sampson.
  Argumentos posicionales:
  - orig: Punto de origen en coordenadas TODO
  - dest: Punto de destino en coordenadas TODO
  - h: Matriz de homografía redimensionada como vector de R⁹

  Devuelve:
  - JJ.T, donde J es la matriz jacobiana de C_H(X)
  """

  orig_e1 = orig.copy()
  orig_e1[0] = orig_e1[0] + 1

  orig_e2 = orig.copy()
  orig_e2[1] += 1

  dest_e3 = dest.copy()
  dest_e3[0] += 1

  dest_e4 = dest.copy()
  dest_e4[1] += 1

  original = C_Hx(orig, dest, h)
  parcial_1 = C_Hx(orig_e1, dest, h) - original
  parcial_2 = C_Hx(orig_e2, dest, h) - original
  parcial_3 = C_Hx(orig, dest_e3, h) - original
  parcial_4 = C_Hx(orig, dest_e4, h) - original

  parcial_1 = parcial_1.dot(np.transpose(parcial_1))
  parcial_2 = parcial_2.dot(np.transpose(parcial_2))
  parcial_3 = parcial_3.dot(np.transpose(parcial_3))
  parcial_4 = parcial_4.dot(np.transpose(parcial_4))

  return (parcial_1 + parcial_2 + parcial_3 + parcial_4)


def error_sampson_corr(orig, dest, h):
  """Calcula el error de Sampson para una correspondencia.
  Argumentos posicionales:
  - orig: Punto de origen en coordenadas TODO
  - dest: Punto de destino en coordenadas TODO
  - h: Matriz de homografía redimensionada como vector de R⁹
  Devuelve:
  - El error de Sampson para la correspondencia"""

  epsilon = C_Hx(orig, dest, h)
  lamb = np.linalg.solve(JJT(orig, dest, h), -epsilon)

  return np.transpose(epsilon).dot(-lamb)


def error_sampson(corr,h):
  """Calcula el error de Sampson para un conjunto de correspondencias.
  Argumentos posicionales:
  - corr: Iterable con pares de puntos origen, destino en coordenadas TODO
  - h: Matriz de homografía redimensionada como vector de R⁹
  Devuelve:
  - El error de Sampson para el conjunto de correspondencias"""

  err = 0
  for orig, dest in corr:
    err += error_sampson_corr(orig, dest, h)
  return err


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


def getHom(corr):
  """Obtiene una homografía a partir de una lista de correspondencias.
  Argumentos posicionales:
  - corr: Lista de pares de puntos en coordenadas inhomogéneas
  Devuelve:
  - Matriz 3x3 que representa la homografía
  - Error residual de la homografía
  """

  f = lambda h: error_sampson(corr, h) # Minimiza error de Sampson
  inicial = inicialHom(corr).reshape((9,)) # Valor inicial dado por DLT
  h, err = iterativo.lm(f, inicial, 0)

  return h.reshape((3,3)), err


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
    puntos = np.random.rand(4,2)

    norm, T = normaliza(puntos, inv = True)
    centroide = np.mean(norm, axis = 0)
    distancia_media = np.mean(np.linalg.norm(norm, axis = 1))
    print(centroide)
    print(distancia_media)

    corr = np.array([[p,p] for p in puntos], float)
    H = inicialHom(puntos, puntos)
    print(H/H[2,2])
    # TODO: Elegir imágenes por defecto
