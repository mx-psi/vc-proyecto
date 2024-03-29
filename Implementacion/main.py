#! /usr/bin/python3
# Authors: Pablo Baeyens y Antonio Checa
# main.py
# Fichero principal para estimación de homografías

import argparse
import math
import numpy as np
import cv2
import iterativo
import auxiliar
from scipy import optimize

def Ai(orig, dest):
  """Calcula ecuaciones que impone una correspondencia sobre las incógnitas de una homografía
  Argumentos posicionales:
  - orig: Punto de origen
  - dest: Punto de destino
  Devuelve:
  - Matriz 2x9 que representa el sistema de ecuaciones
  """

  x, y = dest
  orig = np.append(orig,1)
  zeros = np.zeros(3)
  r1 = np.concatenate((zeros, -orig, y*orig))
  r2 = np.concatenate((orig, zeros, -x*orig))
  return np.vstack((r1,r2))


def C_Hx(orig, dest, h):
  """Calcula C_H(X) donde X = (orig, dest) y H es la homografía asociada al vector h.
  Argumentos posicionales:
  - orig: Punto de origen
  - dest: Punto de destino
  - h: Matriz de homografía redimensionada como vector de R⁹
  Devuelve:
  - Evaluación de C_H(X) con X = (orig, dest)"""
  return Ai(orig, dest).dot(h)


def JJT(orig, dest, h):
  """Función auxiliar para el cálculo del error de Sampson.
  Argumentos posicionales:
  - orig: Punto de origen en coordenadas inhomogéneas
  - dest: Punto de destino en coordenadas inhomogéneas
  - h: Matriz de homografía redimensionada como vector de R⁹

  Devuelve:
  - JJ.T, donde J es la matriz jacobiana de C_H(X)
  """

  JT = iterativo.jacobiana( # Cálculo exacto por ser multilineal
    lambda X: C_Hx(X[:2], X[2:], h), # Función C_H
    np.hstack([orig, dest]), # Punto X en la variedad
    delta = np.ones(4)) # Delta = 1

  return JT.T.dot(JT)

def error_sampson_corr(orig, dest, h):
  """Calcula el error de Sampson para una correspondencia.
  Argumentos posicionales:
  - orig: Punto de origen en coordenadas inhomogéneas
  - dest: Punto de destino en coordenadas inhomogéneas
  - h: Matriz de homografía redimensionada como vector de R⁹
  Devuelve:
  - El error de Sampson para la correspondencia"""

  epsilon = C_Hx(orig, dest, h) # Definimos epsilon como la evaluación de C_H(X)
  lamb = np.linalg.solve(JJT(orig, dest, h), -epsilon)  # Resolvemos la ecuación de lambda
  error_samps = np.transpose(epsilon).dot(-lamb)    # Calculamos el error de Sampson como epsilon*(-lambda)
  return error_samps


def error_sampson(origs, dests, h):
  """Calcula el error de Sampson para un conjunto de correspondencias.
  Argumentos posicionales:
  - corr: Iterable con pares de puntos origen, destino en coordenadas inhomogéneas
  - h: Matriz de homografía redimensionada como vector de R⁹
  Devuelve:
  - El error de Sampson para el conjunto de correspondencias"""

  # Para cada punto de correspondencia sumamos su error de Sampson
  err = 0
  for i in range(len(origs)):
    err += error_sampson_corr(origs[i], dests[i], h)
  return err


def normaliza(puntos, inv = False):
  """Normaliza un conjunto de puntos.
  Argumentos posicionales:
  - puntos: Array Numpy de forma (N,2)
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

  if inv: # Calcula transformación inversa
    S = np.diag([1/escalado,1/escalado,1])
    tx,ty = centroide
    T = np.array([[1,0,tx], [0,1,ty], [0,0,1]])
    M = T.dot(S)
  else: # Calcula transformación directa
    tx,ty = - centroide
    T = np.array([[1,0,tx], [0,1,ty], [0,0,1]])
    S = np.diag([escalado,escalado,1])
    M = S.dot(T)

  # Calcula valores a devolver
  normalizados = escalado*(puntos - centroide)
  return normalizados, M


def inicialHom(origs, dests):
  """Obtiene una estimación inicial de la homografía mediante el algoritmo DLT
  Argumentos posicionales:
  - origs: Lista de orígenes
  - dests: Lista de coordenadas de destinos
  Devuelve:
  - Estimación inicial de homografía que ajusta estas correspondencias mediante algoritmo DLT"""

  orig_n, T_orig = normaliza(origs)
  dest_n, T_dest = normaliza(dests, inv = True)

  A = None

  # Para cada correspondencia concatena la matriz Ai
  for orig, dst in zip(orig_n, dest_n):
    if A is None:
      A = Ai(orig, dst)
    else:
      A = np.vstack((A, Ai(orig,dst)))

  # Halla SVD
  *_, V = np.linalg.svd(A)
  H = V[-1,:].reshape((3,3))

  # Devuelve composición con transformaciones de normalización
  return T_dest.dot(H).dot(T_orig)


def getHom(origs, dests):
  """Obtiene una homografía a partir de una lista de correspondencias.
  Argumentos posicionales:
  - origs, dests: Lista de pares de puntos en coordenadas inhomogéneas
  Devuelve:
  - Matriz 3x3 que representa la homografía
  """

  f = lambda h: error_sampson(origs, dests, h) # Minimiza error de Sampson
  inicial = inicialHom(origs, dests).reshape((9,)) # Valor inicial dado por DLT
  h, err = iterativo.lm(f, inicial, 0) # Aplica Levenberg-Marquadt
  return h.reshape((3,3))


def creaMosaico(archivo1, archivo2, s_x, s_y):
  """Muestra un mosaico a partir de dos imágenes en un canvas de dimensiones dadas.
  Argumentos posicionales:
  - archivo1, archivo2: Archivos de imágenes
  - s_x, s_y: Dimensiones en píxeles del canvas"""

  im1 = auxiliar.lee_imagen(archivo1,1)
  im2 = auxiliar.lee_imagen(archivo2,1)

  # Cogemos los descriptores de las imágenes y sus keypoints
  kp1, des1 = auxiliar.getKpAndDescriptors(im1)
  kp2, des2 = auxiliar.getKpAndDescriptors(im2)

  # Hallamos matches
  matcher12 = auxiliar.getMatchesLowe2NN(des1, des2)

  # Recogemos las listas de keypoints de cada matcher ordenadas
  orderSrcKp1, orderDstKp12 = auxiliar.getOrderedKeypoints(kp1, kp2, matcher12)

  # Se crea el canvas final de tamaño s_x x s_y con tres bandas y uint8 por elemento
  canvas_final = np.zeros((s_x,s_y,3), np.uint8)

  # hi, wi es la altura y anchura de la imagen central
  hi = im2.shape[0]
  wi = im2.shape[1]
  # Punto en el que anclamos la imagen central dentro del canvas
  p_x = 350
  p_y = 50

  # Homografía de la imagen central al canvas. Es solo una traslación de las esquinas
  # a un rectángulo del mismo tamaño centrado dentro del canvas
  # Esta homografía sabemos cuál es, ya que es solo una traslación. La escribimos:
  h_canvas = np.array([[1., 0., p_x*1.0], [0., 1., p_y*1.0], [0., 0., 1.]])


  # Sustituimos encontrar la homografía con el algoritmo de error de Sampson
  ordSrcMod = np.array([orderSrcKp1[i][0] for i in range(len(orderSrcKp1))])
  ordDstMod = np.array([orderDstKp12[i][0] for i in range(len(orderDstKp12))])
  h1 = getHom(ordSrcMod, ordDstMod)

  # Se crea la imagen final haciendo llamadas a warpPerspective de cada imagen
  # con sus transformaciones correspondientes

  canvas_final2 = canvas_final
  # El borderMode del primer warpPerspective se pone a Constant para el fondo
  # El resto de borderMode se ponen a TRANSPARENT para no pisar el resto de imágenes
  canvas_final2 = cv2.warpPerspective(im1,  h_canvas.dot(h1), (s_x, s_y), canvas_final, borderMode = cv2.BORDER_CONSTANT)
  canvas_final2 = cv2.warpPerspective(im2, h_canvas, (s_x, s_y), canvas_final2, borderMode = cv2.BORDER_TRANSPARENT)

  imgOrdSrc = ordSrcMod.copy()

  # Para la mejor visualización de las distancias geométricas entre los puntos y sus imágenes,
  # aquí se colorea con una línea roja la distancia entre la imagen de un punto y su correspondencia.
  for i in range(len(ordSrcMod)):
    imgSrc3 = h_canvas.dot(h1).dot([ordSrcMod[i][0], ordSrcMod[i][1],1])
    imgOrdSrc[i] = [imgSrc3[0]/imgSrc3[2], imgSrc3[1]/imgSrc3[2]]
    imgOrdDst3 = h_canvas.dot([ordDstMod[i][0], ordDstMod[i][1], 1])
    #Se filtran solo los puntos que aparecen cerca del canvas final para crear líneas
    if(imgOrdSrc[i][0] <= 2000 and imgOrdSrc[i][1] <= 2000 and imgOrdDst3[0]/imgOrdDst3[2] <= 2000 and imgOrdDst3[1]/imgOrdDst3[2] <= 2000):
      if(imgOrdSrc[i][0] >= -2000 and imgOrdSrc[i][1] >= -2000 and imgOrdDst3[0]/imgOrdDst3[2] >= -2000 and imgOrdDst3[1]/imgOrdDst3[2] >= -2000):
        #Se crean las líneas rojas de distancias
        cv2.line(canvas_final2, (int(imgOrdSrc[i][0]), int(imgOrdSrc[i][1])), (int(imgOrdDst3[0]/imgOrdDst3[2]),int(imgOrdDst3[1]/imgOrdDst3[2])), (0,0,255), 2)
  # Se pinta el canvas final
  auxiliar.pintaI(canvas_final2, "Mosaico final")


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description =
      "Ejecuta sin argumentos para ver un ejemplo o ejecuta con 'manual' para introducir tu propio ejemplo")
  subparsers = parser.add_subparsers(help="Modo manual: Ejemplifica obtención de homografía \ncon construcción de mosaico entre dos imágenes hallando sus puntos SIFT", dest="modo")
  manual = subparsers.add_parser('manual', description="Ejemplifica obtención de homografía \ncon construcción de mosaico entre dos imágenes hallando sus puntos SIFT")
  manual.add_argument('archivo1', type=str, help="Archivo de imagen 1")
  manual.add_argument('archivo2', type=str, help="Archivo de imagen 2")
  args = parser.parse_args()

  if args.modo == "manual":
    creaMosaico(args.archivo1, args.archivo2, 1300, 700)
  else:
    print("Modo de ejemplo con imágenes Yosemite")
    print("Obtiene puntos SIFT mediante OpenCV y calcula una homografía para los mismos con el error de Sampson")
    print("Ejecuta \n   \"python3 main.py manual -h\" \npara ver instrucciones del modo manual")
    creaMosaico("./imagenes/yosemite1.jpg", "./imagenes/yosemite2.jpg", 1200, 600)
