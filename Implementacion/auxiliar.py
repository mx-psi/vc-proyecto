#! /usr/bin/python3
# Authors: Pablo Baeyens y Antonio Checa
# auxiliar.py
# Fichero auxiliar para mostrar un ejemplo con correspondencias entre imágenes
# Adaptado de la práctica 2 de Antonio

import numpy as np
import cv2
import math
import os.path

def lee_imagen(fileName, flagColor):
  """Lee una imagen.
  Argumentos posicionales:
  - fileName: Nombre del fichero del que leer la imagen
  - flagColor: Indica si se lee a color (1) o en escala de grises (0)
  Devuelve:
  - La imagen leída
  """

  if not os.path.exists(fileName):
    print("[Error irrecuperable] La imagen '{}' no existe.".format(fileName))
    exit(-1)

  return cv2.imread(fileName,flagColor)


def pintaI(im, nombre):
  """Pinta una imagen.
  Argumentos posicionales:
  - im: Imagen
  - nombre: Nombre de la ventana"""

  cv2.namedWindow(nombre, cv2.WINDOW_NORMAL)
  cv2.resizeWindow(nombre, im.shape[1], im.shape[0]) # Hacemos resize de la
  #ventana para que la imagen se vea en su tamaño
  cv2.imshow(nombre, im)
  cv2.waitKey(0)
  cv2.destroyAllWindows()


def getKpAndDescriptors(im):
  """Calcula puntos y descriptores.
  Argumentos posicionales:
  - im: Imagen de la que provienen los puntos
  Devuelve:
  - Puntos clave de la imagen y descriptores SIFT de los mismos"""
  sift = cv2.xfeatures2d.SIFT_create(1000, 3, 0.08, 8, 1.6)
  kp = sift.detect(im)
  return sift.compute(im, kp)


def getMatchesLowe2NN(des1, des2):
  """Realiza un match con el ratio test de Lowe y un knnMatch con k=2
  Argumentos posicionales:
  - des1, des2: Descriptores de origen y destino
  Devuelve:
  - Matches que pasan este filtrado"""

  # Hacemos un matcher normal y realizamos un knnMatch con k=2 de los descriptores
  matcher_2nn = cv2.BFMatcher()
  matches_2nn = matcher_2nn.knnMatch(des1, des2, k=2)

  # Hacemos el ratio test de Lowe
  final_matches = []
  for match1, match2 in matches_2nn:
    # Solo si la divisón de las distancias es menor a 0.8, aceptamos el match
    # De esta forma descartamos aquellos que tienen más de un keypoint muy parecido
    # en los descriptores, y en los que podríamos tener dudas de cuál coger
    if match1.distance/match2.distance <= 0.8:
      final_matches.append(match1)
  #Devolvemos los matches que han pasado el filtro
  return final_matches


def getOrderedKeypoints(kp1, kp2, matcher):
  """Obtiene listas de puntos clave ordenadas por distancia.
  Argumentos posicionales:
  - kp1, kp2: Puntos clave
  - matcher: Correspondencias halladas por el matcher
  Devuelve:
  - Listas de puntos ordenadas por distancia"""

  # srcPoints y dstPoints son los vectores finales que devolveremos
  srcPoints = []
  dstPoints = []

  for m in sorted(matcher, key = lambda x: x.distance)[:100]:
    # Para cada pareja del matcher cogemos el correspondiente keypoint
    # del origen y el del destino, calculamos sus puntos con .pt y los
    # añadimos al vector final
    srcPoints.append(kp1[m.queryIdx].pt)
    dstPoints.append(kp2[m.trainIdx].pt)

  # Como queremos que estos vectores estén de una forma concreta para que
  # findHomography los vea como vectores de puntos, hacemos un reshape(-1,1,2)
  # Esto significa que se convertirá en una matriz con n filas (tantas como
  # keypoints tengamos), 1 columna en la que estarán todos los valores y
  # cada elemento es un vector de 2, es decir, un punto 2D
  srcPoints = np.reshape(srcPoints, (-1,1,2))
  dstPoints = np.reshape(dstPoints, (-1,1,2))
  #Devolvemos el resultado
  return srcPoints, dstPoints


if __name__ == "__main__":
  print("Este es un módulo auxiliar que no puede ejecutarse por sí solo.")
  print("Por favor, ejecuta 'python3 main.py'")
