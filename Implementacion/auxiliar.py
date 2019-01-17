#! /usr/bin/python3
# Authors: Pablo Baeyens y Antonio Checa
# auxiliar.py
# Fichero auxiliar para mostrar un ejemplo con correspondencias entre imágenes
# Adaptado de la práctica 2 de Antonio


import numpy as np
import cv2
import math

#Esta funcion, dado un nombre de fichero fileName y un bool lee la imagen en
# color si flagColor = 1 y en gris si flagColor=0, luego la devuelve
def lee_imagen(fileName, flagColor):
	img = cv2.imread(fileName,flagColor)
	return img

#Dada una imagen im y un string nombre muestra por pantalla la imagen
def pintaI(im, nombre):
	cv2.namedWindow(nombre, cv2.WINDOW_NORMAL)
	cv2.resizeWindow(nombre, im.shape[1], im.shape[0]) # Hacemos resize de la
								#ventana para que la imagen se vea en su tamano
	cv2.imshow(nombre, im)
	cv2.waitKey(0)
	cv2.destroyAllWindows()
	return im

# Función que dada una imagen, un objeto SIFT y unos keypoints devuelve los descriptores
def getDescriptors(im, keypoints, sift):
	kp, des = sift.compute(im, keypoints)
	return des

# Función que dada una imagen hace la detección de keypoints con detect de SIFT
# y de forma separada calcula los descriptores con compute (llamando a getDescriptors)
def getKpAndDescriptors(im):
	sift = cv2.xfeatures2d.SIFT_create(1000, 3, 0.08, 8, 1.6)
	kp = sift.detect(im)
	des = getDescriptors(im, kp, sift)
	return (kp, des)

##Ejercicio 2

# Función que realiza un match con el ratio test de Lowe y un knnMatch con k=2
# Devuelve los matches que pasan este filtrado
def getMatchesLowe2NN(des1, des2):
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


# Función auxiliar que nos da las listas de keypoints kp1 y kp2 ordenadas
# por distancia en el matcher. Se considera que kp1 es de la imagen origen
# mientras que kp2 es de la imagen destino del matcher
def getOrderedKeypoints(kp1, kp2, matcher):
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
	return (srcPoints, dstPoints)
