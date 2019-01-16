#!python3

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

#Dado un vector de imagenes vim y un string nombre, realiza la imagen
# concatenando las diferentes imagenes, lo guarda en full_im, muestra esa
# imagen y la devuelve
def pintaMI(vim, nombre):
	if len(vim[0].shape) < 3:
		height, width = vim[0].shape
		sca = 0 # La variable sca nos guarda el valor escalar con el
				# que rellenaremos los bordes de las imagenes chicas
	else:
		height, width, channel = vim[0].shape
		sca = (0,0,0) # Si la imagen esta en color, el valor escalar debe ser
					  # una tupla de tres valores

	# Conseguimos la maxima altura y anchura de las imagenes para saber de que
	# longitudes finales va a ser la imagen. Lo guardamos en mx_width y mx_height
	mx_width = 0
	mx_height = 0
	for im in vim:
		if mx_width < im.shape[1]:
			mx_width = im.shape[1]
		if mx_height < im.shape[0]:
			mx_height = im.shape[0]

	# Hacemos un nuevo vector que va a ir guardando las imagenes transformadas de
	# forma que todas tengan misma anchura, altura y canales
	nuevo_vim = []
	for im in vim:
		# Cogemos la diferencia de anchura y altura
		dif_wid = mx_width-im.shape[1]
		dif_hei = mx_height-im.shape[0]
		# Rellenamos los bordes, cuidando de rellenar exactamente la diferencia
		# entre esta imagen y el maximo.
		nueva_img = cv2.copyMakeBorder(im, dif_hei//2, (dif_hei+1)//2, dif_wid//2, (dif_wid+1)//2, cv2.BORDER_CONSTANT, sca)
		#Si la imagen tiene un solo canal, rellenamos la imagen con dos
		# canales mas, repitiendo el que hay. De esa forma seguira viendose en
		# gris pero tendra tres canales. Usamos la funcion stack para hacer la
		# traspuesta con axis = -1
		if len(im.shape) < 3:
			nueva_img = np.stack((nueva_img,)*3,-1)
		nuevo_vim.append(nueva_img)

	#Con hconcat metemos todas las imagenes del nuevo vector en una sola
	full_im = []
	full_im = cv2.hconcat(nuevo_vim, len(vim))

	# Mostramos la imagen y la devolvemos
	pintaI(full_im, nombre)
	return full_im

# Función que varia umbrales en el descriptor SIFT de la imagen orig_im
# Pinta los resultados y escribe por pantalla el número de keypoints por cada valor
def siftEj1(orig_im):
	# Variamos el umbral de contraste, contrastThreshold, que se pone en el tercer parámetro
	vec_sift_umbral = []
	# Hacemos un vector con los valores que vamos a probar
	umbral = [0.005, 0.05, 0.08, 0.1, 0.12]
	for i in umbral:
		im = orig_im.copy()
		#El resto de parámetros los dejamos a valores fijos
		sift = cv2.xfeatures2d.SIFT_create(0, 3, i, 10, 1.6)
		keypoints, des = sift.detectAndCompute(im, None)
		#Ponemos por pantalla el número de keypoints por cada valor
		print("Keypoints con el valor del umbral de contraste " + str(i) + ": " + str(len(keypoints)))
		outI = im
		#Creamos una imagen con drawKeypoints
		img=cv2.drawKeypoints(im, keypoints, outI, flags=cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
		vec_sift_umbral.append(img)
	#Pintamos todas las imágenes juntas
	p = pintaMI(vec_sift_umbral, "ejercicio 1 umbral")
	#cv2.imwrite("./resultados/ej1-aContrast.png", p)
	input("Pulse ENTER para continuar.")

	# Variamos el umbral de bordes, edgeThreshold, que se pone en el cuarto parámetro
	vec_sift_edge = []
	# Hacemos un vector con los valores que vamos a probar
	edges_thres = [2, 4, 6, 8, 10]
	for i in edges_thres:
		im = orig_im.copy()
		#El resto de parámetros se deja a un valor fijo
		sift = cv2.xfeatures2d.SIFT_create(0, 3, 0.08, i, 1.6)
		keypoints, des = sift.detectAndCompute(im, None)
		#Ponemos por pantalla el número de keypoints por cada valor
		print("Keypoints con el valor del umbral de bordes " + str(i) + ": " + str(len(keypoints)))
		outI = im
		#Creamos una imagen con drawKeypoints
		img=cv2.drawKeypoints(im, keypoints, outI, flags=cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
		vec_sift_edge.append(img)
	#Pintamos todas las imágenes juntas
	p = pintaMI(vec_sift_edge, "ejercicio 1 edge threshold")
	#cv2.imwrite("./resultados/ej1-aEdges.png", p)
	input("Pulse ENTER para continuar.")

	#Hacemos una imagen con parámetros equilibrados y 1000 keypoints finales
	sift = cv2.xfeatures2d.SIFT_create(1000, 3, 0.08, 8, 1.6)
	keypoints, des = sift.detectAndCompute(im, None)
	p = pintaI(cv2.drawKeypoints(im, keypoints, outI, flags=cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS), "Ejercicio 1a final")
	#cv2.imwrite("./resultados/ej1-afinal.png", p)

#Función que varía el umbral de la Hessiana en el descriptor SURF con la imagen orig_im
# pinta los resultados y escribe por pantalla el número de keypoints por valor
def surfEj1(orig_im):
	#Ponemos en un vector los diferentes valores del umbral de la Hessiana que
	# vamos a probar
	umbral = [400, 600, 800, 1000, 1200]
	vec_surf_umbral = []
	for u in umbral:
		im = orig_im.copy()
		#Creamos un objeto SURF con el umbral deseado y el resto de parámetros por defecto
		s = cv2.xfeatures2d.SURF_create(u, 4, 3, False, False)
		kp, des = s.detectAndCompute(im,None)
		#Ponemos por pantalla el número de keypoints por valor
		print("Keypoints con el valor del umbral de la hessiana " + str(u) + ": " + str(len(kp)))
		outI = im.copy()
		#Disminuimos el tamaño del size de los keypoints para que se pinten más chicos
		for k in kp:
			k.size = k.size/4
		im2 = cv2.drawKeypoints(im,kp, outI, flags=cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
		vec_surf_umbral.append(im2)
	# Pintamos todas las imágenes juntas para poder compararlas
	p = pintaMI(vec_surf_umbral, "ejercicio 1 surf umbral")
	input("Pulse ENTER para continuar.")
	#cv2.imwrite("./resultados/ej1-asurf.png", p)

# Función que llama a siftEj1 y surfEj1 para poder hacer un análisis de ambos métodos
def ejercicio1a(im):
	siftEj1(im)
	surfEj1(im)

#Función de la documentación de OpenCV para SIFT en la que se coge la variable
# octave del keypoint y se saca la octava real
def unpackOctave(kpt):
	octave = kpt.octave & 255
	if(octave > 128):
		octave = (-128 | octave)
	return octave


#Función de la documentación de OpenCV para SIFT en la que se coge la variable
# octave del keypoint y se saca la escala
def unpackScale(kpt):
	octave = unpackOctave(kpt)
	if(octave >= 0):
		scale = 1.0/(1 << octave)
	else:
		scale = (float)(1 << -octave)
	return scale


#Función de la documentación de OpenCV para SIFT en la que se coge la variable
# octave del keypoint y se saca la capa
def unpackLayer(kpt):
	return (kpt.octave >> 8) & 255

# Función que devuelve n colores "diferentes", intenta distribuirlos equiespaciados
def getNColors(n):
	colors_vec = []
	n_colors = (n//6+1)*6
	for j in range(n_colors):
		i = j+1
		# Si nos faltan colores, añadimos seis nuevos
		if(len(colors_vec) < n):
			colors_vec.append((255//i, 0, 0))
			colors_vec.append((255//i, 255//i, 0))
			colors_vec.append((255//i, 0, 255//i))
			colors_vec.append((0, 255//i, 0))
			colors_vec.append((0, 255//i, 255//i))
			colors_vec.append((0, 0, 255//i))
	return colors_vec


#Función que se encarga de llamar a SIFT, ver cuántos keypoints caen en cada capa,
# en cada octava y poner estos resultados por pantalla. Además, pinta una imagen
# en la que los colores de los círculos de cada keypoint varían en función de su
# octava
def ejercicio1bSift(im):
	sift = cv2.xfeatures2d.SIFT_create(1000, 3, 0.08, 8, 1.6)
	keypoints, des = sift.detectAndCompute(im, None)
	# Vamos a guardar en un diccionario los puntos asociados a cada octava y a cada capa
	# La llave será la octava o la capa correspondiente y el valor será el vector
	# de puntos con sus tamaños
	dict_octaves = {}
	dict_layer = {}
	for k in keypoints:
		o = unpackOctave(k)
		l = unpackLayer(k)
		# Rellenamos el diccionario con la octava actual
		# Si está, añadimos el punto, si no creamos la entrada con esa key
		if(o in dict_octaves.keys()):
			dict_octaves[o].append((k.pt, k.size))
		else:
			dict_octaves[o] = [(k.pt, k.size)]
		# Rellenamos el diccionario con la capa actual
		if(l in dict_layer.keys()):
			dict_layer[l].append(k.pt)
		else:
			dict_layer[l] = [k.pt]
	# Cogemos tantos colores como octavas tengamos
	colors = getNColors(len(dict_octaves.keys()))
	im2 = im.copy()
	#Por cada octava hacemos el análisis
	for o in dict_octaves.keys():
		#Calculamos el vector de puntos y tamaños
		v_points = dict_octaves[o]
		#Ponemos por pantalla cuántos keypoints hay en esta octava y los colores
		# Hay que tener cuidado, ya que los colores se poner por pantalla en formato
		# BGR, que es como OpenCV los interpreta
		print("Se han detectado un total de " + str(len(v_points)) + " en la octava " + str(o) +" de color (" + str(colors[o+1]))
		# Pintamos el círculo de cada punto dentro del vector v_points
		# El color que cogemos depende solo de la octava en la que estamos
		for point, size_pt in v_points:
			point = (int(point[0]), int(point[1]))
			im2 = cv2.circle(im2, point, int(size_pt), colors[o+1])
	# Pintamos la imagen resultado
	p = pintaI(im2, "ejercicio1b Sift")
	#cv2.imwrite("./resultados/ej1-bsift.png", p)
	#Ponemos una parada para leer los resultados anteriores
	input("Pulse ENTER para continuar.")

	# Ponemos por pantalla cuántos keypoints por capa tenemos
	for l in dict_layer.keys():
		v_points = dict_layer[l]
		print("Se han detectado un total de " + str(len(v_points)) + " en la capa " + str(l))
	input("Pulse ENTER para continuar")

#Función que se encarga de llamar a SURF, ver cuántos keypoints caen en cada octava
# poner estos resultados por pantalla. Además, pinta una imagen en la que los
# colores de los círculos de cada keypoint varían en función de su octava
def ejercicio1bSurf(im):
	# Ponemos como umbral de la hessiana 800, que era un valor aceptable del 1a
	surf = cv2.xfeatures2d.SURF_create(hessianThreshold = 800)
	kp, des = surf.detectAndCompute(im, None)
	# Esta vez solo necesitaremos un diccionario para la octava
	dict_octaves = {}
	for k in kp:
		# Solo tenemos que acceder a k.octave para conseguir la octava
		o = k.octave
		# Agregamos el punto al diccionario
		if(k.octave in dict_octaves.keys()):
			dict_octaves[o].append((k.pt,k.size))
		else:
			dict_octaves[o] = [(k.pt, k.size)]
	# Cogemos tantos colores como octavas tengamos
	colors = getNColors(len(dict_octaves.keys()))
	im2 = im.copy()
	for o in dict_octaves.keys():
		# Cogemos el vector de puntos y tamaños de la octava
		v_points = dict_octaves[o]
		# Ponemos por pantalla cuántos keypoints hay en esta octava y los colores
		# Hay que tener cuidado, ya que los colores se poner por pantalla en formato
		# BGR, que es como OpenCV los interpreta
		print("Se han detectado un total de " + str(len(v_points)) + " en la octava " + str(o)  +" de color (" + str(colors[o]))
		# Pintamos el círculo de cada keypoint con el tamaño proporcional al k.size
		# pero disminuido, para que la imagen aparezca mejor
		# El color solo depende de la octava en la que estamos
		for point, size_pt in v_points:
			point = (int(point[0]), int(point[1]))
			im2 = cv2.circle(im2, point, int(size_pt)//4, colors[o])

	# Pintamos la imagen final
	p = pintaI(im2, "ejercicio1b Surf")
	#cv2.imwrite("./resultados/ej1-bsurf.png", p)
	# Ponemos una parada para leer los datos de los keypoints
	input("Pulse ENTER para continuar.")

#Función del ejercicio 1b) que solo llama a la versión SIFT y SURF anteriores
# para pintar los keypoints de forma diferente según en qué octava estén y detectar
# poniendo por pantalla cuántos caen en cada una. En SIFT, además, se añaden
# cuántos keypoints por capa caen.
def ejercicio1b(im):
	ejercicio1bSift(im)
	ejercicio1bSurf(im)

# Función que dada una imagen, un objeto SIFT y unos keypoints devuelve los descriptores
def getDescriptors(im, keypoints, sift):
	kp, des = sift.compute(im, keypoints)
	return des

# Función que dada una imagen hace la detección de keypoints con detect de SIFT
# y de forma separada calcula los descriptores con compute (llamando a getDescriptors)
def ejercicio1cSift(im):
	sift = cv2.xfeatures2d.SIFT_create(1000, 3, 0.08, 8, 1.6)
	kp = sift.detect(im)
	des = getDescriptors(im, kp, sift)
	return (kp, des)

# Función que dada una imagen hace la detección de keypoints con detect de SURF
# y de forma separada calcula los descriptores con compute
def ejercicio1cSurf(im):
	surf = cv2.xfeatures2d.SURF_create(hessianThreshold=800)
	kp = surf.detect(im)
	kp, des = surf.compute(im, kp)
	return (kp, des)

# Función del ejercicio 1c) que solo llama a la versión SIFT y SURF para calcular
# con ambos métodos los keypoints y descriptores por separado
def ejercicio1c(im):
	ejercicio1cSift(im)
	ejercicio1cSurf(im)

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

# Función que mediante los descriptores y keypoints de SIFT de dos imágenes
# establece correspondencias o matches entre ellos con dos métodos: fuerza bruta
# + crossCheck y otro que es el de Lowe con un 2nn. Muestra por pantalla las
# imágenes resultantes de ambos métodos, con líneas trazadas entre las correspondencias
def ejercicio2(im1, im2):
	# Creamos los keypoints y descriptores con el 1cSIFT
	kp1, des1 = ejercicio1cSift(im1)
	kp2, des2 = ejercicio1cSift(im2)

	# Para el método fuerza bruta + cross Check solo creamos un BFMatcher con
	# crossCheck = true
	matcher = cv2.BFMatcher(crossCheck = True)
	# Hacemos el match entre los descriptores, lo que nos genera las correspondencias
	matches = matcher.match(des1, des2)
	outI = im1.copy()
	# Pintamos las primeras 100 correspondencias en pantalla con drawMatches
	# Se considera que esto coge 100 al azar porque en principio no están ordenadas
	p = pintaI(cv2.drawMatches(im1, kp1, im2, kp2, matches[:100], outI), "ejercicio 2 BruteForce + CrossCheck")
	#cv2.imwrite("./resultados/ej2BF.png", p)

	# Para el segundo método llamamos a la función auxiliar getMatchesLowe2NN
	# que nos realizaba el matcher con el ratio test de Lowe y 2nn
	final_matches = getMatchesLowe2NN(des1, des2)
	outI = im1.copy()
	# Pintamos las primeras 100 correspondencias en pantalla con drawMatches
	p = pintaI(cv2.drawMatches(im1, kp1, im2, kp2, final_matches[:100], outI), "ejercicio2 2nn")
	#cv2.imwrite("./resultados/ej22nn.png", p)

# Función auxiliar que nos da las listas de keypoints kp1 y kp2 ordenadas
# según salgan en el matcher. Se considera que kp1 es de la imagen origen
# mientras que kp2 es de la imagen destino del matcher
def getOrderedKeypoints(kp1, kp2, matcher):
	# srcPoints y dstPoints son los vectores finales que devolveremos
	srcPoints = []
	dstPoints = []

	for m in matcher:
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

# Función que implementa el ejercicio 3 del guion entero. Genera el mosaico de
# calidad con las imagenes im1, im2 e im3 y lo pone por pantalla. No se ha
# separado en tres apartados porque se ha creído conveniente tener todo junto
# y los apartados del guion no tenían mucho sentido de ser separados, se ha
# entendido que solo era una lista de tareas a realizar
def ejercicio3(im1, im2, im3):
	# Cogemos los descriptores de las tres imágenes y sus keypoints
	kp1, des1 = ejercicio1cSift(im1)
	kp2, des2 = ejercicio1cSift(im2)
	kp3, des3 = ejercicio1cSift(im3)

	# Cogemos los matchers que van en dirección a la imagen central im2
	# Esto se realiza de esta forma para evitar acumulación de errores
	# Los matchers los cogemos con Lowe2nn que era el que tenía más calidad del
	# ejercicio 2
	matcher12 = getMatchesLowe2NN(des1, des2)
	matcher32 = getMatchesLowe2NN(des3, des2)

	# Recogemos las listas de keypoints de cada matcher ordenadas según las
	# correspondencias
	orderSrcKp1, orderDstKp12 = getOrderedKeypoints(kp1, kp2, matcher12)
	orderSrcKp3, orderDstKp32 = getOrderedKeypoints(kp3, kp2, matcher32)

	#(s_x, s_y) es el tamaño final de nuestro canvas
	# Estos valores se han cogido así porque son los que van bien con las imágenes
	# de yosemite
	s_x = 1400
	s_y = 600
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
	# Como va a ser exacta, de cuatro puntos en cuatro puntos, no necesitamos
	# aclarar que utilice cv2.RANSAC
	h_canvas, mask = cv2.findHomography(np.array([[0,0], [wi, 0], [0,hi], [wi, hi]]), np.array([[p_x, p_y], [p_x+wi, p_y], [p_x, hi+p_y], [wi+p_x, hi+p_y]]))

	# Hallamos las otras dos homografías, en dirección a la imagen central
	# Esta vez se aclara que se utilice RANSAC con un máximo error de proyección
	# de 1 para los inliners
	h1, mask = cv2.findHomography(orderSrcKp1, orderDstKp12, cv2.RANSAC, ransacReprojThreshold=1)
	h2, mask = cv2.findHomography(orderSrcKp3, orderDstKp32, cv2.RANSAC, ransacReprojThreshold=1)

	# Se crea la imagen final haciendo llamadas a warpPerspective de cada imagen
	# con sus transformaciones correspondientes
	# Mientras que la central solo es la homografía que hemos hallado antes,
	# para las de los extremos necesitamos componer las homografías que las llevan
	# a la central con la que la lleva a canvas. Como lo primero que se aplica son
	# las que la llevan a la central, el orden al componer es h_canvas * h_12
	# y h_canvas * h_32. Este producto se hace con dot
	canvas_final2 = canvas_final
	# El borderMode del primer warpPerspective se pone a Constant para el fondo
	# El resto de borderMode se ponen a TRANSPARENT para no pisar el resto de imágenes
	canvas_final2 = cv2.warpPerspective(im1,  h_canvas.dot(h1), (s_x, s_y), canvas_final, borderMode = cv2.BORDER_CONSTANT)
	canvas_final2 = cv2.warpPerspective(im2, h_canvas, (s_x, s_y), canvas_final2, borderMode = cv2.BORDER_TRANSPARENT)
	canvas_final2 = cv2.warpPerspective(im3, h_canvas.dot(h2), (s_x, s_y), canvas_final2, borderMode = cv2.BORDER_TRANSPARENT, borderValue = 255)

	#Ponemos por pantalla el mosaico final
	p = pintaI(canvas_final2, "ejercicio 3")
	#cv2.imwrite("./resultados/ej3.png", p)


# Función que implementa el mosaico con n imágenes que se pide en el ejercicio 4
# Es fundamentalmente igual que el 3 pero pasando a vectores cada dato. Pone por
# pantalla el mosaico que se produce de juntar las imágenes en vec_im
def ejercicio4(vec_im):
	# Hallamos los keypoints y descriptores de nuestras imágenes
	vec_kp = []
	vec_des = []
	for i in vec_im:
		kp, des = ejercicio1cSift(i)
		vec_kp.append(kp)
		vec_des.append(des)

	# Conseguimos la imagen central que tiene índice st
	n = len(vec_im)
	st = n//2
	# Hallamos los matches del resto de imágenes a la que tienen en dirección
	# a la imagen central. Todos ellos se guardan en matches, y se hacen con
	# Lowe y 2nn
	matches = []
	for i in range(st):
		# En la primera mitad se va de la i a la i+1
		matches.append(getMatchesLowe2NN(vec_des[i], vec_des[i+1]))
	for i in range(st, n-1):
		# En la segunda mitad se va de la i+1 a la i
		matches.append(getMatchesLowe2NN(vec_des[i+1], vec_des[i]))

	# Se calculan los vectores de keypoints ordenados según digan las correspondencias
	srcKeyp = []
	dstKeyp = []
	for i in range(st):
		# En la primera mitad se va de la i a la i+1
		src, dst = getOrderedKeypoints(vec_kp[i], vec_kp[i+1], matches[i])
		srcKeyp.append(src)
		dstKeyp.append(dst)
	for i in range(st, n-1):
		# En la segunda mitad se va de la i+1 a la i
		src, dst = getOrderedKeypoints(vec_kp[i+1], vec_kp[i], matches[i])
		srcKeyp.append(src)
		dstKeyp.append(dst)

	# Tamaño de nuestra imagen final. 1000x500 es un buen valor para las de mosaico
	s_x = 1000
	s_y = 500

	# Nuestra imagen final es canvas_final
	canvas_final = np.zeros((s_x, s_y, 3), np.uint8)

	# hi, wi es el tamaño de la imagen central
	hi = vec_im[st].shape[0]
	wi = vec_im[st].shape[1]
	# (p_x, p_y) es el punto donde acabará trasladado el (0,0) de la imagen central
	p_x = 400
	p_y = 50

	# Hallamos la homografía de la imagen central, que será solo una traslación
	# así que no añadimos cv2.RANSAC
	h_canvas, mask = cv2.findHomography(np.array([[0,0], [wi, 0], [0,hi], [wi, hi]]), np.array([[p_x, p_y], [p_x+wi, p_y], [p_x, hi+p_y], [wi+p_x, hi+p_y]]))


	# Guardamos en homographies el vector de las homografías que necesitaremos
	homographies = []
	for i in range(n-1):
		# Hallamos las n-1 homografías
		h, mask = cv2.findHomography(srcKeyp[i], dstKeyp[i], cv2.RANSAC, 1)
		homographies.append(h)

	# Llevamos la imagen central al canvas
	canvas_final = cv2.warpPerspective(vec_im[st], h_canvas, (s_x, s_y), borderMode = cv2.BORDER_CONSTANT)

	# h será en el siguiente bucle la transformación que llevaremos hasta el momento
	# así que la iniciamos a la homografía inicial de la imagen central
	h = h_canvas
	# Recorremos la primera mitad del vector al revés, yendo de la imagen central
	# a la primera imagen, para pasarlas al canvas
	for i in range(st-1, -1, -1):
		# En cada iteración componemos con la homografía de la imagen iésima
		h = h.dot(homographies[i])
		# Pintamos la imagen iésima en el canvas con h, que contiene todas
		# las homografías que hacen falta para llevarla a la imagen central
		# y luego al canvas. Ponemos cv2.RANSAC y 1 como máximo error de proyección
		# para los inliners
		canvas_final = cv2.warpPerspective(vec_im[i], h, (s_x, s_y), canvas_final, borderMode = cv2.BORDER_TRANSPARENT)

	# Hacemos el bucle análogo para la segunda mitad, iniciamos h a h_canvas
	h = h_canvas
	for i in range(st, n-1):
		# En cada iteración componemos con la homografía iésima
		h = h.dot(homographies[i])
		# Pasamos la imagen (i+1)-ésima al canvas final
		canvas_final = cv2.warpPerspective(vec_im[i+1], h, (s_x, s_y), canvas_final, borderMode = cv2.BORDER_TRANSPARENT)
	# Pintamos el mosaico final
	p = pintaI(canvas_final, "ejercicio 4, mosaico")
	#cv2.imwrite("./resultados/ej4.png", p)
