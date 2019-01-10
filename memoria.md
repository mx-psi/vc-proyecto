# Introducción

TODO: Describir el problema en términos generales, sus aplicaciones, qué hemos hecho en el trabajo, etc.

En este proyecto seguiremos la notación empleada por el libro TODO.

# El problema de estimación de homografías

:::{.problem name="Estimación de homografías"}
Dado un conjunto de correspondencias $\{\mathbf{x}_i \leftrightarrow \mathbf{x}'_i\}_{i = 1, \dots, N}$ en coordenadas homogéneas,
hallar una homografía dada por su expresión matricial $H \in \mathcal{M}_{3 \times 3}(\mathbb{R})$ que minimice el error entre $H\mathbf{x}_i$ y $\mathbf{x}'_i$.
:::

Como utilizamos coordenadas homogéneas, el caso exacto se dará cuando $H\mathbf{x}_i$ y $\mathbf{x}'_i$ sean proporcionales, que podemos expresar como $\mathbf{x}'_i \times H\mathbf{x}_i = \mathbf{0}$.
Diremos que $\mathbf{x}_i \leftrightarrow \mathbf{x}'_i$ es una *correspondencia exacta* para $H$.

Además, en lo sucesivo utilizaremos $h$ para referirnos a un vector de $\mathbb{R}^9$ que tenga las componentes de $H$ escritas como vector.

- TODO: Describir el cálculo en el caso exacto para 4 correspondencias. 

La expresión matricial por bloques del sistema de ecuaciones que define una correspondencia, notando $\mathbf{x}_i = (x_i, y_i, w_i)$, $\mathbf{x}'_i = (x'_i, y'_i, w'_i)$ y $H = (h^1 h^2 h^3)^T$ por filas es (TODO citar 4.3 p.89)
$$A_iH =
\left(\begin{matrix}\mathbf{0}^T & -w'_i\mathbf{x}_i^T & y'_i\mathbf{x}_i^T \\
w'_i\mathbf{x}_i^T  & \mathbf{0}^T  & -x'_i\mathbf{x}_i^T \end{matrix}\right)
\left(\begin{matrix}h^1 \\ h^2 \\ h^3\end{matrix}\right)$$

Fijada la homografía podemos definir a partir de esta matriz la función de *error algebraico* $\mathcal{C}_H(\mathbf{X}) = A(\mathbf{X})h$.

# Error de Sampson
## Descripción y derivación

El error de Sampson es una función de coste aplicable a problemas generales de estimación en visión por computador.
En el caso de la estimación de homografías aproxima el *error geométrico*, esto es, a la distancia euclídea entre $H\mathbf{x}_i$ y $\mathbf{x}'_i$, dando mejores resultados que reducir en norma el error algebraico.
A continuación hacemos un breve resumen de la derivación teórica del error de Sampson para una sola correspondencia siguiendo el desarrollo y la notación de TODO.

La idea del error de Sampson es: dada una correspondencia expresada en coordenadas inhomogéneas como $\mathbf{X} = (x,y,x',y')$ estimar la correspondencia exacta para $H$ más cercana mediante una aproximación lineal de la función de error algebraico.

Aproximamos la función de error algebraico utilizando un desarrollo de Taylor de primer orden:
$$\mathcal{C}_H(\mathbf{X} + \mathbf{\delta_X}) = \mathcal{C}_H(\mathbf{X}) + \frac{\partial \mathcal{C}_H}{\partial \mathbf{X}}\mathbf{\delta_X}, \text{ donde notaremos } J \overset{\text{def}}{=} \frac{\partial \mathcal{C}_H}{\partial \mathbf{X}}, \mathbf{\epsilon} \overset{\text{def}}{=}  \mathcal{C}_H(\mathbf{X}).$$

Si $\mathbf{X} + \mathbf{\delta_X}$ es una correspondencia exacta para $H$ entonces tendremos que $\mathcal{C}_H(\mathbf{X} + \mathbf{\delta_X}) = \mathbf{0}$
$J\mathbf{\delta_X} = - \mathbf{\epsilon},$
luego encontrar la correspondencia exacta más cercana es encontrar el $\mathbf{\delta_X}$ más pequeño en norma que satisfaga $J\mathbf{\delta_X} = - \mathbf{\epsilon}$.

<!-- TODO: Podría describir los multiplicadores de Lagrange con más detalle -->

Aplicamos multiplicadores de Lagrange, para los que añadimos un vector $\mathbf{\lambda}$.
Tras derivar llegamos a las siguientes dos ecuaciones
$$JJ^T \lambda = -\varepsilon, \quad \mathbf{\delta_X} = J^T\mathbf{\lambda} \implies \mathbf{\delta_X} = -J^T(JJ^T)^{-1}\varepsilon$$

El error de Sampson para una correspondencia será entonces
$$\mathcal{S}_H(\mathbf{X}) \overset{\text{def}}{=} \lVert\mathbf{\delta_X}\rVert^2 = \epsilon^T(JJ^T)^{-1}\epsilon$$.

## Cálculo algorítmico

Para el cálculo algorítmico del error de Sampson para una correspondencia necesitamos una expresión de $J$ en términos de la correspondencia y la homografía.

Siguiendo la notación que utilizamos al plantear el error algebraico, si notamos cada fila de la matriz $A_i$ como $A^{(j)}$ y notamos $\mathbf{X} = (\mathbf{X}_1, \mathbf{X}_2, \mathbf{X}_3, \mathbf{X}_4) = (x,y,x',y')$ entonces tenemos que
$$J = \frac{\partial \mathcal{C}_H}{\mathbf{X}} \overset{\text{def}}{=} \left(\frac{\partial A^{(l)}}{\partial \mathbf{X}_j}\right)_{l,j} = \left
(\begin{matrix}
-w'_ih_{21} + y'_ih_{31} & -w'_ih_{22} + y'_ih_{32} & 0 &  xh_{31} + yh_{32} + h_{33} \\ 
w_ih_{11} - x'_ih_{31} & w_ih_{12}-x'_ih_{32} & xh_{31} + yh_{32} + h_{33} & 0
\end{matrix}\right)$$

Podemos comprobar sin más que sustituir ambos lados de la expresión que, al ser $\mathcal{C}_H$ multilineal en las coordenadas de $\mathbf{X}$, si $\{\mathbf{e}_i\}_{i = 1,\dots, 4}$ es la base usual de $\mathbb{R}^4$ entonces tenemos que 
$$\frac{\partial C_h(\mathbf{X})}{\partial \mathbf{X}_i} = \mathcal{C}_H(\mathbf{X} + \mathbf{e}_i) - \mathcal{C}_H(\mathbf{X})$$
y por tanto podemos calcular
$$JJ^T = \sum_{i = 1}^4 \left(\frac{\partial C_h(\mathbf{X})}{\partial \mathbf{X}_i}\right) \left(\frac{\partial C_h(\mathbf{X})}{\partial \mathbf{X}_i}\right)^T = \sum_{i = 1}^4 (\mathcal{C}_H(\mathbf{X} + \mathbf{e}_i) - \mathcal{C}_H(\mathbf{X})) (\mathcal{C}_H(\mathbf{X} + \mathbf{e}_i) - \mathcal{C}_H(\mathbf{X}))^T$$
lo que nos permite calcular $JJ^T$ directamente a partir de $\mathcal{C}_H$.

## Cálculo para múltiples correspondencias

En el algoritmo general tendremos múltiples correspondencias.
El error de Sampson para las múltiples correspondencias será la suma de los errores de Sampson individuales, esto es
$$\mathcal{S}_H(\{\mathbf{x}_i \leftrightarrow \mathbf{x}'_i\}_{i = 1, \dots, N}) = \sum_{i = 1}^N \mathcal{S}_H(\mathbf{X}_i)$$

# Algoritmo de estimación de homografías

En esta sección describimos el algoritmo iterativo de estimación de homografías 2D a partir de un conjunto de correspondencias haciendo uso del error de Sampson.

## Estimación inicial

El algoritmo iterativo necesita una estimación inicial a partir de la cuál iterar hasta converger en un mínimo.
En principio podría utilizarse una estimación inicial que no dependiera del conjunto de correspondencias.
Sin embargo, en la práctica, el uso de este tipo de estimaciones iniciales no proporciona buenos resultados, ya que el algoritmo tiende a converger a mínimos locales lejos del óptimo.

Por tanto, para la estimación inicial seguiremos la recomendación de [TODO referencia] para utilizar el método *DLT* con normalización, descrito por el siguiente algoritmo.

:::{.algorithm name="\textit{Direct Linear Transform} normalizado"}

$\;$

Entrada
: Un conjunto de correspondencias $\{\mathbf{x}_i \leftrightarrow \mathbf{x}'_i\}_{i = 1, \dots, N}$ en coordenadas homogéneas.

Salida
:  Homografía dada por su expresión matricial $H \in \mathcal{M}_{3 \times 3}(\mathbb{R})$ que se ajusta a las correspondencias.

1. Normalizar $\mathcal{X} = \{\mathbf{x}_i\}_{i = 1,\dots, N}$ mediante una transformación $T$ compuesta de
   a) una traslación que lleve el centroide de $\mathcal{X}$ al origen y
   b) un escalado que lleve la distancia media al origen a $\sqrt{2}$.
2. Normalizar $\mathcal{X}' = \{\mathbf{x}'_i\}_{i = 1,\dots, N}$ mediante una transformación $T'$ compuesta de
   a) una traslación que lleve el centroide de $\mathcal{X}'$ al origen y
   b) un escalado que lleve la distancia media al origen a $\sqrt{2}$.
3. Para cada correspondencia $\mathbf{x}_i \leftrightarrow \mathbf{x}'_i$ calcular la matriz $A_i$,
4. Concatenar las matrices $A_i$ en una única matriz $A$,
5. Calcular el vector singular asociado al menor valor singular de la descomposición en valores de $A$.
6. Transformar ese vector en una matriz $\overset{\sim}{H}$.
7. Devolver la composición de las matrices de transformación con la homografía estimada $T'^{-1}\overset{\sim}{H}T$.
:::

