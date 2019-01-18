\newpage

# Introducción

En este proyecto desarrollamos el problema de la estimación de una homografía dado un conjunto de correspondencias usando la función de error de Sampson. Una homografía queda determinada solo por cuatro puntos no alineados y sus cuatro imágenes, sin embargo, cuando tenemos incertidumbre respecto a las correspondencias podemos tener más de 4 correspondencias para las que no existirá una homografía que las refleje de forma exacta.

Para elegir qué homografía usaremos en estos casos una de las técnicas fundamentales consiste en disminuir el error que obtenemos al comparar las imágenes esperadas de los puntos con las que nos da la homografía. El error conceptualmente más simple es el error geométrico, ya que lo único que tenemos que hacer es comparar las distancias entre las imágenes de la homografía y los datos de destino. De esta idea es de donde sale el estudio de la estimación de homografías a partir de una función de coste basada en una distancia geométrica, que se explica con detalle en [@hartley2003multiple] y sobre el que basaremos la notación a lo largo del proyecto.

En este trabajo vamos a intentar explorar la estimación de una homografía con una función de coste muy concreta, el *error de Sampson*. Esta función de coste es una estimación de la distancia geométrica de la que hablábamos antes, cuya expresión es compleja para su optimización, pero que puede aproximarse de forma efectiva. El error de Sampson nos permite justo esto, encontrar una aproximación eficiente del error geométrico que producimos al realizar la homografía. Desarrollamos esta idea en la sección [Error de Sampson].

Las posibles aplicaciones de este algoritmo son las de mejorar el proceso de estimación de homografías que estudiamos en la segunda práctica. Por un lado, esto añade calidad a la homografía resultante, ya que al reducir el error geométrico de los puntos con sus imágenes por la homografía, obtenemos un resultado visualmente mejor. Por otro, añade consistencia en los casos en los que las correspondencias tengan más ruido, permitiendo resolver una gama más amplia de problemas en los que antes el algoritmo no tenía la suficiente calidad como para ser aceptable. Por último, cabe destacar que las herramientas necesarias para resolver un problema tan abstracto como es el de estimar una homografía bajo ciertas restricciones son importantes por sí mismas, por la utilidad que puedan dar en otros ámbitos y problemas.

El trabajo se divide en el código adjunto en el que implementamos todos estos algoritmos con ejemplos de uso para que sean más fáciles de entender y cuyo funcionamiento se explica en el [Apéndice: Funcionamiento del código adjunto], y en esta memoria, en la que explicamos el fundamento teórico de los mismos. En las siguientes secciones veremos: una introducción al problema a resolver; la descripción y obtención del Error de Sampson junto a su importancia y, por último, la explicación del algoritmo final, que se puede dividir en estimación inicial y el paso iterativo para ir minimizando el error de Sampson. La estimación inicial se hace con el algoritmo *Direct Linear Transform* (DLT) y el método iterativo usa el algoritmo de Levenberg-Marquadt.

# El problema de estimación de homografías

:::{.problem name="Estimación de homografías"}
Dado un conjunto de correspondencias $\{\mathbf{x}_i \leftrightarrow \mathbf{x}'_i\}_{i = 1, \dots, N}$ en coordenadas homogéneas,
hallar una homografía dada por su expresión matricial $H \in \mathcal{M}_{3 \times 3}(\mathbb{R})$ que minimice el error entre $H\mathbf{x}_i$ y $\mathbf{x}'_i$.
:::

Como utilizamos coordenadas homogéneas, el caso exacto se dará cuando $H\mathbf{x}_i$ y $\mathbf{x}'_i$ sean proporcionales, ya que la expresión matricial de una homografía está definida salvo multiplicación por constante.

Podemos expresar esta relación como $\mathbf{x}'_i \times H\mathbf{x}_i = \mathbf{0}$, sabiendo que el producto vectorial de vectores proporcionales es nulo. En este caso diremos que $\mathbf{x}_i \leftrightarrow \mathbf{x}'_i$ es una *correspondencia exacta* para $H$.

Además, en lo sucesivo utilizaremos $h$ para referirnos a un vector de $\mathbb{R}^9$ que tenga las componentes de $H$ escritas como vector.

Cada correspondencia impone tres ecuaciones sobre las componentes de la homografía, de las cuales dos son independientes. En concreto, la expresión matricial por bloques del sistema de ecuaciones que define una correspondencia, notando $\mathbf{x}_i = (x_i, y_i, w_i)$, $\mathbf{x}'_i = (x'_i, y'_i, w'_i)$ y $H = (h^1 h^2 h^3)^T$ por filas es ([@hartley2003multiple 4.3 p.89])
$$A(\mathbf{X})h =  A_ih =
\left(\begin{matrix}\mathbf{0}^T & -w'_i\mathbf{x}_i^T & y'_i\mathbf{x}_i^T \\
w'_i\mathbf{x}_i^T  & \mathbf{0}^T  & -x'_i\mathbf{x}_i^T \end{matrix}\right)
\left(\begin{matrix}h^1 \\ h^2 \\ h^3\end{matrix}\right)$$

Si tenemos 4 correspondencias podemos entonces formar a partir de las matrices $A_i \in \mathcal{M}_{2\times 9}(\mathbb{R})$ una matriz $A \in \mathcal{M}_{8\times 9}(\mathbb{R})$ que, junto con una condición extra como $\lVert h\rVert = 1$ nos da un sistema de ecuaciones que nos permite obtener una homografía para la cual las 4 correspondencias sean exactas.

Esta resolución sin embargo no es útil en la práctica, ya que normalmente las correspondencias tienen ruido y *outliers*, por lo que la homografía obtenida no será útil en la práctica. Veremos un método sencillo para obtener a partir de $N$ correspondencias, una homografía usando el método DLT que nos será útil como estimación inicial para posteriormente refinar con el error de Sampson.

Fijada una homografía, para definir el error de Sampson  Fijada la homografía podemos definir a partir de esta matriz la función de *error algebraico* para una correspondencia $\mathbf{X}$, $\mathcal{C}_H(\mathbf{X}) = A(\mathbf{X})h$.

# Error de Sampson
## Descripción y derivación

El error de Sampson es una función de coste aplicable a problemas generales de estimación en visión por computador.
En el caso de la estimación de homografías aproxima el *error geométrico*, esto es, a la distancia euclídea entre $H\mathbf{x}_i$ y $\mathbf{x}'_i$, dando mejores resultados que reducir en norma el error algebraico.
A continuación hacemos un breve resumen de la derivación teórica del error de Sampson para una sola correspondencia siguiendo el desarrollo y la notación de [@hartley2003multiple].

La idea principal detrás del error de Sampson es: dada una correspondencia expresada en coordenadas inhomogéneas como $\mathbf{X} = (x,y,x',y')$ estimar la correspondencia exacta para $H$ más cercana mediante una aproximación lineal de la función de error algebraico.

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

Por tanto, para la estimación inicial seguiremos la recomendación de [@hartley2003multiple] para utilizar el método *DLT* con normalización, descrito por el siguiente algoritmo.

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

Este algoritmo nos da una buena aproximación inicial.
La normalización nos permite evitar errores derivados de la precisión limitada en el cálculo y mejorar los resultados.
El cálculo del vector singular es el algoritmo utilizado normalmente para la resolución de sistemas de ecuaciones lineales en los que hay más ecuaciones que incógnitas (asumimos por tanto que $N >4$).

## Método iterativo

El método iterativo que se ha usado es, como ya se ha mencionado previamente, el algoritmo de Levenberg-Marquadt. Este algoritmo es una variación de Gauss-Newton en la que se cambia la ecuación de iteración, que, en función de la convergencia del algoritmo oscila entre el método de Gauss-Newton y el método de gradiente descendente.

La ecuación del método de Gauss Newton es
$$J^TJ\Delta = -J^T\epsilon,$$
mientras que la que se utiliza en el método iterativo de Levenberg-Marquadt es
$$(J^TJ+\lambda I)\Delta = -J^T\epsilon$$
donde $\lambda$ es un término que va cambiando por iteración: aumenta si no conseguimos disminuir el error y disminuye si sí lo hacemos. Esta ecuación se llama la ecuación normal aumentada.


La importancia de que $\lambda$ cambie no debe de pasar por alto, ya que es el núcleo de este método iterativo. Imaginemos qué ocurre en las dos posibilidades que hay.

\begin{itemize}
  \item Por un lado, si $\lambda$ se acerca a 0, como ya hemos mencionado antes, la ecuación nos queda exactamente la misma que en Gauss-Newton. En un caso teórico ideal, estaremos en esta situación conforme mejor nos vaya en la exploración de la función a minimizar, ya que es cuando $\lambda$ decrece. En este caso, la convergencia viene heredada de este otro método, por lo que si seguimos decreciendo llegaremos a, al menos, un mínimo local de la función. No tenemos asegurado un mínimo global, pero eso ya lo sabíamos de antemano.

  \item Por otro, si $\lambda$ se hace muy grande, la ecuación que tenemos se ve protagonizada por el propio $\lambda$ y la matriz identidad. Ignorando las matrices que tendrán valores mucho más pequeños, la ecuación que nos queda es

  $$ \lambda I \Delta = -J^T\epsilon $$

  $$ \lambda \Delta = -J^T\epsilon $$

  Notando que $J$ es el jacobiano de la función que estamos minimizando, y que $\epsilon$ es el valor de esta función menos lo que queremos que valga, podemos notar que $J^T\epsilon$ era el vector que aparece en descenso por gradiente con la misma $f$. Esto implica que, si bien en el primer caso heredábamos la convergencia de Gauss-Newton, ahora lo hacemos del descenso por gradiente.
\end{itemize}

Esta dualidad entre Gauss-Newton y descenso por gradiente es el núcleo de la iteración de Levenberg-Marquadt y es lo que la diferencia de estos otros dos métodos. Su flexibilidad le permite aplicar Gauss-Newton cuando la exploración va bien y cada vez se disminuye más el error, y si va mal se acerca al descenso por gradiente para hacer saltos más amplios en la dirección de menor coste.

Una de las funciones que hemos necesitado implementar para realizar este método iterativo ha sido el jacobiano, para poder resolver la ecuación normal aumentada. Como queremos aproximar una derivada, lo que hacemos es simplemente definir un $\delta_i$ que depende de la variable $x_i$ que queramos derivar de $f$ y realizamos la derivada i-esima como $(f(x+\delta_i e_i)-f(x))/\delta_i$ donde e_i denota el i-ésimo vector de la base usual. El resto de la estructura se comenta a continuación:

:::{.algorithm name="\textit{Iteración Levenberg-Marquadt}"}

$\;$

Entrada
:  Como parámetros necesarios, una $f$ que queremos minimizar, un punto 'inicial' desde el que empezamos a explorar y el 'objetivo' al que queremos llegar con $f$ (0, en nuestro caso). Los parámetros opciones son el número máximo de iteraciones en 'max\_iter' y el umbral cerca del objetivo en el que paramos la búsqueda, en 'umbral'.

Salida
:  Un punto $x$ que minimiza la expresión $|f(x)-objetivo|$ y el error, que es justo el valor de $|f(x)-objetivo|$.

1. Inicializa:
   a) $x$ al punto inicial.
   b) iteraciones actuales a 0.
   c) $\lambda$ a 0.0001 (en el código se menciona como 'augment').
   d) $\epsilon$ a $f(x)-objetivo$ y norm a su norma.
2. Mientras 'norm' > 'umbral' y 'iters' < 'max\_iter':
   a) calculamos J como el jacobiano de $f$ en $x$.
   b) conseguimos $\Delta$ resolviendo la ecuación normal aumentada $(J^TJ+\lambda I)\Delta = -J^T\epsilon$.
   c) actualizamos $x$ como $x + \Delta$ si eso disminuye el error, en cuyo caso actualizamos $\epsilon$ y 'norm'.
   d) actualizamos $\lambda$ dividiendo por 10 si hemos disminuido el error, o multiplicando por 10 si no lo hemos hecho.

:::

## Análisis de los resultados

A continuación se exponen varias imágenes resultado de utilizar el algoritmo. Como solo queríamos ver una pequeña muestra de cómo actúa el algoritmo, y de ver que de verdad daba una homografía aceptable, hemos usado las imágenes de yosemite que se nos daban de muestra para las prácticas. Hemos cogido un par de ellas, y cogiendo una lista de correspondencias entre puntos (obtenida de las prácticas anteriores) hemos probado el algoritmo.

Una de las primeras apreciaciones que hemos visto sobre este algoritmo ha sido su efecto sobre outliers en los puntos de la imagen. Hemos observado que el error de Sampson es muy alto cuando hay algunos outliers dentro del dataset de puntos. En el caso que se ve en la figura \ref{malo} podemos ver cómo claramente es posible que hasta empeore la homografía que teníamos hasta ese momento aunque el error de Sampson disminuya. Recordemos que el error de Sampson es una aproximación de la suma geométrica de las distancias. En la figura \ref{malo} se han coloreado unas líneas rojas con estas distancias, empezando en donde está la imagen de un punto frente a su correspondencia. El objetivo es que la suma de longitudes de estas líneas sea lo menor posible. Debido a esto, es posible que el algoritmo esté intentando disminuir esta distancia geométrica acercando los puntos más alejados que son los que afectan más al error, es decir, los outliers. De ahí que cuando intente mejorar el resultado, en realidad empeore nuestra percepción de la homografía final. Podemos fijarnos que hay un corte claro en la ladera de la montaña, donde en la primera imagen se ve de forma continua, aquí encontramos dónde ambas imágenes se mezclan rápidamente. 

\begin{figure}[h!]
\begin{center}
  \includegraphics[width=0.6\textwidth]{Implementacion/resultados/outliers1}

  \includegraphics[width=0.6\textwidth]{Implementacion/resultados/outliers2}

  \caption{Resultado de aplicar el algoritmo con outliers. La primera imagen es el resultado de la homografía inicial (con DLT) y la segunda es la 'mejorada' del error de Sampson, en la que ha empeorado el resultado final.\label{malo}}

  \end{center}

\end{figure}


Esto no es algo que deberíamos dejar pasar en general en la aplicación del algoritmo. Es por esto que antes de pasar la lista de puntos, hacemos un filtro solo con los 100 mejores puntos entre todos los que habíamos encontrado. El filtro lo hacemos mirando la similitud entre descriptores de los puntos. Esto también se menciona en [@hartley2003multiple] , debido justo a que el algoritmo funciona muy bien siempre y cuando los errores y los ruidos en las correspondencias sean pequeños.

Una vez que hacemos esa comprobación, el algoritmo funciona muy bien. A continuación ponemos el ejemplo con la modificación anterior. La podemos ver en la figura \ref{bueno}. Las líneas han pasado a ser puntos, debido a que prácticamente se solapan los puntos respecto a donde deberían estar. Esto es una buena señal.

\begin{figure}[h!]
\begin{center}
  \includegraphics[width=0.6\textwidth]{Implementacion/resultados/bueno1}

  \caption{Resultado de aplicar el algoritmo sin outliers.\label{bueno}}

  \end{center}

\end{figure}

Que el algoritmo devuelva una buena homografía no es lo único que debemos mirar, ya que es necesario también mirar que la solución no es exactamente la misma que la que nos da el algoritmo de inicialización. Es necesario saber si merece la pena todo el análisis y los cálculos adicionales que hemos hecho para la mejora respecto a la inicialización.

Nuestra conclusión, con los ejemplos que hemos visto, es que mejora un poco, pero es prácticamente imperceptible, al menos a ojo. Esto se ha traducido en el error de Samspon también, que cambiaba en un orden de 10^{-6} desde la inicialización al error final, que es una medida algo más concreta y correcta que el cálculo meramente visual. Esto no es culpa del error de Sampson, si no que el algoritmo de inicialización DLT parece ser muy eficiente, al menos en las imágenes que hemos probado (correspondientes al paquete de yosemite).

Esto no quiere decir que todo el trabajo y que el algoritmo sea inútil, ni mucho menos. Este algoritmo tendrá utilidad en aplicaciones concretas en las que la precisión y exactitud de la homografía en cuestión sean de importancia vital para el correcto funcionamiento del programa. Ejemplos en los que esto podría ser cierto son imágenes de medidas de piezas de un taller o una empresa, a partir de una serie de fotos tomadas de ellas y una distancia referencia dentro de la imagen (un metro situado al lado de las piezas, por ejemplo). En este caso, la exactitud de la homografía es lo que permite calcular con precisión el tamaño de las piezas, y podría ser una posible aplicación de esta mejora del cálculo de homografías.

Por otro lado, siempre que tengamos los recursos y tiempo necesarios, es también importante notar que como el resultado del error de Sampson siempre será una mejora del de inicialización, nunca está de más usar este en lugar del otro.


## Propuestas de mejora

El cálculo del error de Sampson es considerablemente lento. Aunque es cierto que el programa tarda en ejecutarse en unos diez segundos, si queremos aplicar este algoritmo en cualquier aplicación de tiempo real (que es algo que es probable, teniendo en cuenta los usos de este tipo de programas) es necesario que la aplicación se quede en dos órdenes de magnitud de segundos de los que tarda ahora. Tras analizar qué era lo que más tiempo necesitaba, hemos visto que era el cálculo del error de Sampson y, en concreto, el cálculo de la matriz C_H(X) mencionada en la sección [Error de Sampson]. Una posible mejora sería hacer más eficiente estos cálculos mediante el uso de métodos numéricos o librerías optimizadas para este tipo de cálculos.

Dentro de las mejoras de partes concretas del algoritmo, no hemos llegado a implementar algunas de las que se comentan en [@hartley2003multiple]. Algunas de ellas son la mejora del algoritmo de Levenberg-Marquadt para mejorar su eficiencia (llamada 'sparse Levenberg-Marquadt'), o la mejora sobre el cálculo del error de Sampson para tener en cuenta una distribución de error en las medidas de las correspondencias que no sea usual y que venga dada por una matriz de covarianza. Ambas mejoras vienen bien explicadas y aportarían una mejor adaptación del algoritmo a otro tipo de situaciones.

Por otro lado, una forma educativa, bonita y visual de ver cómo funciona la optimización del error de Sampson sería ir construyendo una secuencia de imágenes con las homografías de cada iteración. Creemos que esto es una mejora interesante viendo el proyecto como una herramienta de aprendizaje para un alumno que haya cursado Visión por Computador y quiera profundizar un poco más en las estimaciones de homografías. Esta mejora podría ser de utilidad para ese caso concreto.

# Apéndice: Funcionamiento del código adjunto {.unnumbered}
