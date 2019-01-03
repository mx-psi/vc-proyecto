# Introducción

TODO: Describir el problema en términos generales, sus aplicaciones, qué hemos hecho en el trabajo, etc.

# El problema de estimación de homografías

:::{.problem}
Dado un conjunto de correspondencias $\{\mathbf{x}_i \leftrightarrow \mathbf{x}'_i\}_{i = 1, \dots, N}$ en coordenadas homogéneas,
hallar una homografía dada por su expresión matricial $H \in \mathcal{M}_{3 \times 3}(\mathbb{R})$ que minimice el error entre $H\mathbf{x}_i$ y $\mathbf{x}'_i$.
:::

- TODO: Describir el cálculo en el caso exacto para 4 correspondencias. La expresión matricial por bloques del sistema de ecuaciones que define una correspondencias, notando $\mathbf{x}_i = (x_i, y_i, w_i)$, $\mathbf{x}'_i = (x'_i, y'_i, w'_i)$ y $H = (h^1 h^2 h^3)^T$ por filas es (4.3 p.89)
$$A_iH =
\begin{matrix}\mathbf{0}^T & -w'_i\mathbf{x}_i^T & y'_i\mathbf{x}_i^T \\
w'_i\mathbf{x}_i^T  & \mathbf{0}^T  & -x'_i\mathbf{x}_i^T \end{matrix}
\begin{matrix}h^1 \\ h^2 \\ h^3\end{matrix}$$

# Error de Sampson

Si notamos cada fila de la matriz $A_i$ como $A^{(j)}$ y notamos $\mathbf{X} = (x,y,x',y')$ entonces tenemos que:

$$J = \frac{\partial \mathcal{C}_H}{\mathbf{X}} = (\frac{\partial A^{(i)}}{\partial \mathbf{X}_j})_{j = 1 \dots 4, i = 1,2}$$
