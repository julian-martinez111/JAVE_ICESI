---

## title: "Fundamentos Teóricos y de Valoración de los Variance Swaps: Del Lema de Itô al CBOE VIX" author: "Julián Andrés Martínez Ortiz" date: "`r Sys.Date()`" output: github\_document: toc: true toc\_depth: 3 pdf\_document: toc: true

knitr::opts\_chunk\$set(echo \= TRUE)

# Introducción y Evolución Histórica

La evolución teórica y práctica de los *variance swaps* (swaps de varianza) durante la década de 1990 estuvo impulsada por la necesidad institucional de negociar volatilidad pura sin el riesgo direccional inherente a las opciones con cobertura delta (*delta-hedged options*).

* **El hito de Neuberger (1994):** Demostró analíticamente que la varianza realizada podía replicarse de forma perfecta mediante el rebalanceo continuo de una cobertura delta sobre un "contrato logarítmico" (*log contract*), un instrumento cuyo pago es igual al logaritmo del precio del activo \[cite: 1\].  
* **El teorema de replicación de Carr y Madan (1998):** Dado que los contratos logarítmicos puros no se negocian directamente en los mercados organizados, Carr y Madan demostraron matemáticamente que cualquier pago dos veces continuamente diferenciable (incluyendo el pago logarítmico) puede replicarse estáticamente utilizando un portafolio estandarizado de opciones europeas de compra (*calls*) y venta (*puts*) \[cite: 1\].  
* **Formalización de Goldman Sachs (1999):** El grupo de Estrategias Cuantitativas de Goldman Sachs (Demeterfi, Derman, Kamal y Zou) unificó estos conceptos en una metodología estándar para la industria, sentando las bases matemáticas que hoy en día soportan la modernización del índice VIX de la CBOE \[cite: 1\].

---

# Estructura, Mecánica y Liquidación

Un *variance swap* es un contrato forward sobre la varianza anualizada de un activo subyacente específico. Su pago al vencimiento \$T\$ es lineal respecto a la varianza \[cite: 1\]:

\$\$\\text{Payoff} \= (\\sigma\_{R}^{2} \- K\_{\\text{var}}) \\times N\_{\\text{vol}}\$\$

Donde \$\\sigma\_{R}^{2}\$ representa la varianza realizada del activo durante la vida del contrato, \$K\_{\\text{var}}\$ es el precio de entrega (strike de varianza implícita) y \$N\_{\\text{vol}}\$ es el monto nocional por punto de varianza \[cite: 1\].

> **Implicación de liquidación (Settlement al vencimiento):** Los *variance swaps* se liquidan únicamente al vencimiento (estilo europeo) \[cite: 1\]. Si el contrato se liquidara diariamente (*mark-to-market* diario), la parte compradora enfrentaría un severo problema de convexidad. Como la varianza realizada se calcula como un promedio de rendimientos al cuadrado en todo el periodo, un día de alta volatilidad al inicio forzaría un pago masivo que podría no estar justificado si el resto del periodo presenta una volatilidad cercana a cero. La liquidación al vencimiento garantiza que se capture el promedio real independiente de la trayectoria (*path-independent*).

---

# Varianza vs. Volatilidad Pura

Aunque comercialmente se suelen ofrecer como "apuestas de volatilidad pura", esto es impreciso: son instrumentos de **varianza pura** (\$\\sigma^2\$). La razón fundamental por la cual los mercados operan varianza en lugar de volatilidad es una propiedad matemática clave: **la aditividad** \[cite: 1\].

1. **La varianza es aditiva:** La varianza acumulada en un mes es la suma de las varianzas diarias. Por ende, la acumulación diaria de varianza realizada se puede compensar perfectamente con la acumulación diaria y lineal de pérdidas y ganancias (P\&L) de una posición accionaria con cobertura delta \[cite: 1\].  
2. **La volatilidad no es aditiva:** La raíz cuadrada de una suma no es igual a la suma de las raíces cuadradas. Replicar un swap de volatilidad verdadero requiere un pago no lineal sobre la varianza acumulada, exponiendo al creador de mercado (*dealer*) a la temida "volatilidad de la volatilidad" (*vol-of-vol*) \[cite: 1\].

---

# Derivación Matemática Explícita del Strike de Varianza

El objetivo de esta sección es encontrar el precio de entrega justo (\$K\_{\\text{var}}\$), el cual, por la teoría de precios de no-arbitraje, debe ser igual al valor esperado neutral al riesgo de la varianza realizada \[cite: 1\].

## El Proceso de Precios Continuo

El precio de una acción sigue un movimiento browniano geométrico. Esto significa que el precio se mueve por dos fuerzas: una parte “ordenada” o promedio (\$\\mu\$), y otra parte aleatoria (\$\\sigma\_t dW\_t\$) que representa el ruido del mercado. Matemáticamente lo escriben como el retorno porcentual del activo, es decir, \$\\frac{dS\_t}{S\_t}\$, porque en finanzas importa más cuánto cambia en porcentaje que en dólares absolutos.

Donde \$\\mu\_t\$ es la deriva (*drift*), \$\\sigma\_t\$ es la volatilidad estocástica instantánea y \$dW\_t\$ es un proceso de Wiener estándar \[cite: 1\]. La varianza realizada anualizada \$V\$ se define como \[cite: 1\]:

\$\$V \= \\frac{1}{T}\\int\_{0}^{T}\\sigma\_{t}^{2}dt\$\$

## Aplicación del Lema de Itô y Transformación Logarítmica

Luego el objetivo de la investigación cuantitativa es encontrar una manera de expresar la volatilidad realizada (*realized variance*) usando únicamente movimientos observables del precio en el mercado. El problema es que la volatilidad \$\\sigma\_t\$ no se observa directamente en el mercado; lo único observable es el precio de la acción. Por eso buscan una transformación matemática que haga aparecer el término \$\\sigma\_t^2\$, que es la varianza instantánea.

Para lograrlo usan el logaritmo del precio, \$\\ln(S\_t)\$. Esto no es casualidad: el logaritmo y sus dos derivadas hacen que, al aplicarle el Lema de Itô, aparezca automáticamente un término relacionado con la volatilidad. Cuando aplican Itô al logaritmo del precio, sustituyendo las derivadas del logaritmo y expandiendo el \$(dS\_t)^2\$ para obtener:

\$\$\\mu\_t^2 S\_t^2 (dt)^2 \+ 2\\mu\_t \\sigma\_t S\_t^2 dt , dW\_t \+ \\sigma\_t^2 S\_t^2 (dW\_t)^2\$\$

y cancelar los dos primeros términos utilizando las reglas fundamentales del cálculo estocástico, obtenemos una función que simboliza el retorno porcentual menos un ajuste de volatilidad. Ese ajuste es precisamente \$-\\frac{1}{2}\\sigma\_t^2 dt\$, y aparece gracias al término cuadrático del movimiento browniano, es decir, porque en cálculo estocástico \$(dW\_t)^2 \= dt\$.

> **Nota sobre el cálculo estocástico:** Recordemos que \$(dW\_t)^2 \= dt\$ porque el movimiento browniano se comporta de una manera muy distinta a una función normal y suave. En cálculo tradicional, cuando haces un cambio muy pequeño \$dx\$, su cuadrado \$(dx)^2\$ es muchísimo más pequeño todavía, así que se cancela de la expansión de Taylor. Pero el Browniano no cambia “suavemente”; cambia de forma extremadamente zigzagueante y rugosa. Sus pequeños movimientos aleatorios son mucho más grandes de lo que intuitivamente esperarías para intervalos de tiempo diminutos. La clave está en cómo escala el movimiento browniano. Un incremento browniano en un intervalo muy pequeño \$dt\$ no tiene tamaño proporcional a \$dt\$, sino proporcional a la raíz cuadrada del tiempo: \$dW\_t \\sim \\sqrt{dt}\$. Al elevar esta expresión al cuadrado, se cancela la raíz dejándonos con la expresión final (\$dt\$ \= un pequeño paso del tiempo \= diferencial del tiempo).

Después reorganizan la ecuación para dejar sola la volatilidad. Ahí es donde ocurre la parte importante: consiguen escribir \$\\sigma\_t^2 dt\$ en función de cantidades relacionadas únicamente con el precio y su logaritmo.

Finalmente integran toda la expresión desde el tiempo inicial hasta el vencimiento \$T\$. Eso acumula toda la volatilidad a lo largo del periodo y produce la fórmula final de la *realized variance*. El resultado muestra que la varianza realizada puede replicarse usando dos componentes: una estrategia dinámica sobre la acción y una posición sobre un *payoff* logarítmico.

## Esperanza Neutral al Riesgo y el Teorema de Carr-Madan

Posteriormente, para encontrar un término negociable (*tradeable*) para el *payoff* logarítmico, se mueven a la llamada medida neutral al riesgo (\$\\mathbb{Q}\$). Esta es una herramienta muy usada en derivados porque simplifica el *pricing*. Bajo esta medida, se asume que el activo crece en promedio a la tasa libre de riesgo \$r\$. La idea no es que el mundo real funcione exactamente así, sino que bajo esa medida matemática los derivados pueden valorarse como expectativas descontadas. Entonces el *strike* justo del *variance swap* se define como la esperanza neutral al riesgo de la *realized variance* futura.

Entonces el *paper* introduce el teorema de replicación de Carr-Madan. La idea central del teorema es que prácticamente cualquier *payoff* suficientemente suave puede construirse combinando muchas opciones *vanilla* europeas de distintos *strikes*. Matemáticamente utilizan una expansión tipo Taylor con integrales, donde cualquier función puede descomponerse en una parte lineal (que se replica con caja, bonos, futuros o *forwards* del subyacente) más una combinación continua e infinita de *calls* y *puts*.

* La fórmula general de Carr-Madan descompone cualquier función \$f(S\_T)\$ en cuatro partes. La primera parte es un término constante, \$f(S\_*)\$, que simplemente representa el valor de la función en un punto de referencia arbitrario llamado \$S\_*\$, que en la práctica suele ser el precio del *forward* del activo. La segunda parte es un término lineal, \$f'(S\_*)(S\_T \- S\_*)\$, que puede replicarse fácilmente utilizando *forwards* o futuros porque depende linealmente del precio del activo.  
* Las otras dos partes son las más importantes. Una integral utiliza *puts* europeas y la otra utiliza *calls* europeas. Esto ocurre porque los términos \$(K \- S\_T)^+\$ y \$(S\_T \- K)^+\$ son exactamente los *payoffs* de *puts* y *calls* respectivamente. De esta manera, cualquier curvatura o convexidad de la función puede reconstruirse usando opciones distribuidas sobre todos los *strikes* posibles.

Financieramente, esto significa que la volatilidad futura implícita puede extraerse observando precios de muchas opciones diferentes. No basta con mirar una sola opción *ATM*; es necesario integrar información proveniente de toda la superficie de volatilidad. Por eso los *variance swaps* y productos como el VIX utilizan una gran cantidad de *strikes* simultáneamente.

Al aplicar esperanza bajo la medida neutral al riesgo (utilizan probabilidades neutrales al riesgo para los precios del mercado), el precio esperado futuro del activo es:

\$\$\\mathbb{E}^{\\mathbb{Q}}\[S\_T\] \= S\_0 e^{rT}\$\$

Ahora ya no aparecen los *payoffs* \$(K \- S\_T)^+\$ o \$(S\_T \- K)^+\$, sino directamente los precios de las opciones hoy:

* \$P(K)\$: precio de *puts*  
* \$C(K)\$: precio de *calls*

Esto ocurre porque al tomar esperanza neutral al riesgo, el valor esperado descontado de un *payoff* es precisamente el precio actual de la opción \[cite: 1\].

Evaluando la expansión en torno a un precio de referencia \$S\_\*\$, se llega a la fórmula analítica final de valoración \[cite: 1\]:

\$\$K\_{\\text{var}} \= \\frac{2}{T}\\left\[ rT \- \\left(\\frac{S\_0}{S\_*}e^{rT} \- 1\\right) \- \\ln\\left(\\frac{S\_*}{S\_0}\\right) \+ e^{rT}\\int\_{0}^{S\_{*}}\\frac{P(K)}{K^{2}}dK \+ e^{rT}\\int\_{S\_*}^{\\infty}\\frac{C(K)}{K^{2}}dK \\right\]\$\$

Esta formulación demuestra que **la varianza se puede valorar y replicar estáticamente comprando un *strip* infinito de opciones OTM (fuera del dinero) de *puts* y *calls*, ponderadas exactamente por el inverso del cuadrado de su *strike* (\$\\frac{1}{K^2}\$)** \[cite: 1\].

---

# Discretización Práctica y el CBOE VIX

En la práctica de mercado, una integral continua de opciones es inviable debido a que las opciones solo se negocian en intervalos de *strikes* discretos. Los creadores de mercado deben aproximar la integral continua mediante una sumatoria discreta \$\\sum \\Delta K\$ \[cite: 1\]:

\$\$K\_{\\text{var}} \\approx \\frac{2}{T}\\left\[ rT \- \\left(\\frac{S\_0}{S\_*}e^{rT} \- 1\\right) \- \\ln\\left(\\frac{S\_*}{S\_0}\\right) \+ e^{rT}\\sum\_{i}\\frac{\\Delta K\_i}{K\_i^2}Q(K\_i) \\right\]\$\$

* **El Motor del Índice VIX:** Esta sumatoria discreta es el motor matemático detrás del índice VIX de la CBOE. La metodología evalúa un amplio abanico de opciones sobre el S\&P 500, ponderando el precio medio de cada opción por \$\\frac{1}{K^2}\$ sobre un horizonte constante de 30 días, aplicando la raíz cuadrada para expresarlo como un porcentaje de volatilidad implícita anualizada \[cite: 1\].  
* **El Riesgo de Salto (*Jump Risk*):** Teóricamente elegante, la discretización expone a las mesas de dinero a riesgos operativos y de modelo. Si el subyacente experimenta un *gap* de precios violento sobre una región donde no existen *strikes* de opciones líquidos, el *strip* discreto de \$\\frac{1}{K^2}\$ falla en replicar de forma perfecta el pago logarítmico, generando un error de seguimiento (*tracking error*) importante \[cite: 1\].

---

# Estrategia y Prima de Riesgo de Varianza (VRP)

El VRP (*Variance Risk Premium* o Prima de Riesgo de Varianza) se define como la diferencia entre la varianza implícita (lo que el mercado espera que ocurra, \$K\_{\\text{var}}\$) y la varianza realizada:

\$\$\\text{VRP} \= K\_{\\text{var}} \- \\sigma\_R^2\$\$

* **\$\\text{VRP} \> 0\$:** La volatilidad está cara; esto favorece una posición corta (*Short*).  
* **\$\\text{VRP} \< 0\$:** La volatilidad está barata; esto favorece una posición larga (*Long*).

Cuando hay un fuerte *contango* (los futuros mucho más altos que el *spot*), suele ser porque el mercado está relativamente calmado. En estos entornos, los futuros del VIX tienden a bajar con el paso del tiempo (convergiendo hacia el VIX *Spot*). Esto crea un efecto negativo para quien está *short*, y potencialmente positivo para quien entra *long* si la volatilidad sube de repente.

---

# Referencias

1. Carr, Peter, and Dilip Madan. "Towards a Theory of Volatility Trading." *Volatility: New Estimation Techniques for Pricing Derivatives*, edited by R. Jarrow, 1998, pp. 417-427 \[cite: 1\].  
2. Demeterfi, Kresimir, Emanuel Derman, Michael Kamal, and Joseph Zou. "More Than You Ever Wanted To Know About Volatility Swaps." *Goldman Sachs Quantitative Strategies Research Notes*, March 1999 \[cite: 1\].