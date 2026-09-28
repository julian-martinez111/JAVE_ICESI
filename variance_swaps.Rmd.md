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

Asumimos que el precio accionario subyacente \$S\_t\$ sigue un proceso de difusión estándar bajo la medida real \$\\mathbb{P}\$ \[cite: 1\]:

\$\$\\frac{dS\_t}{S\_t} \= \\mu\_t dt \+ \\sigma\_t dW\_t\$\$

Donde \$\\mu\_t\$ es la deriva (*drift*), \$\\sigma\_t\$ es la volatilidad estocástica instantánea y \$dW\_t\$ es un proceso de Wiener estándar \[cite: 1\]. La varianza realizada anualizada \$V\$ se define como \[cite: 1\]:

\$\$V \= \\frac{1}{T}\\int\_{0}^{T}\\sigma\_{t}^{2}dt\$\$

## Aplicación del Lema de Itô y Transformación Logarítmica

Para aislar analíticamente el término de varianza \$\\sigma\_t^2 dt\$, aplicamos el Lema de Itô a la función logarítmica \$f(S\_t) \= \\ln(S\_t)\$, cuya segunda derivada es \$-\\frac{1}{S\_t^2}\$ (lo que permite cancelar el término \$S\_t^2\$ generado por \$(dS\_t)^2\$) \[cite: 1\]:

\$\$d(\\ln S\_t) \= \\frac{dS\_t}{S\_t} \- \\frac{1}{2}\\sigma\_{t}^{2}dt\$\$

Reordenando los términos para despejar el componente de varianza \[cite: 1\]:

\$\$\\frac{1}{2}\\sigma\_{t}^{2}dt \= \\frac{dS\_t}{S\_t} \- d(\\ln S\_t)\$\$

Integrando ambos lados desde \$t=0\$ hasta \$t=T\$ y multiplicando por \$\\frac{2}{T}\$ obtenemos la expresión de la varianza realizada \[cite: 1\]:

\$\$V \= \\frac{2}{T}\\left\[ \\int\_{0}^{T}\\frac{dS\_t}{S\_t} \- \\ln\\left(\\frac{S\_T}{S\_0}\\right) \\right\]\$\$

> **Intuición Financiera:** Esta ecuación demuestra que la varianza realizada puede replicarse manteniendo una posición en acciones continuamente rebalanceada (manteniendo siempre \$\\frac{1}{S\_t}\$ acciones) y tomando una posición corta estática en un contrato logarítmico que paga \$\\ln(S\_T/S\_0)\$ \[cite: 1\].

## Esperanza Neutral al Riesgo y el Teorema de Carr-Madan

Para trasladar la valuación al mercado mediante fijación de precios libres de arbitraje, pasamos a la medida neutral al riesgo \$\\mathbb{Q}\$, donde el activo rinde la tasa libre de riesgo \$r\$ \[cite: 1\]. El strike justo \$K\_{\\text{var}}\$ se define como \[cite: 1\]:

\$\$K\_{\\text{var}} \= \\mathbb{E}^{\\mathbb{Q}}\[V\] \= \\frac{2}{T}\\left( rT \- \\mathbb{E}^{\\mathbb{Q}}\\left\[\\ln\\left(\\frac{S\_T}{S\_0}\\right)\\right\] \\right)\$\$

Dado que el pago logarítmico no se negocia directamente, utilizamos el teorema de expansión de Carr-Madan con opciones europeas de extinción vainilla \[cite: 1\]. Evaluando la expansión en torno a un precio de referencia \$S\_\*\$, se llega a la fórmula analítica final de valoración \[cite: 1\]:

\$\$K\_{\\text{var}} \= \\frac{2}{T}\\left\[ rT \- \\left(\\frac{S\_0}{S\_*}e^{rT} \- 1\\right) \- \\ln\\left(\\frac{S\_*}{S\_0}\\right) \+ e^{rT}\\int\_{0}^{S\_{*}}\\frac{P(K)}{K^{2}}dK \+ e^{rT}\\int\_{S\_*}^{\\infty}\\frac{C(K)}{K^{2}}dK \\right\]\$\$

Esta formulación demuestra que **la varianza se puede valorar y replicar estáticamente comprando un strip infinito de opciones OTM (fuera del dinero) de *puts* y *calls*, ponderadas exactamente por el inverso del cuadrado de su strike (\$\\frac{1}{K^2}\$)** \[cite: 1\].

---

# Discretización Práctica y el CBOE VIX

En la práctica de mercado, una integral continua de opciones es inviable debido a que las opciones solo se negocian en intervalos de *strikes* discretos. Los creadores de mercado deben aproximar la integral continua mediante una sumatoria discreta \$\\sum \\Delta K\$ \[cite: 1\]:

\$\$K\_{\\text{var}} \\approx \\frac{2}{T}\\left\[ rT \- \\left(\\frac{S\_0}{S\_*}e^{rT} \- 1\\right) \- \\ln\\left(\\frac{S\_*}{S\_0}\\right) \+ e^{rT}\\sum\_{i}\\frac{\\Delta K\_i}{K\_i^2}Q(K\_i) \\right\]\$\$

* **El Motor del Índice VIX:** Esta sumatoria discreta es el motor matemático detrás del índice VIX de la CBOE. La metodología evalúa un amplio abanico de opciones sobre el S\&P 500, ponderando el precio medio de cada opción por \$\\frac{1}{K^2}\$ sobre un horizonte constante de 30 días, aplicando la raíz cuadrada para expresarlo como un porcentaje de volatilidad implícita anualizada \[cite: 1\].  
* **El Riesgo de Salto (*Jump Risk*):** Teóricamente elegante, la discretización expone a las mesas de dinero a riesgos operativos y de modelo. Si el subyacente experimenta un *gap* de precios violento sobre una región donde no existen *strikes* de opciones líquidos, el *strip* discreto de \$\\frac{1}{K^2}\$ falla en replicar de forma perfecta el pago logarítmico, generando un error de seguimiento (*tracking error*) importante \[cite: 1\].

---

# Referencias

1. Carr, Peter, and Dilip Madan. "Towards a Theory of Volatility Trading." *Volatility: New Estimation Techniques for Pricing Derivatives*, edited by R. Jarrow, 1998, pp. 417-427 \[cite: 1\].  
2. Demeterfi, Kresimir, Emanuel Derman, Michael Kamal, and Joseph Zou. "More Than You Ever Wanted To Know About Volatility Swaps." *Goldman Sachs Quantitative Strategies Research Notes*, March 1999 \[cite: 1\].