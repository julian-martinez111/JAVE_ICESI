# Variance Swaps alrededor de Earnings: Marco Teórico con Merton (1976) y Bates (1996)

> Valoración y replicación del *fair strike* de un variance swap bajo saltos y volatilidad estocástica, con aplicación a ventanas de earnings.

## Tabla de contenidos

- [Marco Teórico](#marco-teórico)
  - [1. Preliminares de cálculo estocástico](#1-preliminares-de-cálculo-estocástico)
  - [2. El variance swap: definición, replicación y sesgos](#2-el-variance-swap-definición-replicación-y-sesgos)
  - [3. Modelo de Merton (jump-diffusion)](#3-modelo-de-merton-jump-diffusion)
  - [4. Modelo de Heston y modelo de Bates (SV + saltos)](#4-modelo-de-heston-y-modelo-de-bates-sv--saltos)
  - [5. Modelar el evento de earnings y el vol crush](#5-modelar-el-evento-de-earnings-y-el-vol-crush)
  - [6. Valoración de opciones por transformadas de Fourier](#6-valoración-de-opciones-por-transformadas-de-fourier)
  - [7. Calibración](#7-calibración)
  - [8. Simulación Monte Carlo](#8-simulación-monte-carlo)
- [Bibliografía](#bibliografía)

---

## Marco Teórico

**Notación.** $S_t$ precio del subyacente, $r$ tasa libre de riesgo, $q$ dividend yield continuo, $T$ vencimiento (en años de *tiempo de negociación*, es decir, días hábiles/252), $F_0=S_0e^{(r-q)T}$ forward, $X_t=\ln S_t$, $\mathbb{Q}$ medida neutral al riesgo, $\mathbb{P}$ medida física, $t_e$ fecha del anuncio de resultados.

---


### 2. El variance swap: definición, replicación y sesgos
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
#### 2.1 Contrato

El Riesgo de Salto (Jump Risk): Teóricamente elegante, la discretización expone a las mesas de dinero a riesgos operativos y de modelo. Si el subyacente experimenta un gap de precios violento sobre una región donde no existen strikes de opciones líquidos, el strip discreto de  
K 
2
 
1
​	
  falla en replicar de forma perfecta el pago logarítmico, generando un error de seguimiento (tracking error) importante.
Un variance swap paga al vencimiento

$$
\text{Payoff}=N_{var}\,\big(RV-K_{var}\big),\qquad
N_{var}=\frac{N_{vega}}{2\,\sqrt{K_{var}}},
$$

con $RV$ y $K_{var}$ en unidades de varianza anualizada (en la convención de mercado, puntos de volatilidad al cuadrado). La varianza realizada *discreta* estándar es

$$
RV=\frac{A}{n}\sum_{i=1}^{n}\Big(\ln\frac{S_{t_i}}{S_{t_{i-1}}}\Big)^2,\qquad A=252,
$$

**sin restar la media** (convención de mercado), sobre cierres diarios. Los ajustes por dividendos y por "días de mercado" varían según la confirmación del contrato. El swap es una posición larga en varianza si el comprador recibe $RV$ y paga $K_{var}$.

Propiedades que motivan su uso (Demeterfi et al., 1999; Bossu, Strasser y Guichard, 2005):

* La varianza es **aditiva** en el tiempo; la volatilidad no.
* El **dollar gamma** de la tira de replicación es constante: $\Gamma=\frac{2}{TS^2}$, de modo que $S^2\Gamma=2/T$ para cualquier $S$.
* El vega en varianza decae linealmente: $\partial V/\partial\sigma^2=\tau/T$ con $\tau=T-t$.
* Valor de mercado en $t$ (contrato ya iniciado):

$$
V_t=e^{-r(T-t)}\,N_{var}\Big[\tfrac{t}{T}\,RV_{0,t}+\tfrac{T-t}{T}\,K_{t,T}-K_{var}\Big],
$$

donde $K_{t,T}$ es el strike *forward* vigente.

#### 2.2 Del contrato a la esperanza

Bajo $\mathbb{Q}$, con tasas deterministas, $K_{var}=\mathbb{E}^{\mathbb{Q}}[RV]$ (ignorando el descuento del pago único). La cuestión es **qué versión de $RV$**: (i) la variación cuadrática continua $[X]_T/T$, (ii) la suma discreta de cuadrados de log-retornos, (iii) la cantidad replicable por opciones. Bajo difusión pura las tres coinciden; con saltos no.

#### 2.3 Replicación estática (Neuberger, Dupire, Carr–Madan)

**Paso 1 (identidad de Itô).** Con trayectorias continuas, $d\ln S=\frac{dS}{S}-\frac12\sigma^2dt$, luego

$$
\int_0^T\sigma_t^2\,dt=2\Big(\int_0^T\frac{dS_t}{S_t}-\ln\frac{S_T}{S_0}\Big).
$$

La integral estocástica corresponde a una posición **dinámica** de $2/S_t$ acciones (delta-hedge en futuros); el segundo término es un **contrato logarítmico** (*log contract*).

**Paso 2 (spanning de Carr–Madan).** Para cualquier $S^*>0$,

$$
\ln\frac{S_T}{S^*}=\frac{S_T-S^*}{S^*}-\int_0^{S^*}\frac{(K-S_T)^+}{K^2}\,dK-\int_{S^*}^{\infty}\frac{(S_T-K)^+}{K^2}\,dK .
$$

**Paso 3 (esperanza bajo $\mathbb{Q}$).** Con $\mathbb{E}[\int dS/S]=(r-q)T$ y opciones OTM $P(K),C(K)$ con precios *de hoy*:

$$
\boxed{K_{var}(S^*)=\frac{2}{T}\Big[(r-q)T-\Big(\frac{F_0}{S^*}-1\Big)-\ln\frac{S^*}{S_0}+e^{rT}\!\int_0^{S^*}\!\frac{P(K)}{K^2}dK+e^{rT}\!\int_{S^*}^{\infty}\!\frac{C(K)}{K^2}dK\Big]}
$$

Esta es la expresión general de Demeterfi et al. (1999) escrita con dividendos. **Al elegir $S^*=F_0$ todos los términos no integrales se cancelan** y queda la forma limpia (Britten-Jones y Neuberger, 2000):

$$
\boxed{K_{var}=\frac{2}{T}\,e^{rT}\Big[\int_0^{F_0}\frac{P(K)}{K^2}\,dK+\int_{F_0}^{\infty}\frac{C(K)}{K^2}\,dK\Big]}
$$

> **Nota sobre la fórmula del planteamiento original.** Si se escribe la versión con $-\ln(S_0/F_0)$ y sin el factor $e^{rT}$ ni el término lineal, esa expresión solo es correcta bajo supuestos particulares ($r=q=0$ o $S^*=F_0$ con los términos ya cancelados). Usar la forma anterior evita errores de $O(rT)$, que en vencimientos de 5–10 días son pequeños pero no nulos, y errores de $O(\text{dividendos})$ si hay ex-date en la ventana.

#### 2.4 Discretización y truncamiento con strikes de mercado

Con strikes discretos $K_1<\dots<K_m$, la versión operativa (tipo VIX) es

$$
K_{var}\approx\frac{2}{T}\sum_{i=1}^m\frac{\Delta K_i}{K_i^2}\,e^{rT}Q(K_i)-\frac1T\Big(\frac{F_0}{K_0}-1\Big)^2,
$$

donde $Q(K_i)$ es la opción OTM (put si $K_i<K_0$, call si $K_i>K_0$; promedio en $K_0$), $K_0$ es el mayor strike $\le F_0$ y $\Delta K_i=\frac{K_{i+1}-K_{i-1}}{2}$.

Fuentes de error (Jiang y Tian, 2005, 2007):

1. **Truncamiento:** la cadena real no llega a $K\to0,\infty$. Se subestima $K_{var}$; más grave con salto grande de earnings (colas gordas).
2. **Discretización:** malla $\Delta K$ gruesa en near-expiry.
3. **Ruido de precios:** *bid-ask* ancho en wings de vencimientos semanales.

**Práctica recomendada:** (a) interpolar/extrapolar la volatilidad implícita (SVI o splines en delta/log-moneyness, con extrapolación plana en los wings; Gatheral–Jacquier, 2014; Fengler, 2009), (b) reconstruir precios OTM sobre una malla fina con Black–Scholes y (c) integrar numéricamente. Reportar $K_{var}$ con y sin extrapolación como análisis de robustez.

**Varianza realizada discreta bajo un modelo.** Con retornos de paso $\Delta_i$,

$$
\mathbb{E}[RV]=\frac{A}{n}\sum_{i=1}^n\mathbb{E}\big[R_i^2\big],\qquad \mathbb{E}[R_i^2]=-\varphi_{R_i}''(0)=\operatorname{Var}(R_i)+\mathbb{E}[R_i]^2,
$$

que puede calcularse por derivación (analítica o numérica) de la función característica del retorno de paso. Broadie y Jain (2008) muestran que el efecto de la discretización es típicamente pequeño, mientras que el efecto de los saltos puede ser significativo.

#### 2.5 Resultado clave: el sesgo de replicación por saltos

Con saltos, la tira de opciones replica

$$
K^{\text{repl}}=\frac{2}{T}\,\mathbb{E}^{\mathbb{Q}}\Big[\int_0^T\frac{dS_t}{S_{t^-}}-\ln\frac{S_T}{S_0}\Big],
$$

mientras que la variación cuadrática del log-precio es

$$
K^{QV}=\frac1T\,\mathbb{E}^{\mathbb{Q}}\big[[X]_T\big]=\frac1T\,\mathbb{E}^{\mathbb{Q}}\Big[\int_0^T\sigma_t^2dt+\sum_{s\le T}(\Delta X_s)^2\Big].
$$

**Diferencia pathwise.** Como $\ln S_T/S_0=\int dS/S-\tfrac12\int\sigma^2dt-\sum\big(e^{\Delta X}-1-\Delta X\big)$,

$$
[X]_T-2\Big(\int\tfrac{dS}{S}-\ln\tfrac{S_T}{S_0}\Big)=-2\sum_{s\le T}\Big(e^{\Delta X_s}-1-\Delta X_s-\tfrac12\Delta X_s^2\Big)\approx-\tfrac13\sum_{s\le T}(\Delta X_s)^3 .
$$

**Consecuencia:**

$$
\boxed{K^{\text{repl}}-K^{QV}=2\lambda\,\mathbb{E}\Big[e^{J}-1-J-\tfrac12J^2\Big]\approx\tfrac{\lambda}{3}\,\mathbb{E}[J^3]}
$$

* Si los saltos son **negativos en media** ($\mu_J<0$, típico: ajuste bajista en EA), $\mathbb{E}[J^3]<0$ y $K^{\text{repl}}<K^{QV}$: el strike model-free **subestima** la varianza cuadrática esperada.
* El sesgo es de tercer orden en el tamaño del salto: se dispara con saltos grandes como los de earnings.
* En un modelo de saltos **el swap (que paga $RV$ discreta) vale $\approx K^{QV}$, no $K^{\text{repl}}$**. Por eso el análisis debe reportar ambos y no mezclar definiciones al comparar modelo vs. mercado.

Este punto es tratado sistemáticamente por Broadie y Jain (2008) y Carr y Lee (2009); Jarrow, Kchia, Larsson y Protter (2013) estudian la aproximación de la versión discreta por su límite continuo.

#### 2.6 Variance swap vs. volatility swap (convexidad)

Como $\sqrt{\cdot}$ es cóncava, por Jensen $\mathbb{E}[\sqrt{RV}]<\sqrt{\mathbb{E}[RV]}$ y, a segundo orden,

$$
K_{vol}\approx\sqrt{K_{var}}-\frac{\operatorname{Var}(RV)}{8\,K_{var}^{3/2}} .
$$

Con saltos, $\operatorname{Var}(RV)$ es muy grande y esta aproximación puede fallar (Broadie y Jain, 2008); el análisis del proyecto se centra en varianza, que es replicable, y en $\sqrt{K_{var}}$ solo como unidad de reporte en puntos de volatilidad.

#### 2.7 Varianza *forward* (calendar) e implicaciones para earnings

Por aditividad, entre dos vencimientos $T_1<T_2$:

$$
K_{T_1,T_2}=\frac{T_2K_{T_2}-T_1K_{T_1}}{T_2-T_1}.
$$

Es la herramienta model-free para leer la **estructura temporal de varianza** y localizar el salto de varianza que introduce el EA (§5).

---

### 3. Modelo de Merton (jump-diffusion)

#### 3.1 SDE

Bajo $\mathbb{Q}$ (§1.4):

$$
\frac{dS_t}{S_{t^-}}=(r-q-\lambda k)\,dt+\sigma\,dW_t+(e^{J}-1)\,dN_t,\qquad J\sim\mathcal{N}(\mu_J,\sigma_J^2).
$$

**Parámetros (4):** $\sigma$ (volatilidad de difusión), $\lambda$ (intensidad anual de saltos), $\mu_J$ (media del salto en log), $\sigma_J$ (desviación del salto).

#### 3.2 Solución y momentos

$$
\ln\frac{S_T}{S_0}=\Big(r-q-\tfrac12\sigma^2-\lambda k\Big)T+\sigma W_T+\sum_{i=1}^{N_T}J_i .
$$

$$
\mathbb{E}[R]=\Big(r-q-\tfrac12\sigma^2-\lambda k+\lambda\mu_J\Big)T,\qquad
\operatorname{Var}[R]=\big(\sigma^2+\lambda(\mu_J^2+\sigma_J^2)\big)T,
$$

y cumulantes de orden $n\ge3$: $\kappa_n=\lambda T\,\mathbb{E}[J^n]$. La **curtosis en exceso** es $\lambda T\,\mathbb{E}[J^4]/\operatorname{Var}[R]^2$ y decae como $1/(\lambda T)$: para ventanas de earnings con $\lambda T=O(1)$, las colas son muy gruesas.

#### 3.3 Función característica

$$
\varphi_M(u)=\exp\Big\{iu\big[\ln S_0+(r-q-\tfrac12\sigma^2-\lambda k)T\big]-\tfrac12\sigma^2u^2T+\lambda T\big(e^{iu\mu_J-\frac12\sigma_J^2u^2}-1\big)\Big\}.
$$

#### 3.4 Precio en serie de Merton

Condicionando en $N_T=n$ el log-precio es normal, de modo que

$$
C_M=\sum_{n=0}^{\infty}\frac{e^{-\lambda'T}(\lambda'T)^n}{n!}\;C_{BS}\big(S_0,K,T,\sigma_n,r_n,q\big),
$$

$$
\lambda'=\lambda(1+k),\qquad \sigma_n^2=\sigma^2+\frac{n\sigma_J^2}{T},\qquad r_n=r-\lambda k+\frac{n\ln(1+k)}{T}.
$$

Es una **mezcla de Poisson de modelos Black–Scholes**: con $\lambda T\approx1$ pocas iteraciones ($n\le 15$) bastan y da una validación independiente del método de Fourier.

#### 3.5 Fair strike bajo Merton

Con $\mathbb{E}[J]=\mu_J$ y $\mathbb{E}[J^2]=\mu_J^2+\sigma_J^2$:

$$
K^{QV}_M=\sigma^2+\lambda\,(\mu_J^2+\sigma_J^2),\qquad
K^{\text{repl}}_M=\sigma^2+2\lambda\big(k-\mu_J\big)=\sigma^2+2\lambda\big(e^{\mu_J+\frac12\sigma_J^2}-1-\mu_J\big).
$$

**Versión discreta** (paso $\Delta=1/252$, retornos iid):

$$
K^{disc}_M=\sigma^2+\lambda(\mu_J^2+\sigma_J^2)+\Delta\Big(r-q-\tfrac12\sigma^2-\lambda(k-\mu_J)\Big)^2 .
$$

El último término es un efecto de deriva de segundo orden ($\Delta\ll1$).

**Descomposición del riesgo de salto puro.** Definiendo $K^{BS}=\sigma^2$:

$$
\underbrace{K^{\text{repl}}_M-\sigma^2}_{\text{aporte del salto en el strike de mercado}}=2\lambda\big(k-\mu_J\big),\qquad
\text{Sesgo}=K^{\text{repl}}_M-K^{QV}_M=2\lambda\,\mathbb{E}\big[e^J-1-J-\tfrac12J^2\big].
$$

---

### 4. Modelo de Heston y modelo de Bates (SV + saltos)

#### 4.1 Heston (1993)

$$
\begin{aligned}
\frac{dS_t}{S_t}&=(r-q)\,dt+\sqrt{v_t}\,dW_t^{S},\\
dv_t&=\kappa(\theta-v_t)\,dt+\xi\sqrt{v_t}\,dW_t^{v},\qquad d\langle W^S,W^v\rangle_t=\rho\,dt .
\end{aligned}
$$

**Parámetros (5):** $v_0,\kappa,\theta,\xi,\rho$. (Bajo $\mathbb{Q}$, $\kappa=\kappa^{\mathbb{P}}+\lambda_v$ y $\theta=\kappa^{\mathbb{P}}\theta^{\mathbb{P}}/\kappa$ si la prima de riesgo de varianza es $\lambda_v v$.)

* **Condición de Feller:** $2\kappa\theta\ge\xi^2$ garantiza $v_t>0$. Las calibraciones a opciones de corto plazo frecuentemente la violan; es aceptable si el esquema de simulación maneja $v=0$ (QE, §8).
* **Esperanza de la varianza instantánea:** $\mathbb{E}[v_t]=\theta+(v_0-\theta)e^{-\kappa t}$.
* **Strike bajo Heston:**

$$
K^{QV}_H=K^{\text{repl}}_H=\frac1T\,\mathbb{E}\!\int_0^Tv_t\,dt=\theta+\frac{(v_0-\theta)\big(1-e^{-\kappa T}\big)}{\kappa T}.
$$

#### 4.2 Bates (1996)

$$
\begin{aligned}
\frac{dS_t}{S_{t^-}}&=(r-q-\lambda k)\,dt+\sqrt{v_t}\,dW_t^{S}+(e^{J}-1)\,dN_t,\\
dv_t&=\kappa(\theta-v_t)\,dt+\xi\sqrt{v_t}\,dW_t^{v},\qquad d\langle W^S,W^v\rangle_t=\rho\,dt,
\end{aligned}
$$

con $N_t$ y $J\sim\mathcal{N}(\mu_J,\sigma_J^2)$ independientes de $(W^S,W^v)$. **Parámetros (8):** $v_0,\kappa,\theta,\xi,\rho,\lambda,\mu_J,\sigma_J$.

#### 4.3 Función característica (formulación estable, "Little Heston Trap")

$$
\varphi_B(u)=\exp\Big\{iu\big[\ln S_0+(r-q)T\big]+C(u,T)+D(u,T)\,v_0+\lambda T\big(e^{iu\mu_J-\frac12\sigma_J^2u^2}-1-iuk\big)\Big\},
$$

con

$$
\begin{aligned}
b&=\kappa-\rho\,\xi\,iu,\qquad d=\sqrt{b^2+\xi^2(iu+u^2)},\qquad g=\frac{b-d}{b+d},\\
D&=\frac{b-d}{\xi^2}\cdot\frac{1-e^{-dT}}{1-g\,e^{-dT}},\\
C&=\frac{\kappa\theta}{\xi^2}\Big[(b-d)\,T-2\ln\frac{1-g\,e^{-dT}}{1-g}\Big].
\end{aligned}
$$

Esta forma (Albrecher, Mayer, Schoutens y Tistaert, 2007) evita las discontinuidades de la rama del logaritmo del formulismo original de Heston. Proviene de resolver el sistema de Riccati

$$
\partial_\tau D=\tfrac12\xi^2D^2+(\rho\xi iu-\kappa)D-\tfrac12(u^2+iu),\qquad
\partial_\tau C=\kappa\theta D+iu(r-q),\qquad D(0)=C(0)=0 .
$$

#### 4.4 Fair strike bajo Bates

Como la parte de difusión y la parte de salto son independientes,

$$
K^{QV}_B=\underbrace{\theta+\frac{(v_0-\theta)(1-e^{-\kappa T})}{\kappa T}}_{\text{volatilidad estocástica}}+\underbrace{\lambda\,(\mu_J^2+\sigma_J^2)}_{\text{saltos (QV)}},
$$

$$
K^{\text{repl}}_B=\theta+\frac{(v_0-\theta)(1-e^{-\kappa T})}{\kappa T}+2\lambda\big(k-\mu_J\big).
$$

(Demostración de $K^{\text{repl}}_B$: $\mathbb{E}\int dS/S=(r-q)T$ y $\mathbb{E}[\ln S_T/S_0]=(r-q-\lambda k+\lambda\mu_J)T-\tfrac12\mathbb{E}\int v\,dt$, de donde $\tfrac2T\mathbb{E}[\int dS/S-\ln S_T/S_0]=\tfrac1T\mathbb{E}\int v\,dt+2\lambda(k-\mu_J)$.)

**Prueba de consistencia obligatoria:** la integral de replicación de §2.3 evaluada con precios de opciones del modelo (Fourier) debe reproducir $K^{\text{repl}}$ analítico. Es un test unitario de todo el pipeline (verificado numéricamente con error relativo $<10^{-6}$).

**Descomposición por componentes:** dado que $K^{\text{repl}}_B-K^{\text{repl}}_M$ depende de que $\sigma^2$ sea reemplazado por el término de volatilidad estocástica, el experimento anidado recomendado es

$$
\text{BS}\ \to\ \text{Merton (+saltos)}\ \to\ \text{Heston (+SV)}\ \to\ \text{Bates (+SV+saltos)},
$$

reportando, para cada evento, la contribución de cada componente al strike y su cambio marginal.

#### 4.5 Qué aporta cada componente en near-expiry

* **Saltos:** generan skew/curvatura que *no desaparece* cuando $T\to0$ (el smile en wings se hace más pronunciado con vencimientos cortos). Son el motor del smile en 5–10 días.
* **Volatilidad estocástica con $\rho<0$:** aporta skew de nivel casi constante en $T$ y, sobre todo, la **dinámica** de la varianza (mean-reversion) que gobierna la estructura a plazo. En un solo vencimiento corto, $(\kappa,\theta,v_0)$ están débilmente identificados (§7).
* **Limitación conocida:** Bates homogéneo produce curvas de varianza esperada *monótonas* en $T$ (regresan suavemente a $\theta$), por lo que no reproduce un escalón por evento.

---

### 5. Modelar el evento de earnings y el vol crush

#### 5.1 Descomposición de varianza total

Sea $w(T)=\sigma_{imp}^2(T)\,T$ la varianza total implícita ATM. Si el vencimiento cruza el anuncio y la volatilidad base es $\sigma_b^2$ (en tiempo de negociación):

$$
w(T)=\sigma_b^2\,T+\varepsilon,\qquad \varepsilon=\mathbb{E}^{\mathbb{Q}}\big[J_e^2\big]\ \ (\text{varianza de evento}).
$$

**Extracción model-free de $\varepsilon$ con dos vencimientos que cruzan el anuncio** ($t_e<T_1<T_2$) y $\sigma_b$ constante:

$$
\sigma_b^2=\frac{\sigma_{imp}^2(T_2)\,T_2-\sigma_{imp}^2(T_1)\,T_1}{T_2-T_1},\qquad
\varepsilon=\sigma_{imp}^2(T_1)\,T_1-\sigma_b^2\,T_1 .
$$

Si existe un vencimiento previo al anuncio $T_0<t_e$, se usa $\sigma_b^2\approx\sigma_{imp}^2(T_0)$ y $\varepsilon=w(T_1)-\sigma_b^2T_1$. Conceptualmente, esta separación entre volatilidad "normal" e incertidumbre de EA es la que formalizan Dubinsky et al. (2019) con modelos reducidos y estimadores más completos; las identidades de arriba son la versión simple, model-free.

**Movimiento implícito (*implied move*).** Con retorno aproximadamente normal de varianza total $w$, el straddle ATM cuesta $\approx\sqrt{2/\pi}\,S_0\sqrt{w}\approx0.80\,S_0\sqrt{w}$, de modo que $\text{straddle}/S_0\approx\mathbb{E}|R|$.

#### 5.2 Modelo de evento

Modelo general: base (Merton o Bates) más un salto de evento en fecha conocida $t_e$,

$$
\ln S_T=\ln S_0+\big(\text{proceso base}\big)+J_e\,\mathbf{1}_{\{t_e\le T\}}-\ln\mathbb{E}\big[e^{J_e}\big]\,\mathbf{1}_{\{t_e\le T\}},
$$

y, por independencia entre base y evento, la función característica se factoriza:

$$
\varphi(u)=\varphi_{base}(u;T)\cdot\psi_e(u)\,e^{-iu\ln\mathbb{E}[e^{J_e}]},\qquad \psi_e(u)=\mathbb{E}\big[e^{iuJ_e}\big].
$$

Opciones para la ley del salto de evento $J_e$:

| Opción | $J_e$ | Implicancia |
|--------|-------|-------------|
| (i) Gaussiano determinista | $\mathcal{N}(\mu_e,\sigma_e^2)$ | **Sin smile:** $\ln S_T$ es normal con varianza $\sigma^2T+\sigma_e^2$; equivale a Black–Scholes con mayor varianza total. Sirve como *benchmark*, no como modelo del smile. |
| (ii) Poisson de ventana | $N_e\sim\text{Poisson}(\Lambda_e)$, $\Lambda_e=\lambda_e\Delta_e\approx1$ | Mezcla de gaussianas. Es la interpretación de "intensidad alta exactamente en la fecha del anuncio". Genera colas gruesas y smile. |
| (iii) Bernoulli | $N_e\in\{0,1\}$ con prob. $p$ | Permite evento que "no ocurre" (retraso del anuncio). |
| (iv) Doble exponencial | Kou (2002) | Asimetría marcada arriba/abajo. |
| (v) Mezcla bimodal | dos componentes $\pm$ | Reproduce la **bimodalidad** de la densidad neutral al riesgo antes de EA (Kachhara, Markin y Singh, 2023) y las **curvas de volatilidad implícita cóncavas** (Alexiou, Goyal, Kostakis y Rompolis, 2025), que Merton/Bates de una sola componente gaussiana no pueden generar. |

Para la opción (ii), el factor multiplicativo de la función característica es

$$
\psi_e^{(ii)}(u)=\exp\Big\{\Lambda_e\big(e^{iu\mu_e-\frac12\sigma_e^2u^2}-1-iuk_e\big)\Big\},\qquad k_e=e^{\mu_e+\frac12\sigma_e^2}-1 .
$$

**Consecuencia para $K_{var}$:** la contribución del evento al strike es $\varepsilon/T$ en QV o $2\Lambda_e(k_e-\mu_e)/T$ en replicación, con el mismo sesgo de tercer orden de §2.5.

#### 5.3 Parámetros dependientes del tiempo

Para permitir un escalón en la estructura a plazo, se usa intensidad y variancia **por tramos**:

$$
\lambda(t)=\lambda_{pre}\,\mathbf{1}_{t<t_e^-}+\lambda_e\,\mathbf{1}_{|t-t_e|\le\epsilon}+\lambda_{post}\,\mathbf{1}_{t>t_e^+},\qquad
\theta(t),\ \kappa(t)\ \text{análogos}.
$$

Con parámetros constantes por tramos, el sistema de Riccati (§4.3) se **propaga tramo a tramo** en $\tau=T-t$ con la condición terminal de cada tramo (fórmula cerrada con $D(0)\ne0$ o integración RK4). La integral de intensidad $\Lambda(T)=\int_0^T\lambda(t)\,dt$ sustituye a $\lambda T$ en el término de salto (si $\mu_J,\sigma_J$ son constantes).

#### 5.4 Vol crush: definición y lectura en el modelo

Para la opción del mismo vencimiento antes y después del anuncio, y con volatilidad base constante:

$$
\sigma^2_{pre}=\sigma_b^2+\frac{\varepsilon}{T},\qquad \sigma^2_{post}\approx\sigma_b^2,\qquad
CR=\frac{\sigma_{post}}{\sigma_{pre}}\approx\sqrt{\frac{\sigma_b^2}{\sigma_b^2+\varepsilon/T}} .
$$

Mecanismos en los modelos:

1. **Eliminación del salto de evento** ($\lambda(t)$ cae a $\lambda_{post}$): es el mecanismo dominante y natural en Merton/Bates por tramos.
2. **Reversión de $v_t$** hacia $\theta$ (Bates homogéneo): solo produce una caída gradual, no instantánea.
3. **Reset de varianza en $t_e$:** $v_{t_e^+}=c\cdot v_{t_e^-}$ con $c\in(0,1]$, o $\theta(t)$ por tramos. Modelos con saltos *positivos* en varianza (SVJJ; Duffie et al., 2000; Eraker, Johannes y Polson, 2003) generan alzas de varianza, no colapsos, y no son adecuados para el crush sin modificaciones.

**Impacto sobre el swap.** Por §2.1, la posición larga comprada antes del EA gana $R_e^2$ (aporta a $RV_{0,t}$) y pierde la diferencia entre $K_{t,T}$ pre y post: el vol crush ocurre en el término forward $K_{t,T}$, mientras el salto realizado ocurre en $RV_{0,t}$. La **participación del evento** $\eta=\frac{R_e^2}{\sum_iR_i^2}$ es un estadístico clave de cada trayectoria simulada.

---

### 6. Valoración de opciones por transformadas de Fourier

#### 6.1 Gil-Pelaez / Heston (1993)

$$
C(K)=S_0e^{-qT}P_1-Ke^{-rT}P_2,
$$

$$
P_2=\frac12+\frac1\pi\int_0^\infty\operatorname{Re}\Big[\frac{e^{-iu\ln K}\varphi(u)}{iu}\Big]du,\qquad
P_1=\frac12+\frac1\pi\int_0^\infty\operatorname{Re}\Big[\frac{e^{-iu\ln K}\varphi(u-i)}{iu\,\varphi(-i)}\Big]du .
$$

#### 6.2 Carr–Madan (1999): FFT con amortiguamiento

Con factor $\alpha>0$ de amortiguamiento y $k=\ln K$:

$$
C(K)=\frac{e^{-\alpha k}}{\pi}\int_0^\infty e^{-i\nu k}\,\psi(\nu)\,d\nu,\qquad
\psi(\nu)=\frac{e^{-rT}\,\varphi\big(\nu-(\alpha+1)i\big)}{\alpha^2+\alpha-\nu^2+i(2\alpha+1)\nu}.
$$

Se evalúa con FFT para toda una malla de strikes por vencimiento. Existe además la formulación de Lewis (2001) y el análisis de error de Lee (2004).

#### 6.3 Aspectos numéricos para near-expiry

* Para $T\approx1/50$ años el decaimiento de $\varphi(u)$ en $u$ puede ser lento (más aún en Bates por $\xi$); use límites de integración amplios y **verifique convergencia** doblando el límite superior.
* Para una calibración rápida, **cuadratura de Gauss–Legendre fija** sobre la integral de Lewis/Gil-Pelaez, común a todos los strikes de un vencimiento (Lee, 2004; Cui, del Baño Rollin y Germano, 2017, para Heston).
* Precios de puts OTM muy pequeños: evitar cancelación numérica al usar paridad put-call en strikes muy bajos; integrar $P(K)$ vía Fourier directamente o usar Black–Scholes sobre la vol interpolada.

---

### 7. Calibración

#### 7.1 Datos requeridos

Por evento y fecha de observación:

* Cadena completa de calls y puts a 2–3 vencimientos: uno **pre-EA** (si existe), los **near-expiry que cruzan** el anuncio (5, 7, 10 días) y uno posterior (para anclar $\kappa,\theta$).
* Bid/ask, open interest, volumen; spot y ex-dividendos dentro de la ventana; tasa libre de riesgo y costo de préstamo.
* Fecha/hora exacta del anuncio (*before market open* vs. *after market close*), pues determina **qué retorno diario** contiene el salto.

#### 7.2 Limpieza y preparación

1. **Forward implícito** por paridad put-call en strikes ATM: $F_0=K+e^{rT}(C-P)$, promediado en varios strikes; deduce $q$ efectivo (dividendos + *borrow*).
2. **Opciones americanas.** Las opciones sobre acciones individuales son americanas. Convertirlas a equivalentes europeas (binomial/Bjerksund–Stensland con dividendos discretos) o usar solo **OTM** (donde la prima de ejercicio anticipado es pequeña) para inferir la vol implícita europea.
3. **Filtros:** descartar bid$=0$, spreads relativos excesivos, precios que violen límites de no-arbitraje, y strikes fuera de un rango de delta razonable.
4. **No arbitraje:** comprobar convexidad en $K$ (mariposas) y monotonía calendario; corregir con SVI/Fengler antes de calibrar.
5. Trabajar en **tiempo de negociación** (días hábiles/252): en tiempo calendario los fines de semana diluyen la varianza y distorsionan $\sigma_b$ y $\varepsilon$.

#### 7.3 Función objetivo

Sean $\Theta$ los parámetros y $i$ el índice de opción:

$$
\hat\Theta=\arg\min_{\Theta}\ \sum_{i=1}^{M}w_i\Big(\sigma^{model}_i(\Theta)-\sigma^{mkt}_i\Big)^2+\alpha\,\big\|\Theta-\Theta_{prior}\big\|^2 .
$$

* Equivalente en precios: $\sum_iw_i'\big(C_i^{model}-C_i^{mkt}\big)^2$ con $w_i'=1/\text{vega}_i^2$ (o $1/(\text{ask}-\text{bid})^2$), porque $\Delta C\approx\text{vega}\cdot\Delta\sigma$. **El RMSE de precios sin ponderar** sobreajusta las opciones caras (ITM/largas) y falla en wings baratos, que son los que determinan el sesgo por saltos.
* $w_i$: inverso del spread bid–ask o de la varianza del error de medición.
* El término $\alpha\|\Theta-\Theta_{prior}\|^2$ (regularización de Tikhonov) estabiliza parámetros débilmente identificados ($\kappa,\theta$).

#### 7.4 Restricciones y optimizador

| Parámetro | Restricción sugerida |
|-----------|----------------------|
| $\sigma,\ v_0,\ \theta,\ \xi$ | $>0$ (cotas superiores razonables) |
| $\kappa$ | $[0.1,\,20]$ |
| $\rho$ | $(-1,\,1)$ |
| $\lambda$ | $[0,\ \lambda_{max}]$ con $\lambda_{max}$ acorde a $\lambda T\lesssim5$ |
| $\mu_J$ | $[-0.5,\,0.5]$ |
| $\sigma_J$ | $(0,\,1]$ |

* **Global → local:** `scipy.optimize.differential_evolution` (o *multi-start*), seguido de `scipy.optimize.least_squares` (método `trf`, con *bounds*).
* Monitorear múltiples mínimos locales; reportar la **dispersión de parámetros** entre inicializaciones.
* Verificar la condición de **martingala** $\varphi(-i)=F_0$ tras cada calibración.

#### 7.5 Identificabilidad

* Con **un solo vencimiento**, $(v_0,\theta,\kappa)$ son casi colineales: en $T\to0$ el modelo solo depende de $v_0$ y de la combinación $(\rho,\xi)$; $\kappa,\theta$ influyen apenas. Estrategias: calibrar **conjuntamente** varios vencimientos, fijar $\kappa$ y $\theta$ con vencimientos largos, o regularizar hacia priors históricos.
* Bates tiene 8 parámetros y Merton 4: comparar con criterios penalizados (AIC/BIC) y **validación fuera de muestra** (calibrar en strikes impares y testear en pares; o calibrar en $t-1$ y evaluar en $t$).
* Los saltos y la vol estocástica compiten para explicar el mismo smile: monitorear la correlación entre estimadores.

---

### 8. Simulación Monte Carlo

#### 8.1 Merton: esquema exacto

Para paso $\Delta$ (típicamente 1 día = $1/252$):

$$
\ln S_{t+\Delta}=\ln S_t+\Big(r-q-\tfrac12\sigma^2-\lambda k\Big)\Delta+\sigma\sqrt\Delta\,Z_1+N\mu_J+\sigma_J\sqrt N\,Z_2,\qquad N\sim\text{Poisson}(\lambda\Delta),
$$

con $Z_1,Z_2\sim\mathcal{N}(0,1)$ independientes. La suma de $N$ saltos normales es $\mathcal{N}(N\mu_J,N\sigma_J^2)$, por lo que **no hay error de discretización** (el esquema es exacto en distribución).

#### 8.2 Bates: esquema QE de Andersen (2008)

**Varianza** (Quadratic-Exponential). Dado $v_t$:

$$
m=\theta+(v_t-\theta)e^{-\kappa\Delta},\qquad
s^2=\frac{v_t\xi^2e^{-\kappa\Delta}(1-e^{-\kappa\Delta})}{\kappa}+\frac{\theta\xi^2(1-e^{-\kappa\Delta})^2}{2\kappa},\qquad \psi=\frac{s^2}{m^2}.
$$

* Si $\psi\le\psi_c$ (típicamente $\psi_c=1.5$): $v_{t+\Delta}=a(b+Z_v)^2$, con $b^2=\frac2\psi-1+\sqrt{\frac2\psi}\sqrt{\frac2\psi-1}$ y $a=\frac{m}{1+b^2}$.
* Si $\psi>\psi_c$: con $p=\frac{\psi-1}{\psi+1}$, $\beta=\frac{1-p}{m}$ y $U\sim\mathcal{U}(0,1)$: $v_{t+\Delta}=0$ si $U\le p$; en otro caso $v_{t+\Delta}=\frac1\beta\ln\frac{1-p}{1-U}$.

**Log-precio** (con $\gamma_1=\gamma_2=\tfrac12$):

$$
\ln S_{t+\Delta}=\ln S_t+(r-q-\lambda k)\Delta+K_0+K_1v_t+K_2v_{t+\Delta}+\sqrt{K_3v_t+K_4v_{t+\Delta}}\;Z+N\mu_J+\sigma_J\sqrt N\,Z_2,
$$

$$
K_0=-\frac{\rho\kappa\theta\Delta}{\xi},\quad
K_1=\gamma_1\Delta\Big(\frac{\kappa\rho}{\xi}-\frac12\Big)-\frac\rho\xi,\quad
K_2=\gamma_2\Delta\Big(\frac{\kappa\rho}{\xi}-\frac12\Big)+\frac\rho\xi,\quad
K_3=\gamma_1\Delta(1-\rho^2),\quad K_4=\gamma_2\Delta(1-\rho^2).
$$

Alternativa más simple, aunque sesgada: **full truncation Euler** (Lord, Koekkoek y van Dijk, 2010). Para $\Delta=1/252$ y $\xi$ moderado, QE es notablemente más preciso cerca de $v=0$.

#### 8.3 Ventana de earnings

Se simula sobre la malla diaria hasta el vencimiento ($n$ pasos) con parámetros por tramos:

1. **Pre-EA:** $\lambda_{pre}$ (bajo) y $v_0$ calibrado.
2. **Día del EA** (índice de la fila que contiene el salto según BMO/AMC): $\Lambda_e=\lambda_e\Delta_e\approx1$ (p. ej. $\lambda_e=252$ en un paso diario) **o** salto único con distribución elegida en §5.2.
3. **Post-EA:** $\lambda_{post}$ y, opcionalmente, reset de varianza $v\leftarrow c\,v$ (vol crush).

El compensador del salto en el paso de evento es $-\Lambda_ek_e$ (no $-\lambda k\Delta$ con la intensidad base).

#### 8.4 Estadísticos por trayectoria y estimadores

* $RV^{(m)}=\frac{252}{n}\sum_{i=1}^{n}\big(\ln S^{(m)}_{t_i}/S^{(m)}_{t_{i-1}}\big)^2$ (misma convención que el contrato).
* $\hat K=\frac1M\sum_mRV^{(m)}$ con error estándar $\hat s/\sqrt M$. Contrastar $\hat K$ con: $K^{disc}$ (Merton, exacto), $K^{QV}$ y $K^{\text{repl}}$ analíticos.
* $P\&L^{(m)}=N_{var}\big(RV^{(m)}-K_{var}\big)$: media, desviación, asimetría, cuantiles al 1/5%.
* $\eta^{(m)}=R_{e}^2/\sum_iR_i^2$ (participación del día del evento) y $RV^{(m)}_{ex}$ (RV excluyendo el día del EA).
* Chequeo de **martingala**: $\mathbb{E}[S_T]e^{-(r-q)T}\approx1$.

#### 8.5 Reducción de varianza y semillas

* Variables antitéticas en $Z$; *control variate* usando una variable con esperanza conocida (por ejemplo, $\int v_t\,dt$ aproximada con la regla trapezoidal, cuya esperanza es la de §4.1).
* $M\ge10^5$–$10^6$ trayectorias; fijar semilla (`numpy.random.default_rng(seed)`) y reportar el intervalo de confianza.

---

## Bibliografía

> **Nota:** las notas de research de bancos (Goldman Sachs, J.P. Morgan) son documentos de *sell-side*; su disponibilidad en línea varía. Verifique DOI/páginas y accesibilidad antes de someter el paper.

### Research de bancos (variance swaps y volatilidad)

* Demeterfi, K., Derman, E., Kamal, M. y Zou, J. (1999). *More Than You Ever Wanted to Know About Volatility Swaps*. Goldman Sachs Quantitative Strategies Research Notes, marzo de 1999. Versión publicada: *A Guide to Volatility and Variance Swaps*, **Journal of Derivatives**, 6(4), 9–32.
* Derman, E., Kamal, M., Kani, I., McClure, J., Pirasteh, C. y Zou, J. (1998). *Investing in Volatility*. Futures & Options World. (Goldman Sachs).
* Bossu, S., Strasser, E. y Guichard, R. (2005). *Just What You Need to Know About Variance Swaps*. J.P. Morgan, Equity Derivatives Investor Marketing / Quantitative Research and Development, Londres, mayo de 2005.
* Allen, P., Einchcomb, S. y Granger, N. (2006). *Variance Swaps*. J.P. Morgan Securities Ltd., European Equity Derivatives Research.
* Bergomi, L. (2016). *Stochastic Volatility Modeling*. Chapman & Hall/CRC. (Práctica cuantitativa de Société Générale).

### Modelos y replicación

* Black, F. y Scholes, M. (1973). The Pricing of Options and Corporate Liabilities. **Journal of Political Economy**, 81(3), 637–654.
* Merton, R. C. (1976). Option Pricing When Underlying Stock Returns Are Discontinuous. **Journal of Financial Economics**, 3(1–2), 125–144.
* Heston, S. L. (1993). A Closed-Form Solution for Options with Stochastic Volatility with Applications to Bond and Currency Options. **Review of Financial Studies**, 6(2), 327–343.
* Bates, D. S. (1996). Jumps and Stochastic Volatility: Exchange Rate Processes Implicit in Deutsche Mark Options. **Review of Financial Studies**, 9(1), 69–107.
* Bates, D. S. (2000). Post-'87 Crash Fears in the S&P 500 Futures Option Market. **Journal of Econometrics**, 94(1–2), 181–238.
* Duffie, D., Pan, J. y Singleton, K. (2000). Transform Analysis and Asset Pricing for Affine Jump-Diffusions. **Econometrica**, 68(6), 1343–1376.
* Eraker, B., Johannes, M. y Polson, N. (2003). The Impact of Jumps in Volatility and Returns. **Journal of Finance**, 58(3), 1269–1300.
* Pan, J. (2002). The Jump-Risk Premia Implicit in Options: Evidence from an Integrated Time-Series Study. **Journal of Financial Economics**, 63(1), 3–50.
* Kou, S. G. (2002). A Jump-Diffusion Model for Option Pricing. **Management Science**, 48(8), 1086–1101.
* Breeden, D. y Litzenberger, R. (1978). Prices of State-Contingent Claims Implicit in Option Prices. **Journal of Business**, 51(4), 621–651.
* Neuberger, A. (1994). The Log Contract: A New Instrument to Hedge Volatility. **Journal of Portfolio Management**, invierno de 1994.
* Dupire, B. (1993). Model Art. **Risk**, septiembre de 1993.
* Carr, P. y Madan, D. (1998). Towards a Theory of Volatility Trading. En *Volatility: New Estimation Techniques for Pricing Derivatives*, Risk Books, 417–427.
* Britten-Jones, M. y Neuberger, A. (2000). Option Prices, Implied Price Processes, and Stochastic Volatility. **Journal of Finance**, 55(2), 839–866.
* Bakshi, G., Cao, C. y Chen, Z. (1997). Empirical Performance of Alternative Option Pricing Models. **Journal of Finance**, 52(5), 2003–2049.
* Bakshi, G., Kapadia, N. y Madan, D. (2003). Stock Return Characteristics, Skew Laws, and the Differential Pricing of Individual Equity Options. **Review of Financial Studies**, 16(1), 101–143.
* Cont, R. y Tankov, P. (2004). *Financial Modelling with Jump Processes*. Chapman & Hall/CRC.
* Gatheral, J. (2006). *The Volatility Surface: A Practitioner's Guide*. Wiley.

### Variance/volatility swaps con saltos y discretización

* Broadie, M. y Jain, A. (2008). The Effect of Jumps and Discrete Sampling on Volatility and Variance Swaps. **International Journal of Theoretical and Applied Finance**, 11(8), 761–797.
* Broadie, M. y Jain, A. (2008). Pricing and Hedging Volatility Derivatives. **Journal of Derivatives**, 15(3), 7–24.
* Carr, P. y Lee, R. (2007). Realized Volatility and Variance: Options via Swaps. **Risk**, 20(5), 76–83.
* Carr, P. y Lee, R. (2009). Volatility Derivatives. **Annual Review of Financial Economics**, 1, 319–339.
* Carr, P. y Lewis, K. (2004). Corridor Variance Swaps. **Risk**, 17(2), 67–72.
* Jarrow, R., Kchia, Y., Larsson, M. y Protter, P. (2013). Discretely Sampled Variance and Volatility Swaps Versus Their Continuous Approximations. **Finance and Stochastics**, 17(2), 305–324.
* Jiang, G. J. y Tian, Y. S. (2005). The Model-Free Implied Volatility and Its Information Content. **Review of Financial Studies**, 18(4), 1305–1342.
* Jiang, G. J. y Tian, Y. S. (2007). Extracting Model-Free Volatility from Option Prices: An Examination of the VIX Methodology. **Journal of Derivatives**, 14(3), 35–60.
* Barndorff-Nielsen, O. y Shephard, N. (2004). Power and Bipower Variation with Stochastic Volatility and Jumps. **Journal of Financial Econometrics**, 2(1), 1–37.
* Carr, P. y Wu, L. (2009). Variance Risk Premiums. **Review of Financial Studies**, 22(3), 1311–1341.
* Bollerslev, T., Tauchen, G. y Zhou, H. (2009). Expected Stock Returns and Variance Risk Premia. **Review of Financial Studies**, 22(11), 4463–4492.

### Earnings announcements y opciones

* Dubinsky, A., Johannes, M., Kaeck, A. y Seeger, N. J. (2019). Option Pricing of Earnings Announcement Risks. **Review of Financial Studies**, 32(2), 646–687.
* Barth, M. E. y So, E. C. (2014). Non-Diversifiable Volatility Risk and Risk Premiums at Earnings Announcements. **The Accounting Review**, 89(5), 1579–1607.
* Gao, C., Xing, Y. y Zhang, X. (2018). Anticipating Uncertainty: Straddles Around Earnings Announcements. **Journal of Financial and Quantitative Analysis**, 53(6), 2587–2617.
* Alexiou, L., Goyal, A., Kostakis, A. y Rompolis, L. (2025). Pricing Event Risk: Evidence from Concave Implied Volatility Curves. **Review of Finance**, 29(4), 963–1007.
* Kachhara, D., Markin, J. K. E. y Singh, A. (2023). *Option Smile Volatility and Implied Probabilities: Implications of Concavity in IV Curves*. arXiv:2307.15718.
* Patell, J. y Wolfson, M. (1979). Anticipated Information Releases Reflected in Call Option Prices. **Journal of Accounting and Economics**, 1(2), 117–140.
* Beaver, W. (1968). The Information Content of Annual Earnings Announcements. **Journal of Accounting Research**, 6 (Supplement), 67–92.
* Isakov, D. y Périgon, C. (2001). Evolution of Market Uncertainty Around Earnings Announcements. **Journal of Banking & Finance**, 25(9), 1769–1788.

### Métodos numéricos, simulación y calibración

* Carr, P. y Madan, D. (1999). Option Valuation Using the Fast Fourier Transform. **Journal of Computational Finance**, 2(4), 61–73.
* Lewis, A. L. (2001). *A Simple Option Formula for General Jump-Diffusion and Other Exponential Lévy Processes*. SSRN.
* Lee, R. W. (2004). Option Pricing by Transform Methods: Extensions, Unification and Error Control. **Journal of Computational Finance**, 7(3), 51–86.
* Albrecher, H., Mayer, P., Schoutens, W. y Tistaert, J. (2007). The Little Heston Trap. **Wilmott Magazine**, enero de 2007, 83–92.
* Andersen, L. (2008). Simple and Efficient Simulation of the Heston Stochastic Volatility Model. **Journal of Computational Finance**, 11(3), 1–42.
* Lord, R., Koekkoek, R. y van Dijk, D. (2010). A Comparison of Biased Simulation Schemes for Stochastic Volatility Models. **Quantitative Finance**, 10(2), 177–194.
* Cui, Y., del Baño Rollin, S. y Germano, G. (2017). Full and Fast Calibration of the Heston Stochastic Volatility Model. **European Journal of Operational Research**, 263(2), 625–638.
* Gatheral, J. y Jacquier, A. (2014). Arbitrage-Free SVI Volatility Surfaces. **Quantitative Finance**, 14(1), 59–71.
* Fengler, M. R. (2009). Arbitrage-Free Smoothing of the Implied Volatility Surface. **Quantitative Finance**, 9(4), 417–428.

---

**Licencia y uso:** documento de investigación con fines académicos; no constituye asesoramiento financiero.
