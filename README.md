# Variance Swaps alrededor de Earnings: Marco Teórico con Merton (1976) y Bates (1996)

> Valoración y replicación del *fair strike* de un variance swap bajo saltos y volatilidad estocástica, con aplicación a ventanas de earnings.

## Tabla de contenidos

- [Marco Teórico](#marco-teórico)
  - [1. El variance swap: definición y replicación](#1-el-variance-swap-definición-y-replicación)
  - [2. Modelo de Merton (jump-diffusion)](#2-El-modelo-de-merton-jump-diffusion)
  - [3. Modelo de Heston y modelo de Bates (SV + saltos)](#3-modelo-de-heston-y-modelo-de-bates-sv--saltos)
  - [4. Modelar el evento de earnings y el vol crush](#4-modelar-el-evento-de-earnings-y-el-vol-crush)
  - [5. Valoración de opciones por transformadas de Fourier](#5-valoración-de-opciones-por-transformadas-de-fourier)
  - [6. Calibración](#6-calibración)
  - [7. Simulación Monte Carlo](#7-simulación-monte-carlo)
- [Bibliografía](#bibliografía)

---

## Marco Teórico

**Notación.** $S_t$ precio del subyacente, $r$ tasa libre de riesgo, $q$ dividend yield continuo, $T$ vencimiento (en años de *tiempo de negociación*, es decir, días hábiles/252), $F_0=S_0e^{(r-q)T}$ forward, $X_t=\ln S_t$, $\mathbb{Q}$ medida neutral al riesgo, $\mathbb{P}$ medida física, $t_e$ fecha del anuncio de resultados.

---

# 1. El variance swap: definición y replicación

## Estructura, Mecánica y Liquidación

Un *variance swap* es un contrato forward sobre la varianza anualizada de un activo subyacente específico. Su pago al vencimiento $T$ es lineal respecto a la varianza realizada:

$$\text{Payoff} = (\sigma_{R}^{2} - K_{\text{var}}) \times N_{\text{vol}}$$

Donde $\sigma_{R}^{2}$ denota la varianza realizada del activo subyacente durante el horizonte del contrato, $K_{\text{var}}$ representa el *strike* de varianza implícita y $N_{\text{vol}}$ corresponde al valor nocional expresado en unidades monetarias por punto de varianza.

---

## Derivación Matemática del Strike de Varianza ($K_{\text{var}}$)

El objetivo central de la modelización cuantitativa consiste en determinar el precio de ejercicio justo ($K_{\text{var}}$), el cual, bajo condiciones de ausencia de arbitraje, debe igualar la expectativa matemática de la varianza realizada bajo la medida neutral al riesgo.

### El Proceso de Precios Continuo
Se asume que la trayectoria del precio del subyacente $S_t$ sigue un movimiento browniano geométrico bajo la medida histórica, expresado estocásticamente en términos de su retorno instantáneo:

$$\frac{dS_t}{S_t} = \mu_t dt + \sigma_t dW_t$$

Donde $\mu_t$ denota el *drift*, $\sigma_t$ representa la volatilidad instantánea y $dW_t$ corresponde al incremento de un proceso de Wiener estándar. La varianza realizada anualizada $V$ se formaliza como:

$$V = \frac{1}{T}\int_{0}^{T}\sigma_{t}^{2}dt$$

### Aplicación del Lema de Itô y Transformación Logarítmica
El desafío de esta valoración radica en aislar la varianza instantánea $\sigma_t^2$ utilizando únicamente variables de mercado observables (en este caso, el precio spot del activo). Puesto que $\sigma_t$ no es una variable directamente observable en el mercado, se recurre a una transformación logarítmica mediante la aplicación del Lema de Itô sobre la función $f(S_t) = \ln(S_t)$.

Para aislar este término del proceso del precio del activo, partimos de las derivadas de la función: $f'(S) = \frac{1}{S}$ y $f''(S) = -\frac{1}{S^2}$, las cuales permiten cancelar el término $S^2$ al aplicar el lema de Itô:

$$d(\ln S_t) = f'(S_t)dS_t + \frac{1}{2}f''(S_t)(dS_t)^2$$

Desarrollando el término $(dS_t)^2$ a partir del movimiento browniano geométrico del precio:

$$(dS_t)^2 = (\mu_t S_t \ dt + \sigma_t S_t \ dW_t)^2$$
$$(dS_t)^2 = \mu_t^2 S_t^2 (dt)^2 + 2\mu_t \sigma_t S_t^2 dt \ dW_t + \sigma_t^2 S_t^2 (dW_t)^2$$

Al aplicar las reglas del cálculo estocástico para cancelar los términos de mayor orden:

$$(dt)^2 = 0, \quad dW_t \ dt = 0, \quad (dW_t)^2 = dt$$

Obtenemos el desarrollo diferencial del término cuadrático:

$$(dS_t)^2 = \sigma_t^2 S_t^2 \ dt$$

> **Nota explicativa sobre el cálculo estocástico empleado:**  
> La identidad del cálculo de Itô $(dW_t)^2 = dt$ difiere del cálculo clásico porque los incrementos del movimiento browniano son tan impredecibles y rápidos que su variación cuadrática en un intervalo pequeño no se anula, sino que converge exactamente al paso del tiempo $dt$.

Sustituyendo todo en la ecuación original del lema de Itô:

$$d(\ln S_t) = \frac{1}{S_t}dS_t - \frac{1}{2}\left(\frac{1}{S_t^2}\right)(\sigma_t^2 S_t^2 \ dt)$$

Simplificando la expresión, se obtiene el retorno porcentual menos un ajuste por convexidad de $-\frac{1}{2}\sigma_t^2 dt$:

$$d(\ln S_t) = \frac{1}{S_t}dS_t - \frac{1}{2}\sigma_t^2 \ dt$$

---

Reorganizando la ecuación resultante, se aísla el término de la varianza $\sigma_t^2 dt$. Tras integrar la expresión a lo largo del intervalo $[0, T]$, se demuestra que la varianza realizada se refleja exactamente manteniendo una posición en acciones rebalanceada continuamente (manteniendo siempre $\frac{1}{S_t}$ acciones) y tomando una posición corta estática en un instrumento que pague $\ln\left(\frac{S_T}{S_0}\right)$:

$$V = \frac{2}{T}\left[\int_{0}^{T}\frac{dS_t}{S_t} - \ln\left(\frac{S_T}{S_0}\right)\right]$$
## Esperanza Neutral al Riesgo y el Teorema de Carr-Madan

Para encontrar un instrumento negociable en el mercado que sea equivalente al pago logarítmico, la valoración se traslada al marco de la probabilidad neutral al riesgo ($\mathbb{Q}$). Bajo esta medida, el activo subyacente crece en promedio a la tasa libre de riesgo $r$, permitiendo descontar flujos esperados a valor presente sin sesgos de preferencias de riesgo subjetivas. Así, el *strike* teórico del *variance swap* se define como:

$$K_{\text{var}} = \mathbb{E}^{\mathbb{Q}}\left[ \frac{1}{T}\int_{0}^{T}\sigma_{t}^{2}dt \right]$$

Para eliminar la necesidad de negociar contratos logarítmicos sintéticos, se implementa el teorema de replicación estática desarrollado por Peter Carr y Dilip Madan (1998). Este teorema demuestra que cualquier pago dos veces diferenciable $f(S_T)$ puede descomponerse analíticamente mediante una expansión de Taylor ponderada sobre una cadena infinita de opciones europeas *vanilla* (*calls* y *puts*).

### Desglose analítico del Teorema de Carr-Madan

* **Término de posición constante:** 

$$
f(S_*)
$$

evaluado en un precio de referencia preestablecido $S_*$. Usualmente, $S_*$ se fija como el precio *forward* vigente del activo. Este término representa simplemente un valor constante que no depende del precio futuro $S_T$.

* **Término de exposición lineal:**

$$
f'(S_{\ast})(S_T-S_{\ast})
$$

Este componente representa la parte lineal del payoff, ya que depende directamente de $S_T$. Debido a esta característica, puede replicarse mediante una posición estática en contratos *forward* o futuros sobre el activo subyacente.

* **Componentes convexos de opciones:** Son los términos encargados de representar la curvatura del payoff y están dados por dos integrales.

La primera corresponde a una combinación de opciones *put* europeas:

$$
\int_0^{S_*} f''(K)(K-S_T)^+\, dK
$$

La segunda corresponde a una combinación de opciones *call* europeas:

$$
\int_{S_*}^{\infty} f''(K)(S_T-K)^+\, dK
$$

Esto se debe a que:

$$
(K-S_T)^+
$$

es exactamente el *payoff* de una *put* europea, mientras que:

$$
(S_T-K)^+
$$

es exactamente el *payoff* de una *call* europea.

Por lo tanto, la representación completa de Carr-Madan puede expresarse como:

$$
f(S_T) = f(S_{\ast}) + f'(S_{\ast})(S_T-S_{\ast}) + \int_0^{S_{\ast}} f''(K)(K-S_T)^+\, dK + \int_{S_{\ast}}^{\infty} f''(K)(S_T-K)^+\, dK
$$

En conjunto, la idea es que un payoff complejo puede descomponerse en una posición constante, una exposición lineal mediante *forwards* o futuros y una combinación de *puts* y *calls* con diferentes *strikes* para reproducir la curvatura del payoff.
En conjunto, la idea es que un payoff complejo puede descomponerse en una posición constante, una exposición lineal mediante *forwards* o futuros y una combinación de *puts* y *calls* con diferentes *strikes* para reproducir la curvatura del payoff.

Con $f(S)=\ln S$ se tiene $f''(K)=-1/K^2$, de modo que

$$
\ln\frac{S_T}{S_{\ast}}=\frac{S_T-S_{\ast}}{S_{\ast}}-\int_0^{S_{\ast}}\frac{(K-S_T)^+}{K^2}\, dK-\int_{S_{\ast}}^{\infty}\frac{(S_T-K)^+}{K^2}\, dK .
$$

Al aplicar el operador de esperanza bajo $\mathbb{Q}$ se usan dos hechos: (i) $\mathbb{E}^{\mathbb{Q}}[S_T]=S_0e^{rT}$ y (ii) $\mathbb{E}^{\mathbb{Q}}\!\left[\int_0^T \frac{dS_t}{S_t}\right]=rT$. Además, el valor esperado no descontado de un pago de opción es su prima capitalizada, $e^{rT}P(K)$ o $e^{rT}C(K)$, donde $P(K)$ y $C(K)$ denotan las primas de las opciones *put* y *call* de *strike* $K$. Sustituyendo en (2) se llega a la fórmula de valoración general:

$$
K_{var} = \frac{2}{T} \left[ rT - \left( \frac{S_0}{S_{\ast}} e^{rT} - 1 \right) - \ln\left(\frac{S_{\ast}}{S_0}\right) + e^{rT} \int_0^{S_{\ast}} \frac{P(K)}{K^2} dK + e^{rT} \int_{S_{\ast}}^\infty \frac{C(K)}{K^2} dK \right] \tag{3}
$$

Si se elige $S_{\ast}=F_0=S_0e^{rT}$ (el *forward*), todos los términos no integrales se cancelan y queda la forma limpia de Britten-Jones y Neuberger (2000):

$$
K_{var}=\frac{2}{T}\,e^{rT}\left[\int_0^{F_0}\frac{P(K)}{K^2}\, dK+\int_{F_0}^{\infty}\frac{C(K)}{K^2}\, dK\right] \tag{4}
$$

(Con un dividendo continuo $q$, se reemplaza $rT$ por $(r-q)T$ en el primer término y $F_0=S_0e^{(r-q)T}$.)

Esta formulación demuestra que, **si el precio tiene trayectorias continuas**, la varianza se puede valorar y replicar estáticamente comprando un *strip* infinito de opciones OTM (fuera del dinero) de *puts* y *calls*, ponderadas exactamente por el inverso del cuadrado de su *strike* ($\frac{1}{K^2}$). La sección siguiente muestra qué se rompe cuando esa hipótesis falla.

---

## Límite de la réplica: saltos

### ¿Qué varianza paga realmente el contrato?

Hay tres objetos que bajo difusión pura coinciden, pero con saltos no:

1. la variación cuadrática continua del log-precio, $[X]_T/T$ con $X_t=\ln S_t$;
2. la suma discreta de cuadrados de log-retornos diarios (lo que el contrato paga, con $A=252$ y sin restar la media);
3. la cantidad que replica la tira de opciones, $\frac{2}{T}\left(\int_0^T \frac{dS_t}{S_{t^-}}-\ln\frac{S_T}{S_0}\right)$.

La suma discreta (2) converge a $[X]_T/T$ (1) al refinar el muestreo, y su diferencia es típicamente pequeña (Broadie y Jain, 2008). El punto crítico es la diferencia entre (1) y (3).

### Dinámica con saltos e identidad trayectoria por trayectoria

Se extiende el modelo con saltos de actividad finita. Bajo $\mathbb{Q}$,

$$
\frac{dS_t}{S_{t^-}}=(r-\lambda\kappa)\, dt+\sigma_t\, dW_t+\left(e^{J}-1\right)dN_t,\qquad \kappa=\mathbb{E}[e^{J}-1],
$$

donde $N_t$ es un proceso de Poisson de intensidad $\lambda$ y $J$ es el salto del log-precio. El término $-\lambda\kappa\, dt$ compensa los saltos, de modo que $\mathbb{E}^{\mathbb{Q}}\!\left[\int_0^T \frac{dS_t}{S_{t^-}}\right]=rT$ sigue valiendo.

Aplicando Itô con saltos, en cada instante de salto el log-precio cambia en $\, delta X=J$ mientras que $S$ cambia en $e^{\, delta X}-1$ en términos relativos. Por tanto

$$
\ln\frac{S_T}{S_0}=\int_0^T\frac{dS_t}{S_{t^-}}-\frac12\int_0^T\sigma_t^2\, dt-\sum_{s\le T}\left(e^{\, delta X_s}-1-\, delta X_s\right).
$$

La cantidad que replica la tira de opciones es $R_T:=2\left(\int_0^T \frac{dS_t}{S_{t^-}}-\ln\frac{S_T}{S_0}\right)$, y la variación cuadrática es $[X]_T=\int_0^T\sigma_t^2dt+\sum_{s\le T}(\, delta X_s)^2$. Restando:

$$
R_T-[X]_T=2\sum_{s\le T} g(\, delta X_s),\qquad g(x):=e^{x}-1-x-\tfrac12x^2 .
$$

La función $g$ tiene $g(0)=g'(0)=g''(0)=0$ y $g'''(x)=e^{x}>0$, de donde:

* $g(x)=\tfrac16x^3+\tfrac1{24}x^4+\cdots\approx\tfrac16x^3$ para saltos moderados (el error es de **tercer orden**);
* $g(x)$ tiene siempre el **mismo signo que $x$**: un salto bajista siempre resta, uno alcista siempre suma.

### El sesgo de replicación

Tomando esperanza bajo $\mathbb{Q}$ (con $\mathbb{E}\sum g(\, delta X_s)=\lambda T\,\mathbb{E}[g(J)]$):

$$
\boxed{K^{\text{repl}}-K^{QV}=2\lambda\,\mathbb{E}\!\left[e^{J}-1-J-\tfrac12J^2\right]\approx\frac{\lambda}{3}\,\mathbb{E}[J^3]}
$$

donde $K^{\text{repl}}=\mathbb{E}^{\mathbb{Q}}[R_T]/T$ es el valor que entrega la tira de opciones en las ecuaciones (3)-(4) (depende solo de precios de opciones observados) y $K^{QV}=\mathbb{E}^{\mathbb{Q}}\big[[X]_T\big]/T$ es la esperanza de la varianza que el contrato paga en el límite de muestreo fino.

* Si los saltos son negativos en media, $\mathbb{E}[J^3]<0$ y $K^{\text{repl}}<K^{QV}$: el *strike* *model-free* de las ecuaciones (3)-(4) **subestima** la varianza cuadrática esperada.
* Lo que el swap paga es $\approx K^{QV}$, **no** $K^{\text{repl}}$. Al comparar modelo contra mercado hay que reportar ambas cantidades y no mezclar definiciones.
* Para quien vende el swap y se cubre con la réplica, el descalce es $R_T-[X]_T$: pierde con saltos bajistas y gana con alcistas, con un efecto cúbico en el tamaño del salto. Es consistente con lo que señalan Carr y Lee (2009) sobre la exposición de estos contratos a retornos al cubo y de orden superior en movimientos bruscos.
* El resultado no es una curiosidad teórica: Broadie y Jain (2008) encuentran que el efecto del muestreo discreto suele ser pequeño mientras que el de los saltos puede ser significativo. Carr, Lee y Lorig (2021) desarrollan réplicas robustas para difusiones con saltos de actividad finita y tamaño acotado.
---

## Varianza forward y estructura temporal

La varianza total esperada no se distribuye de manera uniforme a lo largo del tiempo, sino que es estrictamente aditiva. Con el fin de aislar la expectativa de volatilidad que ocurrirá de forma exclusiva entre dos vencimientos futuros $T_1$ y $T_2$, se emplea la estructura temporal de la varianza *forward* (Demeterfi et al., 1999). Esta relación matemática se define como:

$$
K_{T_1,T_2}=\frac{T_2K_{T_2}-T_1K_{T_1}}{T_2-T_1}.
$$

La ecuación anterior se construye ponderando la varianza total de cada vencimiento respecto a su respectivo horizonte temporal y dividiendo el resultado entre la diferencia de ambos plazos. El principal uso de esta estrategia se emplea mediante la adquisición de un *variance swap* con vencimiento en $T_2$ y la venta simultánea de otro con vencimiento en $T_1$, permitiendo asi identificar en qué segmento del calendario se concentra la incertidumbre.

---

## Saltos genéricos vs. saltos de earnings

La sección de saltos supone que estos llegan en instantes aleatorios (Poisson). Un anuncio de resultados (*earnings announcement*, EA) es distinto: la **fecha $\tau$ se conoce con antelación**, aunque no se conozca la respuesta del precio. Esto cambia qué se puede identificar y cómo se lee el sesgo.

### Descomposición unificada

Para comprender cómo se compone la volatilidad total de un activo, la descomposición unificada integra la varianza cuadrática esperada ($QV$) separándola en tres fuentes fundamentales de aleatoriedad del precio:

$$_{T}K_T^{QV} = \int_0^T \mathbb{E}^\mathbb{Q}[\sigma_t^2] dt + \lambda_c T \mathbb{E}[J_c^2] + \mathbb{1}_{\tau \le T} \mathbb{E}[J_{EA}^2]$$

Esta formulación se fundamenta en los modelos clásicos de difusión con saltos (Merton, 1976; Cont & Tankov, 2004). El primer componente de esta expresión corresponde a la difusión continua, modelada mediante la integral esperada de la varianza en la trayectoria del activo (El movimiento normal diario del precio). El segundo componente agrupa los saltos aleatorios, siguiendo un proceso de Poisson. Finalmente, el tercer componente incorpora los saltos programados asociados a los anuncios de resultados (*earnings*), cuya fecha de ocurrencia $\tau$ si se conoce, diferenciándose así de los choques estocásticos del mercado.

---

### Replicación con opciones y el sesgo

No obstante, cuando se negocian opciones en los mercados financieros, la cadena no mide la varianza cuadrática de forma directa, sino que incorpora una prima adicional debido a la asimetría de mercado, las colas gordas y el skew por miedo a desplomes en los precios. Esta discrepancia se formaliza mediante la ecuación de replicación (Neuberger, 1994; Carr & Madan, 1998), la cual ajusta la varianza realizada sumando los términos representados por la función especial de sesgo para cada tipo de salto:

$$_{T}K_T^{\text{repl}} = {}_{T}K_T^{QV} + 2\lambda_c T \mathbb{E}[g(J_c)] + 2\mathbb{1}_{\tau \le T} \mathbb{E}[g(J_{EA})]$$

Dicha función ($g(x) = e^x - 1 - x - \frac{1}{2}x^2$) mide exactamente la distorsión introducida por el contrato logarítmico frente a los movimientos extremos del subyacente. En ese caso, $e^x$ se dispara en valor, incrementando la pendiente de la funcion.

---

### 4. Aislamiento del salto de *earnings*

Para aislar el impacto exclusivo de un reporte de resultados (*earnings*), se resta la varianza calculada por opciones ($K^{\text{repl}}$) de la varianza realizada pura ($QV$), lo que hace que el ruido normal del mercado y los saltos aleatorios se cancelen. Al aplicar una expansión matemática (el desarrollo en serie de Taylor) sobre la función $g(x)$, los términos iniciales se eliminan y sobrevive el término cúbico. Al combinarlo con los factores de la fórmula, surge directamente un tercio del momento cúbico ($\mathbb{E}[J_{EA}^3]$)

$$K_T^{\text{repl}} - K_T^{QV}|_{EA} = \frac{2}{T} \mathbb{E}[g(J_{EA})] \approx \frac{1}{3T} \mathbb{E}[J_{EA}^3]$$

Esto demuestra matemáticamente que el sesgo de los *earnings* refleja la asimetría o el riesgo direccional extremo que el mercado anticipa ante el anuncio.

### Comparación conceptual

| | Salto genérico (Poisson) | Salto de earnings (programado) |
|---|---|---|
| Momento | Aleatorio, intensidad $\lambda_c$ | Fecha $\tau$ conocida |
| Tamaño | Aleatorio; típicamente sesgo negativo (riesgo de caída) | Aleatorio y grande; signo *a priori* incierto |
| Aporte a $K^{QV}$ | $\lambda_c\,\mathbb{E}[J_c^2]$, repartido de forma uniforme en el horizonte | $\mathbb{E}[J_{EA}^2]/T$, concentrado en una fecha |
| Sesgo de la réplica | $2\lambda_c\,\mathbb{E}[g(J_c)]$ | $\tfrac{2}{T}\,\mathbb{E}[g(J_{EA})]$ |
| Cómo se identifica con opciones | Solo vía *smile* (asimetría y curtosis), mezclado con la difusión | Vía estructura temporal: escalón de varianza entre vencimientos que rodean a $\tau$ |
| Efecto observable | Suaviza la superficie | La volatilidad implícita sube antes del anuncio y cae al conocerse; estructura temporal decreciente a corto plazo |


# 2. El modelo de Merton (*jump-diffusion*)

La replicación estática que desarrollamos antes parte de un supuesto clave que usamos: el paso del Lema de Itô para transformar $\ln S_T$ en una integral de $\sigma_t^2$ exige que las trayectorias del precio sean **continuas**. Sin embargo, un evento como un anuncio de resultados rompe esta regla por completo. El precio no evoluciona suavemente durante la noche en que se publican los resultados; por el contrario, *salta* de golpe de un nivel a otro entre el cierre y la apertura del mercado.

Para capturar este comportamiento, Merton (1976) propuso un modelo que mantiene la difusión browniana para modelar el ruido de alta frecuencia del día a día y le añade un **proceso de saltos de Poisson**. En nuestra investigación, este modelo cumple tres propósitos esenciales:

1. Genera de forma natural **colas pesadas y el característico *smile*** en vencimientos cortos, que es justo cuando un movimiento drástico domina la distribución de $S_T$.
2. Funciona como un **laboratorio analítico** muy limpio, permitiéndonos ver con exactitud matemática qué partes de la fórmula (4) dejan de cumplirse cuando hay saltos.
3. El término de los saltos se interpreta de inmediato como la **varianza de evento** asociada al anuncio de resultados.

---

## Estructura del modelo

### El proceso de Poisson y los saltos

Sea $N_t$ un proceso de Poisson con intensidad $\lambda$, que representa el número esperado de saltos por año. Sus propiedades básicas se conocen:

$$
\mathbb{P}(N_t=n)=\frac{e^{-\lambda t}(\lambda t)^n}{n!},\qquad \mathbb{E}[N_t]=\text{Var}[N_t]=\lambda t
$$

con incrementos independientes a lo largo del tiempo. En un intervalo infinitesimal $dt$, la probabilidad de que ocurra más de un salto es muy baja, de modo que hay a lo sumo un salto:

$$
\mathbb{P}(dN_t = 1) = \lambda dt, \qquad \mathbb{P}(dN_t = 0) = 1 - \lambda dt
$$

> **Nota explicativa sobre el cálculo estocástico con saltos:**
> Del mismo modo que en el movimiento browniano estándar se cumple $(dW_t)^2=dt$, en un proceso de saltos el incremento $dN_t$ solo puede valer $0$ o $1$, lo que implica que $(dN_t)^2=dN_t$. Por otro lado, los términos cruzados se anulan: $dN_t\, dt=0$ y $dW_t\, dN_t=0$, ya que los saltos y el recorrido de difusión del precio son independientes y no ocurren al mismo tiempo. Como consecuencia, la variación cuadrática de los saltos no se acumula de forma continua con el tiempo; lo hace **únicamente** en el momento exacto del salto, y su valor es simplemente el cuadrado del tamaño de dicho salto.

Cada vez que un salto de ese estilo ocurre, el precio se multiplica por un factor estocástico $e^{J}$, donde el **tamaño logarítmico del salto** sigue una distribución normal:

$$
J_i\sim\mathcal{N}(\mu_J,\sigma_J^2)\ \text{iid, independiente de }W\text{ y de }N .
$$

Por ende, el cambio porcentual que experimenta el precio en el salto es $e^{J}-1$. Su valor esperado se calcula con la función generadora de momentos de la normal, sabiendo que $\mathbb{E}[e^{J}]=e^{\mu_J+\frac12\sigma_J^2}$:

$$
k:=\mathbb{E}[e^{J}-1]=e^{\mu_J+\frac12\sigma_J^2}-1 .
$$

Finalmente, si sumamos $N_t$ saltos consecutivos (lo que se conoce como un **proceso de Poisson compuesto**), aplicando la ley de la varianza total obtenemos:

$$
\mathbb{E}\Big[\sum_{i=1}^{N_t}J_i\Big]=\lambda t\,\mu_J,\qquad
\text{Var}\Big[\sum_{i=1}^{N_t}J_i\Big]=\underbrace{\mathbb{E}[N_t]\text{Var}[J]}_{\lambda t\,\sigma_J^2}+\underbrace{\text{Var}[N_t]\,\mathbb{E}[J]^2}_{\lambda t\,\mu_J^2}=\lambda t\,(\mu_J^2+\sigma_J^2)
$$

El primer termino es la incertidumbre debida al tamaño de los saltos individuales multiplicada por el número esperado de saltos. El segundo termino es la incertidumbre debida a que el número mismo de saltos también varía con el tiempo. Si sumamos ambos terminos y factorizamos $$\lambda t\$$, obtenemos $$\lambda t\(\mu_J^2+\sigma_J^2)$$.


Aquí, el término $\mu_J^2+\sigma_J^2=\mathbb{E}[J^2]$ aparecerá más adelante de forma natural como la contribución directa de los saltos a la varianza realizada.

### La dinámica del precio y la compensación del salto

Bajo la medida histórica $\mathbb{P}$, el precio del activo evoluciona según la expresión:

$$
\frac{dS_t}{S_{t^-}} = (\mu - \lambda k) dt + \sigma dW_t + (e^J - 1) dN_t
$$

donde $S_{t^-}$ representa el nivel del precio justo antes de que se produzca un salto. El término $-\lambda k\, dt$ cumple un papel fundamental: actúa como el **compensador**. Dado que los saltos por sí mismos generan un valor esperado de crecimiento ($\lambda k\, dt$), restarlos dentro de la tendencia (*drift*) funciona como un freno que evita que el retorno promedio esperado del browniano suba en exceso y se desvíe, manteniéndose exactamente en su tendencia base $\mu$, es decir, que $\mathbb{E}[dS_t/S_{t^-}]=\mu\, dt$.

Al cambiar a la medida neutral al riesgo $\mathbb{Q}$, el *drift* total del proceso pasa a ser la tasa libre de riesgo $r$ (si hubiese un rendimiento por dividendos continuo $q$, simplemente sustituiríamos $r$ por $r-q$):

$$
\boxed{\frac{dS_t}{S_{t^-}} = (r - \lambda k) dt + \sigma dW_t + (e^J - 1) dN_t} \tag{5}
$$

---

## El log-precio: Itô con saltos

De forma idéntica al caso puramente continuo, nos interesa estudiar la función $f(S)=\ln S$. La gran diferencia ahora es que $S_t$ combina dos fuentes de movimiento distintas, por lo que el Lema de Itô se aplica analizando ambos componentes por separado:

- **En los periodos entre saltos**, el precio se comporta exactamente igual que en un movimiento browniano geométrico con un *drift* de $(r-\lambda k)$. Aplicando la misma derivada parcial que antes, con $f'(S)=1/S$ y $f''(S)=-1/S^2$, surge de manera natural el ajuste por convexidad $-\frac12\sigma^2dt$.
- **En el momento exacto de un salto**, el precio salta bruscamente de $S_{t^-}$ a $S_{t^-}e^{J}$. Al aplicar el logaritmo, el cambio que experimenta el precio es $$\ln\big(S_{t^-}e^{J}\big)-\ln S_{t^-}=J ,$$ y **no** $e^{J}-1$ (que era el cambio porcentual del precio al momento del salto).

Juntando ambos efectos en una sola expresión, tenemos:

$$
\boxed{d(\ln S_t)=\Big(r-\lambda k-\tfrac12\sigma^2\Big)dt+\sigma\, dW_t+J\, dN_t} \tag{6}
$$

> Vale la pena notar una distinción crucial: mientras que el retorno simple $dS/S$ experimenta un salto de magnitud $e^{J}-1$, el retorno logaritmico $d\ln S$ salta exactamente en $J$. En un día común del mercado sin saltos abruptos, ambos coinciden a primer orden y la única diferencia es el término de convexidad $\frac12\sigma^2dt$. Sin embargo, al analizar los activos en el mercado real aparece una discrepancia **adicional** de $e^{J}-1-J$ en cada evento de salto. Esta diferencia (que suele ser positiva ya que la magnitud de $e^{J}-1$ es mayor que la magnitud $J$) es la raíz del sesgo de replicación que demostraremos más adelante y explica el por que es casi imposible replicar el variance swap con opciones que tienen retornos logaritmicos (Carr-Madan).

Si integramos la ecuación (6) a lo largo del intervalo temporal $[0,T]$, obtenemos:

$$
\ln\frac{S_T}{S_0}=\Big(r-\lambda k-\tfrac12\sigma^2\Big)T+\sigma W_T+\sum_{i=1}^{N_T}J_i \tag{7}
$$

### Distribución: una mezcla de Poisson de normales

Si supiéramos el numero exacto de saltos que ocurrirán $N_T=n$, tanto el segundo como el tercer sumando de (7) se comportarian como variables con un comportamiento normal. Esto implica que el log-retorno $\ln(S_T/S_0)$ sigue una distribución normal condicionada con:

$$
m_n=\Big(r-\lambda k-\tfrac12\sigma^2\Big)T+n\,\mu_J,\qquad v_n=\sigma^2T+n\,\sigma_J^2 ,
$$

Al no conocer el numero fijo de saltos, la distribución real resulta ser una mezcla ponderada de normales donde cada campana está pesada por la probabilidad de que ocurran n saltos: $\sum_{n}\frac{e^{-\lambda T}(\lambda T)^n}{n!}\,\mathcal{N}(m_n,v_n)$.

### Momentos y cumulantes

A partir de la media y la varianza de cada componente, definimos $R=\ln(S_T/S_0)$ y deducimos fácilmente sus propiedades globales (Ambas analizadas a lo largo del tiempo $T$:

$$
\mathbb{E}[R]=\Big(r-\tfrac12\sigma^2-\lambda k+\lambda\mu_J\Big)T,\qquad
\text{Var}[R]=\Big(\sigma^2+\lambda(\mu_J^2+\sigma_J^2)\Big)T .
$$

Para analizar órdenes superiores, la teoría nos muestra que los cumulantes de un proceso de Poisson compuesto se calculan como $\kappa_n = \lambda T \mathbb{E}[J^n]$. Estos valores se obtienen expandiendo en serie el logaritmo de la función característica, donde la difusión tradicional solo aporta a los dos primeros momentos. Sabiendo que $\mathbb{E}[J^3] = \mu_J^3 + 3\mu_J \sigma_J^2$ y $\mathbb{E}[J^4] = \mu_J^4 + 6\mu_J^2 \sigma_J^2 + 3\sigma_J^4$, podemos calcular:

$$
\text{Asimetría}=\frac{\lambda T\,\mathbb{E}[J^3]}{\text{Var}[R]^{3/2}},\qquad
\text{Curtosis en exceso}=\frac{\lambda T\,\mathbb{E}[J^4]}{\text{Var}[R]^{2}} .
$$

> El signo de la asimetría (*skewness*) depende del parámetro $\mu_J$ (por ejemplo, un "crash" bajista tras un reporte de resultados genera una asimetría negativa marcada), mientras que la curtosis disminuye inversamente proporcional a $1/(\lambda T)$. En una ventana corta de 5 a 10 días alrededor de los anuncios, el número esperado de saltos $\lambda T$ ronda el valor de 1, lo que provoca que las colas de la distribución aumenten en una gran proporsion. Esto explica con claridad por qué la dinámica de saltos domina la forma del *smile* de volatilidad en vencimientos de corto plazo.

---

## Función característica

Para implementar modelos de valuación basados en métodos de Fourier, necesitamos calcular la función característica $\varphi(u)=\mathbb{E}^{\mathbb{Q}}[e^{iu\ln S_T}]$. Esto se logra condicionando sobre $N_T$ y aprovechando que la difusión y los saltos son independientes entre sí. La parte correspondiente a la difusión aporta un factor $\exp\!\big(iu\,m\,T-\tfrac12\sigma^2u^2T\big)$, con $m=r-\lambda k-\tfrac12\sigma^2$. Por otro lado, utilizando la función generadora de los saltos $\psi(u)=\mathbb{E}[e^{iuJ}]=e^{iu\mu_J-\frac12\sigma_J^2u^2}$, resolvemos la esperanza del proceso de Poisson:

$$
\mathbb{E}\Big[e^{iu\sum_{i=1}^{N_T}J_i}\Big]=\sum_{n=0}^{\infty}\frac{e^{-\lambda T}(\lambda T)^n}{n!}\,\psi(u)^n=e^{-\lambda T}e^{\lambda T\psi(u)}=\exp\Big(\lambda T\big(\psi(u)-1\big)\Big)
$$

Al multiplicar ambas contribuciones, llegamos a la expresión compacta:

$$
\boxed{\varphi_M(u)=\exp \Big\{ iu\big[\ln S_0+(r-\lambda k-\tfrac{1}{2}\sigma^2)T\big] - \tfrac{1}{2}\sigma^2u^2T + \lambda T \Big(e^{iu\mu_J - \frac{1}{2}\sigma_J^2u^2} - 1\Big) \Big\}} \tag{8}
$$

> **Relevancia para la investigación:**
> La fórmula (8) es la pieza central para aplicar técnicas de inversión de Fourier (como los algoritmos de Carr–Madan o Gil-Pelaez), mediante las cuales calculamos con gran eficiencia los precios teóricos de opciones *call* y *put* bajo el modelo de Merton y, a partir de ellos, evaluamos numéricamente la integral de replicación (4). Además, verificar que $\varphi_M(-i)=F_0$ nos sirve como una robusta prueba de consistencia numérica en cada rutina de calibración.

---

## Precio de opciones europeas: la serie de Merton

Si condicionamos el sistema a que ocurran exactamente $n$ saltos ($N_T=n$), el log-precio resultante adopta una distribución perfectamente normal, lo que significa que cada término puede resolverse mediante una fórmula clásica de Black–Scholes, ajustando adecuadamente sus parámetros equivalentes $(\sigma_n,r_n)$:

- **Varianza total ajustada:** $\sigma_n^2T=v_n=\sigma^2T+n\sigma_J^2$, lo que se traduce en una volatilidad equivalente $\sigma_n^2=\sigma^2+\, dfrac{n\sigma_J^2}{T}$.
- **Tasa *forward* ajustada:** $\mathbb{E}[S_T \mid N_T=n]=S_0\exp\big((r-\lambda k)T+n(\mu_J+\tfrac{1}{2}\sigma_J^2)\big)=S_0e^{r_nT}$, lo que define una tasa de interés modificada $r_n=r-\lambda k+\, dfrac{n\ln(1+k)}{T}$.
  
Mientras que el modelo de Black–Scholes estándar descuenta los flujos usando la tasa $r_n$ a través de $e^{-r_nT}$, el descuento correcto en este marco exige utilizar la tasa libre de riesgo $r$ mediante $e^{-rT}$. El factor de corrección necesario es $e^{(r_n-r)T}=e^{-\lambda kT}(1+k)^n$, el cual, al combinarse con las probabilidades del peso de Poisson, nos deja con:

$$
e^{-\lambda T}\frac{(\lambda T)^n}{n!}\,e^{-\lambda kT}(1+k)^n=\frac{e^{-\lambda' T}(\lambda' T)^n}{n!},\qquad \lambda'=\lambda(1+k).
$$

De este desarrollo surge la elegante serie analítica de Merton para opciones europeas:

$$
\boxed{C_M(K)=\sum_{n=0}^{\infty}\frac{e^{-\lambda' T}(\lambda' T)^n}{n!}\;C_{BS}\big(S_0,K,T,\sigma_n,r_n\big)} \tag{9}
$$

Los valores para las opciones *put* se derivan directamente aplicando la paridad *put-call* (la cual se cumple sin problemas gracias a que la condición $\mathbb{E}[S_T]=F_0$ se preserva).


> La ecuación (9) representa una **mezcla ponderada de distribuciones de Poisson sobre modelos Black–Scholes**: cada escenario con $n$ saltos cuenta con su propia volatilidad efectiva $\sigma_n$, y la combinación de estas distintas volatilidades es lo que genera de forma natural el *smile* observado. Dado que en las ventanas de resultados corporativos el parámetro $\lambda' T\approx1$, los pesos de la serie decrecen rápidamente de forma factorial, por lo que bastan entre 15 y 20 términos para lograr una convergencia perfecta. Esta serie funciona como un método de validación cruzada independiente frente al enfoque de Fourier: ambos procedimientos deben arrojar exactamente los mismos precios.

---

## El variance swap bajo Merton: qué se rompe con saltos

Una vez establecida la arquitectura matemática del modelo de saltos, podemos analizar con precisión quirúrgica qué partes de los resultados anteriores siguen vigentes. Para ello, compararemos el pago real que efectúa el contrato frente a la cantidad que logra replicar la cartera estática de opciones.

### Lo que paga el contrato: la variación cuadrática

Por definición, la varianza realizada que liquida el contrato al vencimiento corresponde a la suma de los cuadrados de los log-retornos observados. Al refinar la partición temporal de la muestra, esta sumatoria converge analíticamente hacia la **variación cuadrática** del proceso del log-precio, denotada por $[X]_T$, la cual absorbe de forma natural el impacto de los saltos. Aplicando las reglas de cálculo estocástico $(dW_t)^2=dt$, $(dN_t)^2=dN_t$ y $dW_t\, dN_t=0$ sobre el diferencial (6):

$$
(d\ln S_t)^2=\sigma^2\, dt+J^2\, dN_t\ \ \Longrightarrow\ \ [X]_T=\sigma^2T+\sum_{i=1}^{N_T}J_i^2 . \tag{10}
$$

Si tomamos la esperanza matemática de esta expresión, recordando que $\mathbb{E}[N_T]=\lambda T$ y $\mathbb{E}[J^2]=\mu_J^2+\sigma_J^2$, obtenemos:

$$
\boxed{K^{QV}=\frac{1}{T}\mathbb{E}^{\mathbb{Q}}\big[[X]_T\big]=\sigma^2+\lambda\,(\mu_J^2+\sigma_J^2)} \tag{11}
$$

> La diferencia $K^{QV}-\sigma^2=\lambda(\mu_J^2+\sigma_J^2)$, multiplicada por el horizonte $T$, representa exactamente la **varianza de evento** asociada al anuncio de earnings: $\varepsilon=\mathbb{E}\big[\sum_{i=1}^{N_T}J_i^2\big]=\lambda T(\mu_J^2+\sigma_J^2)$.

### La tira de opciones

La fórmula de replicación se fundamento en dos pilares: el teorema de Carr–Madan y la esperanza neutral al riesgo $\mathbb{E}^{\mathbb{Q}}[\int dS_t/S_t]=rT$ (que también se mantiene intacta gracias a la acción del compensador $-\lambda k$). Por consiguiente, la cartera de opciones continúa replicando exactamente la misma cantidad teórica:

$$
K^{\text{repl}}=\frac{2}{T}\,\mathbb{E}^{\mathbb{Q}}\Big[\int_0^T\frac{dS_t}{S_{t^-}}-\ln\frac{S_T}{S_0}\Big].
$$

Lo único que deja de cumplirse es la equivalencia del Lema de Itô que antes permitía simplificar esa integral directamente en $\int\sigma_t^2dt$. Para comprender qué representa ahora esta cantidad, restamos los diferenciales del retorno simple $dS/S$ en (5) y del log-retorno $d\ln S$ en (6):

$$
\frac{dS_t}{S_{t^-}}-d\ln S_t=\tfrac12\sigma^2\, dt+\big(e^{J}-1-J\big)\, dN_t .
$$

Integrando a lo largo del periodo y multiplicando por el factor 2, llegamos a:

$$
2\Big[\int_0^T\frac{dS_t}{S_{t^-}}-\ln\frac{S_T}{S_0}\Big]=\sigma^2T+2\sum_{i=1}^{N_T}\big(e^{J_i}-1-J_i\big). \tag{12}
$$

Si contrastamos la ecuación (12) con la de la variación cuadrática (10), notamos que en un régimen puramente continuo ambos lados coinciden perfectamente; sin embargo, en presencia de saltos, cada evento aporta una contribución diferente: $2(e^{J}-1-J)$ para el lado replicable frente a $J^2$ para la variación cuadrática. Al tomar la esperanza matemática en (12), sabiendo que $\mathbb{E}[e^{J}-1-J]=k-\mu_J$, deducimos:

$$
\boxed{K^{\text{repl}}=\sigma^2+2\lambda\,(k-\mu_J)=\sigma^2+2\lambda\big(e^{\mu_J+\frac12\sigma_J^2}-1-\mu_J\big)} \tag{13}
$$

# 3. Modelo de Heston y modelo de Bates (SV + saltos)

El modelo de Merton resolvió el problema de los saltos, pero dejó un supuesto incómodo: la volatilidad de difusión $\sigma$ es **constante**. Como el variance swap paga justamente la varianza realizada, y esta cambia con el tiempo (hay periodos calmados y periodos turbulentos), conviene que la varianza también sea aleatoria. Heston (1993) propone esto, y Bates (1996) lo combina con los saltos de Merton. Para la investigación sirve para tres cosas:

1. Darle **dinámica a la varianza**, lo que permite hablar de una estructura temporal y no solo de un nivel.
2. **Separar el efecto de la volatilidad estocástica del de los saltos**, comparando directamente con Merton.
3. Cambiar la lectura del *strike*: la parte continua deja de ser $\sigma^2$ y pasa a ser la varianza esperada promedio.

---

## Estructura del modelo

### La varianza como proceso aleatorio (Heston)

Heston modela la varianza instantánea $v_t=\sigma_t^2$ con el siguiente proceso:

$$
dv_t=\kappa(\theta-v_t)\,dt+\xi\sqrt{v_t}\,dW_t^{v} \tag{14}
$$

* $\theta$: nivel de largo plazo al que la varianza tiende a volver.
* $\kappa$: velocidad con la que vuelve a $\theta$ (reversión a la media).
* $\xi$: volatilidad de la volatilidad, es decir, cuánto fluctúa la varianza al azar.
* $v_0$: varianza de hoy.

La raíz cuadrada hace que la varianza se mueva más cuando ya es alta (los periodos turbulentos se agrupan), esto impide que $v_t$ llegue a ser negativa.

El ruido de la varianza está correlacionado con el del precio: $d\langle W^S,W^v\rangle_t=\rho\,dt$. Normalmente $\rho<0$: cuando el precio cae, la varianza sube (efecto apalancamiento), y eso produce la asimetría negativa del *smile* en la parte de difusión.

> **Nota explicativa sobre el cálculo estocástico:**
> Con dos brownianos correlacionados se cumple $(dW_t^S)^2=(dW_t^v)^2=dt$, pero el término cruzado ya no se anula: $dW_t^S\,dW_t^v=\rho\,dt$. Los saltos siguen siendo independientes de ambos: $dN_t\,dW_t^S=dN_t\,dW_t^v=0$ y $(dN_t)^2=dN_t$, como en Merton.

**¿Cuánta varianza esperamos en promedio?** Esto es lo que necesita el variance swap, porque el contrato paga el promedio de $v_t$. Al tomar esperanza en (14), la parte aleatoria desaparece (la esperanza de una integral de Itô es cero) y queda una ecuación simple para $\mathbb{E}[v_t]$, cuya solución es

$$
\mathbb{E}[v_t]=\theta+(v_0-\theta)\,e^{-\kappa t},\qquad
\bar v(T):=\frac1T\,\mathbb{E}\!\int_0^T v_t\,dt=\theta+\frac{(v_0-\theta)\big(1-e^{-\kappa T}\big)}{\kappa T} . \tag{15}
$$

En esta ecuacion, $\bar v(T)$ la varianza esperada es el nivel de largo plazo $v_0$ mas la distancia actual entre la varianza promedio y el valor de largo plazo $(v_0-\theta)$. Como los contratos de opciones o *Variance Swaps* duran un periodo determinado $T$ (por ejemplo, 30 días, 3 meses o 1 año), lo que realmente nos interesa es calcular el promedio de la varianza esperada desde $t = 0$ hasta $t = T$:

En este caso, * **A plazos muy cortos ($T \to 0$):** El término del exponente decae poco, por lo que la varianza promedio $\bar{v}(T)$ estará sumamente influenciada por el nivel actual del mercado ($v_0$).

* **A plazos muy largos ($T \to \infty$):** El factor $(1 - e^{-\kappa T}) / (\kappa T)$ tiende a cero, haciendo que la varianza promedio converja hacia el nivel de largo plazo ($\theta$), ya que la velocidad de reversión a la media ($\kappa$) habrá absorbido todos los movimientos durante este periodo.




Es el equivalente de $\sigma^2$ en Merton.

### Bates: Heston más los saltos de Merton

Bates simplemente le agrega a Heston el mismo proceso de saltos de la sección anterior. Bajo $\mathbb{Q}$ (con dividendo continuo $q$ se cambia $r$ por $r-q$):

$$
\boxed{\frac{dS_t}{S_{t^-}}=(r-\lambda k)\,dt+\sqrt{v_t}\,dW_t^{S}+(e^{J}-1)\,dN_t} \tag{16}
$$

$$
\boxed{dv_t=\kappa(\theta-v_t)\,dt+\xi\sqrt{v_t}\,dW_t^{v},\qquad d\langle W^S,W^v\rangle_t=\rho\,dt} \tag{17}
$$

donde $N_t$, $J\sim\mathcal{N}(\mu_J,\sigma_J^2)$ y $k=e^{\mu_J+\frac12\sigma_J^2}-1$ son los mismos de Merton. El modelo tiene **8 parámetros**: cinco de Heston $(v_0,\kappa,\theta,\xi,\rho)$ y tres de saltos $(\lambda,\mu_J,\sigma_J)$. Si apagamos ingredientes, recuperamos los modelos anteriores:

| Si... | El modelo se reduce a... |
|-------|--------------------------|
| $\lambda=0$ | Heston |
| $\xi=0$ y $v_0=\theta$ (varianza constante) | Merton, con $\sigma^2=\theta$ |
| $\lambda=0$, $\xi=0$ y $v_0=\theta$ | Black–Scholes |

Esto es la base del experimento anidado: Black–Scholes $\to$ Merton $\to$ Heston $\to$ Bates, agregando **un solo ingrediente a la vez**, para poder atribuir cada cambio en el *strike* a una causa concreta. Como en Merton, los parámetros calibrados son bajo $\mathbb{Q}$ (incluyen primas de riesgo), lo que los hace comparables con el mercado.

---

## El log-precio: Itô con volatilidad estocástica y saltos

El procedimiento es el mismo de Merton, cambiando $\sigma^2$ por $v_t$: entre saltos aparece el ajuste de convexidad $-\frac12v_t\,dt$, y en un salto el log-precio cambia en $J$.

$$
\boxed{d(\ln S_t)=\Big(r-\lambda k-\tfrac12 v_t\Big)dt+\sqrt{v_t}\,dW_t^{S}+J\,dN_t} \tag{18}
$$

Lo importante es que al integrar aparece $\int_0^T v_t\,dt$, que es justamente la **varianza realizada** de la primera sección. En Merton esa integral valía $\sigma^2T$ (fija); ahora es aleatoria, y lo que nos interesa de ella es su esperanza.

---

## Función característica

Igual que en Merton, para valorar opciones con Fourier necesitamos $\varphi_B(u)=\mathbb{E}^{\mathbb{Q}}[e^{iu\ln S_T}]$. Como la parte de Heston y los saltos son independientes, la función característica es el **producto** de ambas. La parte de saltos ya la tenemos de (8); la de Heston sale así:

1. La esperanza condicional $\mathbb{E}[e^{iu\ln S_T}\mid S_t,v_t]$ es una martingala, así que al aplicarle Itô su *drift* debe ser cero. Eso da una ecuación en derivadas parciales.
2. Heston es un modelo **afín**, entonces se propone una solución de la forma $\exp\{iu\ln S_t+C(\tau)+D(\tau)v_t\}$, con $\tau=T-t$. Al sustituirla y separar los términos que contienen $v$ de los que no, quedan dos ecuaciones diferenciales ordinarias (tipo Riccati):

$$
\partial_\tau D=\tfrac12\xi^2D^2+(\rho\xi iu-\kappa)D-\tfrac12(u^2+iu),\qquad
\partial_\tau C=\kappa\theta D+iu(r-q),\qquad D(0)=C(0)=0 . \tag{19}
$$

3. La primera es cuadrática en $D$ y se resuelve de forma explícita. Con $b=\kappa-\rho\xi\,iu$, $d=\sqrt{b^2+\xi^2(iu+u^2)}$ y $g=\frac{b-d}{b+d}$:

$$
D=\frac{b-d}{\xi^2}\cdot\frac{1-e^{-dT}}{1-g\,e^{-dT}},\qquad
C=\frac{\kappa\theta}{\xi^2}\Big[(b-d)\,T-2\ln\frac{1-g\,e^{-dT}}{1-g}\Big]+iu(r-q)T . \tag{20}
$$

(Esta es la versión estable de Albrecher et al. (2007), *Little Heston Trap*, que evita problemas numéricos con el logaritmo complejo.) Multiplicando por la parte de saltos, incluido el compensador $-iu\lambda kT$:

$$
\boxed{\varphi_B(u)=\exp\Big\{iu\ln S_0+C(u,T)+D(u,T)\,v_0+\lambda T\Big(e^{iu\mu_J-\frac12\sigma_J^2u^2}-1-iuk\Big)\Big\}} \tag{21}
$$

**Verificación de martingala.** En $u=-i$ debe resultar el *forward*: la parte de saltos vale $\lambda T(1+k-1-k)=0$, y en la de Heston $d=b$, por lo que $g=0$, $D=0$ y $C=(r-q)T$. Así $\varphi_B(-i)=S_0e^{(r-q)T}=F_0$.

> **Relevancia para la investigación:**
> A diferencia de Merton, aquí no hay una serie como (9): al fijar el número de saltos, lo que queda es un Heston y no un Black–Scholes. Por eso en Bates los precios de opciones se calculan solo por Fourier, y $\varphi_B(-i)=F_0$ es la primera prueba numérica en cada calibración.

---

## El variance swap bajo Bates

Repetimos el análisis de Merton: comparamos lo que **paga** el contrato con lo que **replica** la tira de opciones. Todo es igual, solo que $\sigma^2$ se reemplaza por $v_t$.

**Lo que paga el contrato.** Al elevar (18) al cuadrado con las reglas de la nota, la variación cuadrática es $[X]_T=\int_0^Tv_t\,dt+\sum_{i=1}^{N_T}J_i^2$. Su esperanza, con (15), $\mathbb{E}[N_T]=\lambda T$ y $\mathbb{E}[J^2]=\mu_J^2+\sigma_J^2$, es

$$
\boxed{K^{QV}_B=\bar v(T)+\lambda\,(\mu_J^2+\sigma_J^2)} \tag{22}
$$

**Lo que replica la tira de opciones.** Como en Merton, Carr–Madan es una identidad válida con saltos y $\mathbb{E}[\int dS_t/S_{t^-}]=rT$ se mantiene gracias al compensador. Restando (16) y (18) se obtiene $\frac{dS_t}{S_{t^-}}-d\ln S_t=\frac12v_t\,dt+(e^{J}-1-J)\,dN_t$, y al integrar y tomar esperanza:

$$
\boxed{K^{\text{repl}}_B=\bar v(T)+2\lambda\,(k-\mu_J)} \tag{23}
$$

Es la ecuación (13) de Merton con $\sigma^2$ reemplazada por $\bar v(T)$.

**El sesgo no cambia.** Al restar (22) a (23), la parte continua $\bar v(T)$ se cancela:

$$
\boxed{K^{\text{repl}}_B-K^{QV}_B=2\lambda\,\mathbb{E}\Big[e^{J}-1-J-\tfrac12J^2\Big]\approx\frac{\lambda}{3}\,\mathbb{E}[J^3]} \tag{24}
$$

Es el mismo resultado de la sección 1 y de Merton. La conclusión es que **la volatilidad estocástica, por sí sola, no genera sesgo de replicación**: en un Heston puro las trayectorias son continuas y $K^{QV}_H=K^{\text{repl}}_H=\bar v(T)$. El sesgo viene únicamente de los saltos.

**Resumen de los cuatro modelos:**

| Modelo | $K^{QV}$ | $K^{\text{repl}}$ |
|--------|----------|-------------------|
| Black–Scholes | $\sigma^2$ | $\sigma^2$ |
| Merton | $\sigma^2+\lambda(\mu_J^2+\sigma_J^2)$ | $\sigma^2+2\lambda(k-\mu_J)$ |
| Heston | $\bar v(T)$ | $\bar v(T)$ |
| Bates | $\bar v(T)+\lambda(\mu_J^2+\sigma_J^2)$ | $\bar v(T)+2\lambda(k-\mu_J)$ |

> **Relevancia para la investigación:**
> En la práctica los parámetros de salto calibrados serán distintos en Merton y en Bates, porque parte del *smile* que Merton explica con saltos grandes puede quedar absorbida por $\rho$ y $\xi$ en Bates. Por eso lo relevante no es solo cuánto cambia $K$, sino **cuánto del $K$ viene de los saltos** en cada modelo: $2\lambda(k-\mu_J)/K^{\text{repl}}$.

**Ejemplo numérico** (parámetros sintéticos: $T=7/252$, $\kappa=3$, $\theta=0.09$, $\xi=0.5$, $\rho=-0.6$, y saltos $\lambda=36$, $\mu_J=-1\%$, $\sigma_J=6\%$, de modo que $\lambda T=1$). Variamos solo $v_0$; los valores de $K^{\text{repl}}_B$ coinciden, con error relativo del orden de $10^{-6}$, con la integral (4) evaluada sobre precios de Fourier calculados con (21):

| $v_0$ | $\bar v(T)$ | $K^{\text{repl}}_B$ | $K^{QV}_B$ |
|:-----:|:-----------:|:-------------------:|:----------:|
| 0.06 | 0.0612 | 0.1932 | 0.1944 |
| 0.09 | 0.0900 | 0.2220 | 0.2232 |
| 0.12 | 0.1188 | 0.2508 | 0.2520 |

Dos observaciones: la brecha $K^{QV}_B-K^{\text{repl}}_B\approx0.0012$ es **igual en las tres filas** (solo depende de los saltos, como dice (24)); y la fila $v_0=\theta$ coincide con Merton con $\sigma^2=0.09$ y los mismos saltos.

Además, a 7 días la reversión a la media casi no actúa: $\bar v(T)\approx v_0+\tfrac{\kappa(\theta-v_0)}{2}T$, y el peso de $v_0-\theta$ en (15) es $0.96$ a 7 días frente a $0.32$ a un año. Es decir, en near-expiry la parte continua del *strike* es casi $v_0$, y por eso $\kappa$ y $\theta$ quedan mal identificados con un solo vencimiento corto.

---

## Varianza forward y estructura temporal

Con parámetros constantes, $T\,K^{\text{repl}}_B(T)=\int_0^T\mathbb{E}[v_t]\,dt+2\lambda T(k-\mu_J)$. Usando la varianza *forward* de la sección 1, $K_{T_1,T_2}=\frac{T_2K_{T_2}-T_1K_{T_1}}{T_2-T_1}$, queda

$$
K^{\text{repl}}_{B}(T_1,T_2)=\theta+\frac{(v_0-\theta)\big(e^{-\kappa T_1}-e^{-\kappa T_2}\big)}{\kappa\,(T_2-T_1)}+2\lambda\,(k-\mu_J). \tag{25}
$$

La curva $\mathbb{E}[v_t]$ de (15) es **monótona** (va suave de $v_0$ a $\theta$) y la parte de saltos es constante en el tiempo. Por eso un Bates con parámetros constantes **no puede producir un escalón** de varianza en la fecha del anuncio, que es justo lo que se observa en los vencimientos que lo rodean.

Para recuperar el escalón se separa, como en la sección de saltos genéricos y de *earnings*, entre los saltos aleatorios $(\lambda_c,J_c)$ y el salto programado $J_{EA}$ en la fecha conocida $\tau$:

$$
T\,K^{\text{repl}}(T)=\int_0^T\mathbb{E}[v_t]\,dt+2\lambda_cT\,\mathbb{E}\big[e^{J_c}-1-J_c\big]+\mathbb{1}_{\tau\le T}\,2\,\mathbb{E}\big[e^{J_{EA}}-1-J_{EA}\big] . \tag{26}
$$

Para la variación cuadrática es igual, con $J^2$ en lugar de $2(e^{J}-1-J)$. Esta es la descomposición unificada de la sección 1 con $\sigma_t^2$ reemplazada por $v_t$. En la varianza *forward* entre $T_1<\tau\le T_2$, el último término no se cancela y queda concentrado en ese intervalo; así aparece el escalón.

> **Relevancia para la investigación:**
> Un Bates homogéneo solo explica una caída gradual de la varianza (vía reversión a $\theta$). Para reproducir el *vol crush* se necesitan parámetros que dependan del tiempo, por ejemplo una intensidad $\lambda(t)$ concentrada en $\tau$ o un reajuste de $v_t$ en esa fecha. Eso se desarrolla en la sección del evento de *earnings*.

---

## Importancia en la investigación y límites del modelo

**Qué aporta:**

* Separar volatilidad estocástica y saltos, gracias a la tabla de los cuatro modelos y al experimento anidado.
* Dejar claro que el sesgo de replicación (24) viene solo de los saltos, así que un error de pricing de Bates frente al mercado no se puede atribuir a la volatilidad estocástica por sí misma.
* Conectar con la varianza *forward* de la sección 1 a través de (25) y (26).

**Límites:**

* **Identificación:** 8 parámetros y un vencimiento corto; $\kappa$ y $\theta$ casi no influyen a 7 días, así que hay que calibrar varios vencimientos a la vez o fijar algunos parámetros.
* **Un solo salto gaussiano:** igual que Merton, no reproduce densidades bimodales ni curvas de volatilidad implícita cóncavas observadas antes de algunos anuncios (Alexiou, Goyal, Kostakis y Rompolis, 2025; Kachhara, Markin y Singh, 2023).
* **Parámetros constantes en el tiempo:** impiden reproducir el escalón y el *vol crush*, como se vio en (25).
* **Saltos independientes de la varianza:** un salto no modifica la varianza posterior. Los modelos con saltos también en $v_t$ (Duffie, Pan y Singleton, 2000; Eraker, Johannes y Polson, 2003) relajan esto, pero tienden a generar alzas de varianza, no colapsos.




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
| (ii) Poisson de ventana | $N_e\sim\text{Poisson}(\Lambda_e)$, $\Lambda_e=\lambda_e\, delta_e\approx1$ | Mezcla de gaussianas. Es la interpretación de "intensidad alta exactamente en la fecha del anuncio". Genera colas gruesas y smile. |
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

Con parámetros constantes por tramos, el sistema de Riccati (§4.3) se **propaga tramo a tramo** en $\tau=T-t$ con la condición terminal de cada tramo (fórmula cerrada con $D(0)\ne0$ o integración RK4). La integral de intensidad $\Lambda(T)=\int_0^T\lambda(t)\, dt$ sustituye a $\lambda T$ en el término de salto (si $\mu_J,\sigma_J$ son constantes).

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
C(K)=\frac{e^{-\alpha k}}{\pi}\int_0^\infty e^{-i\nu k}\,\psi(\nu)\, d\nu,\qquad
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

* Equivalente en precios: $\sum_iw_i'\big(C_i^{model}-C_i^{mkt}\big)^2$ con $w_i'=1/\text{vega}_i^2$ (o $1/(\text{ask}-\text{bid})^2$), porque $\, delta C\approx\text{vega}\cdot\, delta\sigma$. **El RMSE de precios sin ponderar** sobreajusta las opciones caras (ITM/largas) y falla en wings baratos, que son los que determinan el sesgo por saltos.
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

Para paso $\, delta$ (típicamente 1 día = $1/252$):

$$
\ln S_{t+\, delta}=\ln S_t+\Big(r-q-\tfrac12\sigma^2-\lambda k\Big)\, delta+\sigma\sqrt\, delta\,Z_1+N\mu_J+\sigma_J\sqrt N\,Z_2,\qquad N\sim\text{Poisson}(\lambda\, delta),
$$

con $Z_1,Z_2\sim\mathcal{N}(0,1)$ independientes. La suma de $N$ saltos normales es $\mathcal{N}(N\mu_J,N\sigma_J^2)$, por lo que **no hay error de discretización** (el esquema es exacto en distribución).

#### 8.2 Bates: esquema QE de Andersen (2008)

**Varianza** (Quadratic-Exponential). Dado $v_t$:

$$
m=\theta+(v_t-\theta)e^{-\kappa\, delta},\qquad
s^2=\frac{v_t\xi^2e^{-\kappa\, delta}(1-e^{-\kappa\, delta})}{\kappa}+\frac{\theta\xi^2(1-e^{-\kappa\, delta})^2}{2\kappa},\qquad \psi=\frac{s^2}{m^2}.
$$

* Si $\psi\le\psi_c$ (típicamente $\psi_c=1.5$): $v_{t+\, delta}=a(b+Z_v)^2$, con $b^2=\frac2\psi-1+\sqrt{\frac2\psi}\sqrt{\frac2\psi-1}$ y $a=\frac{m}{1+b^2}$.
* Si $\psi>\psi_c$: con $p=\frac{\psi-1}{\psi+1}$, $\beta=\frac{1-p}{m}$ y $U\sim\mathcal{U}(0,1)$: $v_{t+\, delta}=0$ si $U\le p$; en otro caso $v_{t+\, delta}=\frac1\beta\ln\frac{1-p}{1-U}$.

**Log-precio** (con $\gamma_1=\gamma_2=\tfrac12$):

$$
\ln S_{t+\, delta}=\ln S_t+(r-q-\lambda k)\, delta+K_0+K_1v_t+K_2v_{t+\, delta}+\sqrt{K_3v_t+K_4v_{t+\, delta}}\;Z+N\mu_J+\sigma_J\sqrt N\,Z_2,
$$

$$
K_0=-\frac{\rho\kappa\theta\, delta}{\xi},\quad
K_1=\gamma_1\, delta\Big(\frac{\kappa\rho}{\xi}-\frac12\Big)-\frac\rho\xi,\quad
K_2=\gamma_2\, delta\Big(\frac{\kappa\rho}{\xi}-\frac12\Big)+\frac\rho\xi,\quad
K_3=\gamma_1\, delta(1-\rho^2),\quad K_4=\gamma_2\, delta(1-\rho^2).
$$

Alternativa más simple, aunque sesgada: **full truncation Euler** (Lord, Koekkoek y van Dijk, 2010). Para $\, delta=1/252$ y $\xi$ moderado, QE es notablemente más preciso cerca de $v=0$.

#### 8.3 Ventana de earnings

Se simula sobre la malla diaria hasta el vencimiento ($n$ pasos) con parámetros por tramos:

1. **Pre-EA:** $\lambda_{pre}$ (bajo) y $v_0$ calibrado.
2. **Día del EA** (índice de la fila que contiene el salto según BMO/AMC): $\Lambda_e=\lambda_e\, delta_e\approx1$ (p. ej. $\lambda_e=252$ en un paso diario) **o** salto único con distribución elegida en §5.2.
3. **Post-EA:** $\lambda_{post}$ y, opcionalmente, reset de varianza $v\leftarrow c\,v$ (vol crush).

El compensador del salto en el paso de evento es $-\Lambda_ek_e$ (no $-\lambda k\, delta$ con la intensidad base).

#### 8.4 Estadísticos por trayectoria y estimadores

* $RV^{(m)}=\frac{252}{n}\sum_{i=1}^{n}\big(\ln S^{(m)}_{t_i}/S^{(m)}_{t_{i-1}}\big)^2$ (misma convención que el contrato).
* $\hat K=\frac1M\sum_mRV^{(m)}$ con error estándar $\hat s/\sqrt M$. Contrastar $\hat K$ con: $K^{disc}$ (Merton, exacto), $K^{QV}$ y $K^{\text{repl}}$ analíticos.
* $P\&L^{(m)}=N_{var}\big(RV^{(m)}-K_{var}\big)$: media, desviación, asimetría, cuantiles al 1/5%.
* $\eta^{(m)}=R_{e}^2/\sum_iR_i^2$ (participación del día del evento) y $RV^{(m)}_{ex}$ (RV excluyendo el día del EA).
* Chequeo de **martingala**: $\mathbb{E}[S_T]e^{-(r-q)T}\approx1$.

#### 8.5 Reducción de varianza y semillas

* Variables antitéticas en $Z$; *control variate* usando una variable con esperanza conocida (por ejemplo, $\int v_t\, dt$ aproximada con la regla trapezoidal, cuya esperanza es la de §4.1).
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
