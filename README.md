# Variance Swaps alrededor de Earnings: Marco Teórico con Merton (1976) y Bates (1996)

> Valoración y replicación del *fair strike* de un variance swap bajo saltos y volatilidad estocástica, con aplicación a ventanas de earnings.

## Tabla de contenidos

- [Marco Teórico](#marco-teórico)
  - [1. El variance swap: definición y replicación](#1-el-variance-swap-definición-y-replicación)
  - [2. Modelo de Merton (jump-diffusion)](#2-modelo-de-merton-jump-diffusion)
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

$$(dS_t)^2 = (\mu_t S_t \, dt + \sigma_t S_t \, dW_t)^2$$
$$(dS_t)^2 = \mu_t^2 S_t^2 (dt)^2 + 2\mu_t \sigma_t S_t^2 dt \, dW_t + \sigma_t^2 S_t^2 (dW_t)^2$$

Al aplicar las reglas del cálculo estocástico para cancelar los términos de mayor orden:

$$(dt)^2 = 0, \quad dW_t \, dt = 0, \quad (dW_t)^2 = dt$$

Obtenemos el desarrollo diferencial del término cuadrático:

$$(dS_t)^2 = \sigma_t^2 S_t^2 \, dt$$

> **Nota explicativa sobre el cálculo estocástico empleado:**  
> La identidad del cálculo de Itô $(dW_t)^2 = dt$ difiere del cálculo clásico porque los incrementos del movimiento browniano son tan impredecibles y rápidos que su variación cuadrática en un intervalo pequeño no se anula, sino que converge exactamente al paso del tiempo $dt$.

Sustituyendo todo en la ecuación original del lema de Itô:

$$d(\ln S_t) = \frac{1}{S_t}dS_t - \frac{1}{2}\left(\frac{1}{S_t^2}\right)(\sigma_t^2 S_t^2 \, dt)$$

Simplificando la expresión, se obtiene el retorno porcentual menos un ajuste por convexidad de $-\frac{1}{2}\sigma_t^2 dt$:

$$d(\ln S_t) = \frac{1}{S_t}dS_t - \frac{1}{2}\sigma_t^2 \, dt$$

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
\int_0^{S_*} f''(K)(K-S_T)^+\,dK
$$

La segunda corresponde a una combinación de opciones *call* europeas:

$$
\int_{S_*}^{\infty} f''(K)(S_T-K)^+\,dK
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
f(S_T) = f(S_{\ast}) + f'(S_{\ast})(S_T-S_{\ast}) + \int_0^{S_{\ast}} f''(K)(K-S_T)^+\,dK + \int_{S_{\ast}}^{\infty} f''(K)(S_T-K)^+\,dK
$$

En conjunto, la idea es que un payoff complejo puede descomponerse en una posición constante, una exposición lineal mediante *forwards* o futuros y una combinación de *puts* y *calls* con diferentes *strikes* para reproducir la curvatura del payoff.
En conjunto, la idea es que un payoff complejo puede descomponerse en una posición constante, una exposición lineal mediante *forwards* o futuros y una combinación de *puts* y *calls* con diferentes *strikes* para reproducir la curvatura del payoff.

Con $f(S)=\ln S$ se tiene $f''(K)=-1/K^2$, de modo que

$$
\ln\frac{S_T}{S_{\ast}}=\frac{S_T-S_{\ast}}{S_{\ast}}-\int_0^{S_{\ast}}\frac{(K-S_T)^+}{K^2}\,dK-\int_{S_{\ast}}^{\infty}\frac{(S_T-K)^+}{K^2}\,dK .
$$

Al aplicar el operador de esperanza bajo $\mathbb{Q}$ se usan dos hechos: (i) $\mathbb{E}^{\mathbb{Q}}[S_T]=S_0e^{rT}$ y (ii) $\mathbb{E}^{\mathbb{Q}}\!\left[\int_0^T \frac{dS_t}{S_t}\right]=rT$. Además, el valor esperado no descontado de un pago de opción es su prima capitalizada, $e^{rT}P(K)$ o $e^{rT}C(K)$, donde $P(K)$ y $C(K)$ denotan las primas de las opciones *put* y *call* de *strike* $K$. Sustituyendo en (2) se llega a la fórmula de valoración general:

$$
K_{var} = \frac{2}{T} \left[ rT - \left( \frac{S_0}{S_{\ast}} e^{rT} - 1 \right) - \ln\left(\frac{S_{\ast}}{S_0}\right) + e^{rT} \int_0^{S_{\ast}} \frac{P(K)}{K^2} dK + e^{rT} \int_{S_{\ast}}^\infty \frac{C(K)}{K^2} dK \right] \tag{3}
$$

Si se elige $S_{\ast}=F_0=S_0e^{rT}$ (el *forward*), todos los términos no integrales se cancelan y queda la forma limpia de Britten-Jones y Neuberger (2000):

$$
K_{var}=\frac{2}{T}\,e^{rT}\left[\int_0^{F_0}\frac{P(K)}{K^2}\,dK+\int_{F_0}^{\infty}\frac{C(K)}{K^2}\,dK\right] \tag{4}
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
\frac{dS_t}{S_{t^-}}=(r-\lambda\kappa)\,dt+\sigma_t\,dW_t+\left(e^{J}-1\right)dN_t,\qquad \kappa=\mathbb{E}[e^{J}-1],
$$

donde $N_t$ es un proceso de Poisson de intensidad $\lambda$ y $J$ es el salto del log-precio. El término $-\lambda\kappa\,dt$ compensa los saltos, de modo que $\mathbb{E}^{\mathbb{Q}}\!\left[\int_0^T \frac{dS_t}{S_{t^-}}\right]=rT$ sigue valiendo.

Aplicando Itô con saltos, en cada instante de salto el log-precio cambia en $\Delta X=J$ mientras que $S$ cambia en $e^{\Delta X}-1$ en términos relativos. Por tanto

$$
\ln\frac{S_T}{S_0}=\int_0^T\frac{dS_t}{S_{t^-}}-\frac12\int_0^T\sigma_t^2\,dt-\sum_{s\le T}\left(e^{\Delta X_s}-1-\Delta X_s\right).
$$

La cantidad que replica la tira de opciones es $R_T:=2\left(\int_0^T \frac{dS_t}{S_{t^-}}-\ln\frac{S_T}{S_0}\right)$, y la variación cuadrática es $[X]_T=\int_0^T\sigma_t^2dt+\sum_{s\le T}(\Delta X_s)^2$. Restando:

$$
R_T-[X]_T=2\sum_{s\le T} g(\Delta X_s),\qquad g(x):=e^{x}-1-x-\tfrac12x^2 .
$$

La función $g$ tiene $g(0)=g'(0)=g''(0)=0$ y $g'''(x)=e^{x}>0$, de donde:

* $g(x)=\tfrac16x^3+\tfrac1{24}x^4+\cdots\approx\tfrac16x^3$ para saltos moderados (el error es de **tercer orden**);
* $g(x)$ tiene siempre el **mismo signo que $x$**: un salto bajista siempre resta, uno alcista siempre suma.

### El sesgo de replicación

Tomando esperanza bajo $\mathbb{Q}$ (con $\mathbb{E}\sum g(\Delta X_s)=\lambda T\,\mathbb{E}[g(J)]$):

$$
\boxed{K^{\text{repl}}-K^{QV}=2\lambda\,\mathbb{E}\!\left[e^{J}-1-J-\tfrac12J^2\right]\approx\frac{\lambda}{3}\,\mathbb{E}[J^3]}
$$

donde $K^{\text{repl}}=\mathbb{E}^{\mathbb{Q}}[R_T]/T$ es el valor que entrega la tira de opciones en las ecuaciones (3)-(4) (depende solo de precios de opciones observados) y $K^{QV}=\mathbb{E}^{\mathbb{Q}}\big[[X]_T\big]/T$ es la esperanza de la varianza que el contrato paga en el límite de muestreo fino.

**Lectura.**

* Si los saltos son negativos en media, $\mathbb{E}[J^3]<0$ y $K^{\text{repl}}<K^{QV}$: el *strike* *model-free* de las ecuaciones (3)-(4) **subestima** la varianza cuadrática esperada.
* Lo que el swap paga es $\approx K^{QV}$, **no** $K^{\text{repl}}$. Al comparar modelo contra mercado hay que reportar ambas cantidades y no mezclar definiciones.
* Para quien vende el swap y se cubre con la réplica, el descalce es $R_T-[X]_T$: pierde con saltos bajistas y gana con alcistas, con un efecto cúbico en el tamaño del salto. Es consistente con lo que señalan Carr y Lee (2009) sobre la exposición de estos contratos a retornos al cubo y de orden superior en movimientos bruscos.
* El resultado no es una curiosidad teórica: Broadie y Jain (2008) encuentran que el efecto del muestreo discreto suele ser pequeño mientras que el de los saltos puede ser significativo. Carr, Lee y Lorig (2021) desarrollan réplicas robustas para difusiones con saltos de actividad finita y tamaño acotado.

### Dos fuentes de error que no deben mezclarse

1. **Sesgo por saltos (conceptual).** Ocurre aunque existiera una cadena continua de *strikes* $K\in(0,\infty)$: la réplica de la ecuación (2) deja de ser exacta trayectoria por trayectoria.
2. **Truncamiento y discretización de la cadena (operativo).** Si el subyacente hace un *gap* sobre una región sin *strikes* líquidos, o la cadena real no llega a $K\to0,\infty$, el *strip* discreto tampoco replica bien el contrato logarítmico. Esto ocurre incluso bajo difusión pura, y es peor cuando hay un salto grande (colas gordas) y vencimientos cortos.

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

No obstante, cuando se negocian opciones en los mercados financieros, la cadena no mide la varianza cuadrática de forma directa, sino que incorpora una prima adicional debido a la asimetría de mercado, las colas gordas y el riesgo de caídas abruptas. Esta discrepancia se formaliza mediante la ecuación de replicación, la cual ajusta la varianza realizada sumando los términos de convexidad representados por la función especial de sesgo para cada tipo de salto:

$$_{T}K_T^{\text{repl}} = {}_{T}K_T^{QV} + 2\lambda_c T \mathbb{E}[g(J_c)] + 2\mathbb{1}_{\tau \le T} \mathbb{E}[g(J_{EA})]$$

Dicha función ($g(x) = e^x - 1 - x - \frac{1}{2}x^2$) mide exactamente la distorsión introducida por el contrato logarítmico frente a los movimientos extremos del subyacente[cite: 2].

---

### 4. Aislamiento del salto de *earnings*

Para el salto programado, la intensidad temporal desaparece al tratarse de un evento único en una fecha conocida[cite: 2]. Al aislar el impacto específico de un reporte de resultados, los componentes de difusión y los saltos aleatorios se neutralizan[cite: 2]. Dado que un anuncio corporativo se concentra en un único retorno discreto entre el cierre previo y la apertura posterior, su efecto estructural se aproxima mediante una expansión matemática expresada como[cite: 2]:

$$K_T^{\text{repl}} - K_T^{QV}|_{EA} = \frac{2}{T} \mathbb{E}[g(J_{EA})] \approx \frac{1}{3T} \mathbb{E}[J_{EA}^3]$$

Esto demuestra analíticamente que el sesgo generado por los *earnings* en la estructura de volatilidad está íntimamente relacionado con la asimetría direccional (momento cúbico) del precio ante la revelación de la información[cite: 2].

### Comparación conceptual

| | Salto genérico (Poisson) | Salto de earnings (programado) |
|---|---|---|
| Momento | Aleatorio, intensidad $\lambda_c$ | Fecha $\tau$ conocida |
| Tamaño | Aleatorio; típicamente sesgo negativo (riesgo de caída) | Aleatorio y grande; signo *a priori* incierto |
| Aporte a $K^{QV}$ | $\lambda_c\,\mathbb{E}[J_c^2]$, repartido de forma uniforme en el horizonte | $\mathbb{E}[J_{EA}^2]/T$, concentrado en una fecha |
| Sesgo de la réplica | $2\lambda_c\,\mathbb{E}[g(J_c)]$ | $\tfrac{2}{T}\,\mathbb{E}[g(J_{EA})]$ |
| Cómo se identifica con opciones | Solo vía *smile* (asimetría y curtosis), mezclado con la difusión | Vía estructura temporal: escalón de varianza entre vencimientos que rodean a $\tau$ |
| Efecto observable | Suaviza la superficie | La volatilidad implícita sube antes del anuncio y cae al conocerse; estructura temporal decreciente a corto plazo |

La consecuencia central: **el salto de earnings se puede aislar con la estructura temporal de varianza** (sección anterior) mientras que el salto genérico solo se distingue de la difusión a través de la forma del *smile*. La literatura de earnings aprovecha esto: Dubinsky, Johannes, Kaeck y Seeger (2019) separan la incertidumbre del anuncio de la volatilidad diaria normal comparando volatilidad antes y después del anuncio o usando la estructura temporal de volatilidades implícitas.

### Extracción del salto de earnings con la varianza forward

Sea $T_1<\tau<T_2$ y sea $\bar v$ la varianza forward "base" (sin evento), estimada por ejemplo interpolando las varianzas forward de las ventanas contiguas que no contienen $\tau$. Como $w^{QV}(T_2)-w^{QV}(T_1)=(T_2-T_1)\bar v+\mathbb{E}[J_{EA}^2]$:

$$
\mathbb{E}^{\mathbb{Q}}[J_{EA}^2]\approx(T_2-T_1)\left(K^{\text{repl}}_{T_1,T_2}-\bar v\right)-\underbrace{2\,\mathbb{E}[g(J_{EA})]}_{\approx\,\mathbb{E}[J_{EA}^3]/3}
$$

El último término corrige el sesgo de replicación cuando $K_{T_1,T_2}$ se calcula con las ecuaciones (3)-(4); si la ventana es corta, este salto domina el resultado. Esta es una medida bajo $\mathbb{Q}$, por lo que incluye la prima de riesgo por el salto: no es directamente la esperanza física del movimiento.

### Magnitud del sesgo (ilustración con cálculo propio)

El sesgo relativo al aporte del salto a la varianza es $\;2\,\mathbb{E}[g(J)]\,/\,\mathbb{E}[J^2]$, independiente de $T$:

| Caso | Sesgo relativo |
|---|---|
| Earnings simétrico con $\sigma_J=6\%$ y asimetría $-0.5$ (aprox. cúbica: $\text{asim.}\cdot\sigma_J/3$) | $\approx -1.0\%$ |
| Salto fijo $J=-10\%$ | $\approx -3.3\%$ |
| Salto fijo $J=-20\%$ | $\approx -6.3\%$ |

Lo que sí depende de $T$ es cuánto pesa el salto en el *strike*: entra como $\mathbb{E}[J^2]/T$. Por ejemplo, con $T=10$ días hábiles ($10/252$), un salto de earnings con $\mathbb{E}[J_{EA}^2]=(6\%)^2$ suma $\approx0.091$ en varianza, de modo que un *strike* de 25% de volatilidad base pasa a $\approx39\%$; el sesgo por asimetría en ese caso es del orden de $0.1$ puntos de volatilidad.

### Implicaciones para el análisis

1. **El salto de earnings entra completo en lo que paga el swap** ($\mathbb{E}[J_{EA}^2]$), y la tira de opciones lo captura, salvo por el error de tercer orden $\tfrac{1}{3T}\mathbb{E}[J_{EA}^3]$. En vencimientos cortos que contienen el anuncio, el componente de earnings domina el *strike* y el sesgo relativo, aunque pequeño en porcentaje, deja de ser despreciable en términos absolutos.
2. **Los saltos genéricos y los de earnings no deben mezclarse** en el modelo: los primeros se estiman con el *smile*, los segundos con el escalón de varianza en la estructura temporal.
3. **Reportar siempre** $K^{\text{repl}}$ (lo que dan las opciones) y $K^{QV}$ (lo que paga el swap bajo el modelo), con y sin extrapolación de las colas, para separar el sesgo por saltos del error operativo de truncamiento.
4. Para vencimientos que **no** contienen el anuncio ($T<\tau$) el término de earnings es cero y la fórmula (4) es una buena aproximación salvo por los saltos genéricos.

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
