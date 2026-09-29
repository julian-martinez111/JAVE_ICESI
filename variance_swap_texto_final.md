# 1. El variance swap: definición y replicación

## Estructura, mecánica y liquidación

Un *variance swap* es un contrato forward sobre la varianza anualizada de un activo subyacente específico. Su pago al vencimiento $T$ es lineal respecto a la varianza realizada:

$$\text{Payoff} = (\sigma_{R}^{2} - K_{\text{var}}) \times N_{\text{vol}}$$

Donde $\sigma_{R}^{2}$ denota la varianza realizada del activo subyacente durante el horizonte del contrato, $K_{\text{var}}$ representa el *strike* de varianza implícita y $N_{\text{vol}}$ corresponde al valor nocional expresado en unidades monetarias por punto de varianza.

---

## Derivación matemática del strike de varianza ($K_{\text{var}}$)

El objetivo central de la modelización cuantitativa consiste en determinar el precio de ejercicio justo ($K_{\text{var}}$), el cual, bajo condiciones de ausencia de arbitraje, debe igualar la expectativa matemática de la varianza realizada bajo la medida neutral al riesgo.

### El proceso de precios continuo

Se asume que la trayectoria del precio del subyacente $S_t$ sigue una difusión de Itô bajo la medida histórica $\mathbb{P}$, expresada en términos de su retorno instantáneo (el movimiento browniano geométrico es el caso particular con $\mu$ y $\sigma$ constantes; aquí ambos pueden ser procesos adaptados):

$$\frac{dS_t}{S_t} = \mu_t dt + \sigma_t dW_t$$

Donde $\mu_t$ denota el *drift*, $\sigma_t$ representa la volatilidad instantánea y $dW_t$ corresponde al incremento de un proceso de Wiener estándar. La varianza realizada anualizada $V$ se formaliza como:

$$V = \frac{1}{T}\int_{0}^{T}\sigma_{t}^{2}dt$$

### Aplicación del Lema de Itô y transformación logarítmica

El desafío de esta valoración radica en aislar la varianza instantánea $\sigma_t^2$ utilizando únicamente variables de mercado observables (en este caso, el precio spot del activo). Puesto que $\sigma_t$ no es una variable directamente observable en el mercado, se recurre a una transformación logarítmica mediante la aplicación del Lema de Itô sobre la función $f(S_t) = \ln(S_t)$.

Para aislar este término del proceso del precio del activo, partimos de las derivadas de la función: $f'(S) = \frac{1}{S}$ y $f''(S) = -\frac{1}{S^2}$, las cuales permiten cancelar el término $S^2$ al aplicar el lema de Itô:

$$d(\ln S_t) = f'(S_t)dS_t + \frac{1}{2}f''(S_t)(dS_t)^2$$

Desarrollando el término $(dS_t)^2$ a partir de la dinámica del precio:

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

$$d(\ln S_t) = \frac{1}{S_t}dS_t - \frac{1}{2}\sigma_t^2 \, dt \tag{1}$$

---

Reorganizando (1), se aísla el término de la varianza $\sigma_t^2 dt$. Tras integrar la expresión a lo largo del intervalo $[0, T]$, se demuestra que la varianza realizada se refleja exactamente manteniendo una posición en acciones rebalanceada continuamente (con $\frac{2}{T S_t}$ acciones en cada instante) y tomando una posición corta estática en un instrumento que pague $\frac{2}{T}\ln\left(\frac{S_T}{S_0}\right)$ (un *contrato logarítmico*):

$$V = \frac{2}{T}\left[\int_{0}^{T}\frac{dS_t}{S_t} - \ln\left(\frac{S_T}{S_0}\right)\right] \tag{2}$$

**Importante:** la identidad (2) es válida *trayectoria por trayectoria* solo si el precio tiene trayectorias continuas. Esta hipótesis se relaja en la sección "Límite de la réplica: saltos".

---

## Esperanza neutral al riesgo y el Teorema de Carr-Madan

Para encontrar un instrumento negociable en el mercado que sea equivalente al pago logarítmico, la valoración se traslada al marco de la probabilidad neutral al riesgo ($\mathbb{Q}$). Bajo esta medida, el activo subyacente crece en promedio a la tasa libre de riesgo $r$, permitiendo descontar flujos esperados a valor presente sin sesgos de preferencias de riesgo subjetivas. Así, el *strike* teórico del *variance swap* se define como:

$$K_{\text{var}} = \mathbb{E}^{\mathbb{Q}}\left[ \frac{1}{T}\int_{0}^{T}\sigma_{t}^{2}dt \right]$$

Para eliminar la necesidad de negociar contratos logarítmicos sintéticos, se implementa el teorema de replicación estática desarrollado por Peter Carr y Dilip Madan (1998). Este teorema demuestra que cualquier pago dos veces diferenciable $f(S_T)$ puede descomponerse analíticamente mediante una expansión de Taylor ponderada sobre una cadena infinita de opciones europeas *vanilla* (*calls* y *puts*).

### Desglose analítico del Teorema de Carr-Madan

* **Término de posición constante:**

$$
f(S_{\ast})
$$

evaluado en un precio de referencia preestablecido $S_{\ast}$. Usualmente, $S_{\ast}$ se fija como el precio *forward* vigente del activo. Este término representa simplemente un valor constante que no depende del precio futuro $S_T$.

* **Término de exposición lineal:**

$$
f'(S_{\ast})(S_T-S_{\ast})
$$

Este componente representa la parte lineal del payoff, ya que depende directamente de $S_T$. Debido a esta característica, puede replicarse mediante una posición estática en contratos *forward* o futuros sobre el activo subyacente.

* **Componentes convexos de opciones:** son los términos encargados de representar la curvatura del payoff y están dados por dos integrales.

La primera corresponde a una combinación de opciones *put* europeas:

$$
\int_0^{S_{\ast}} f''(K)(K-S_T)^+\,dK
$$

La segunda corresponde a una combinación de opciones *call* europeas:

$$
\int_{S_{\ast}}^{\infty} f''(K)(S_T-K)^+\,dK
$$

Esto se debe a que $(K-S_T)^+$ es exactamente el *payoff* de una *put* europea, mientras que $(S_T-K)^+$ es exactamente el *payoff* de una *call* europea.

Por lo tanto, la representación completa de Carr-Madan puede expresarse como:

$$
f(S_T) = f(S_{\ast}) + f'(S_{\ast})(S_T-S_{\ast}) + \int_0^{S_{\ast}} f''(K)(K-S_T)^+\,dK + \int_{S_{\ast}}^{\infty} f''(K)(S_T-K)^+\,dK
$$

En conjunto, la idea es que un payoff complejo puede descomponerse en una posición constante, una exposición lineal mediante *forwards* o futuros y una combinación de *puts* y *calls* con diferentes *strikes* para reproducir la curvatura del payoff.

### Cierre de la derivación: aplicación a $f(S)=\ln S$

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

Como la varianza (no la volatilidad) es aditiva en el tiempo, la varianza total esperada $w(T):=T\,K_T$ es aditiva y, entre dos vencimientos $T_1<T_2$, el *strike* de un swap con inicio a futuro es

$$
K_{T_1,T_2}=\frac{T_2K_{T_2}-T_1K_{T_1}}{T_2-T_1}=\frac{w(T_2)-w(T_1)}{T_2-T_1}.
$$

Se replica comprando un variance swap a $T_2$ con nocional en varianza $\frac{T_2}{T_2-T_1}$ y vendiendo uno a $T_1$ con nocional $\frac{T_1}{T_2-T_1}$: la diferencia de pagos es $RV_{T_1,T_2}-K_{T_1,T_2}$. Es la herramienta *model-free* para leer la **estructura temporal de varianza** y localizar en el calendario dónde se concentra la varianza esperada.

---

## Saltos genéricos vs. saltos de earnings

La sección de saltos supone que estos llegan en instantes aleatorios (Poisson). Un anuncio de resultados (*earnings announcement*, EA) es distinto: la **fecha $\tau$ se conoce con antelación**, aunque no se conozca la respuesta del precio. Esto cambia qué se puede identificar y cómo se lee el sesgo.

### Descomposición unificada

Sea $J_c$ el salto genérico (Poisson, intensidad $\lambda_c$) y $J_{EA}$ el salto programado en la fecha $\tau$. Para un vencimiento $T$:

$$
T\,K^{QV}_T=\underbrace{\int_0^T \mathbb{E}^{\mathbb{Q}}[\sigma_t^2]\,dt}_{\text{difusión}}+\underbrace{\lambda_c T\,\mathbb{E}[J_c^2]}_{\text{saltos genéricos}}+\underbrace{\mathbb{1}_{\{\tau\le T\}}\,\mathbb{E}[J_{EA}^2]}_{\text{earnings}}
$$

y la tira de opciones entrega lo mismo más el sesgo de la sección anterior:

$$
T\,K^{\text{repl}}_T=T\,K^{QV}_T+2\lambda_cT\,\mathbb{E}[g(J_c)]+2\,\mathbb{1}_{\{\tau\le T\}}\,\mathbb{E}[g(J_{EA})].
$$

Para el salto programado desaparece $\lambda T$: es un único evento en fecha conocida, con $g(x)=e^x-1-x-\tfrac12x^2$ como antes. Su sesgo en unidades anualizadas es

$$
K^{\text{repl}}_T-K^{QV}_T\Big|_{EA}=\frac{2}{T}\,\mathbb{E}[g(J_{EA})]\approx\frac{1}{3T}\,\mathbb{E}[J_{EA}^3].
$$

Como el swap discreto mide retornos diarios, el salto de earnings queda contenido en **un solo retorno** (el del cierre previo a la apertura posterior al anuncio) y entra completo, con su $J_{EA}^2$, en la varianza realizada.

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

## Referencias

* Britten-Jones, M., y Neuberger, A. (2000). Option prices, implied price processes, and stochastic volatility. *Journal of Finance*, 55(2).
* Broadie, M., y Jain, A. (2008). The effect of jumps and discrete sampling on volatility and variance swaps. *International Journal of Theoretical and Applied Finance*, 11(8), 761-797.
* Carr, P., y Lee, R. (2009). Volatility derivatives. *Annual Review of Financial Economics*, 1, 319-339.
* Carr, P., Lee, R., y Lorig, M. (2021). Robust replication of volatility and hybrid derivatives on jump diffusions. *Mathematical Finance*.
* Carr, P., y Madan, D. (1998). Towards a theory of volatility trading. En *Volatility: New Estimation Techniques for Pricing Derivatives*.
* Demeterfi, K., Derman, E., Kamal, M., y Zou, J. (1999). More than you ever wanted to know about volatility swaps. Goldman Sachs Quantitative Strategies Research Notes.
* Dubinsky, A., Johannes, M., Kaeck, A., y Seeger, N. J. (2019). Option pricing of earnings announcement risks. *Review of Financial Studies*, 32(2), 646-687.
* Dupire, B. (1993). Model art. *Risk*, 6(9).
* Jarrow, R., Kchia, Y., Larsson, M., y Protter, P. (2013). Discretely sampled variance and volatility swaps versus their continuous approximations. *Finance and Stochastics*, 17(2), 305-324.
* Leung, T., y Santoli, M. (2014). Accounting for earnings announcements in the pricing of equity options. arXiv:1412.8414.
* Neuberger, A. (1994). The log contract. *Journal of Portfolio Management*, 20(2).
* Patell, J., y Wolfson, M. (1979). Anticipated information releases reflected in call option prices. *Journal of Accounting and Economics*, 1(2), 117-140.
