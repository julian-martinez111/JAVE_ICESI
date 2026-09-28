# Variance Swaps alrededor de Earnings: Marco Teórico con Merton (1976) y Bates (1996)

> **Fair strike teórico vs. mercado en contratos near-expiry con anuncio de resultados.**
> Replicación estática (log-contract), calibración por transformadas de Fourier, simulación Monte Carlo de la ventana de earnings y análisis del *vol crush*.

---

## Tabla de contenidos

1. [Resumen, preguntas de investigación e hipótesis](#1-resumen-preguntas-de-investigación-e-hipótesis)
2. [Preliminares de cálculo estocástico](#2-preliminares-de-cálculo-estocástico)
3. [El variance swap: definición, replicación y sesgos](#3-el-variance-swap-definición-replicación-y-sesgos)
4. [Modelo de Merton (jump-diffusion)](#4-modelo-de-merton-jump-diffusion)
5. [Modelo de Heston y modelo de Bates (SV + saltos)](#5-modelo-de-heston-y-modelo-de-bates-sv--saltos)
6. [Modelar el evento de earnings y el vol crush](#6-modelar-el-evento-de-earnings-y-el-vol-crush)
7. [Valoración de opciones por transformadas de Fourier](#7-valoración-de-opciones-por-transformadas-de-fourier)
8. [Calibración](#8-calibración)
9. [Simulación Monte Carlo](#9-simulación-monte-carlo)
10. [Diseño empírico y métricas](#10-diseño-empírico-y-métricas)
11. [Checklist de consideraciones prácticas](#11-checklist-de-consideraciones-prácticas)
12. [Ejemplo numérico sintético (reproducible)](#12-ejemplo-numérico-sintético-reproducible)
13. [Referencia de implementación en Python](#13-referencia-de-implementación-en-python)
14. [Estructura sugerida del repositorio](#14-estructura-sugerida-del-repositorio)
15. [Limitaciones y extensiones](#15-limitaciones-y-extensiones)
16. [Bibliografía](#16-bibliografía)

**Notación.** $S_t$ precio del subyacente, $r$ tasa libre de riesgo, $q$ dividend yield continuo, $T$ vencimiento (en años de *tiempo de negociación*, ver §11), $F_0=S_0e^{(r-q)T}$ forward, $X_t=\ln S_t$, $\mathbb{Q}$ medida neutral al riesgo, $\mathbb{P}$ medida física, $t_e$ fecha del anuncio de resultados.

---

## 1. Resumen, preguntas de investigación e hipótesis

### 1.1 Motivación

Un anuncio de resultados (*earnings announcement*, EA) es un evento de información **programado** cuya fecha se conoce, pero cuyo contenido no. Los precios de opciones de vencimiento corto que "cruzan" el anuncio incorporan una varianza extra concentrada en un solo día (la *varianza de evento*), que desaparece al resolverse la incertidumbre (*vol crush*). Dubinsky, Johannes, Kaeck y Seeger (2019) muestran que esa incertidumbre anticipada es cuantitativamente grande, varía en el tiempo y es informativa sobre la volatilidad futura, y proponen modelos reducidos para separarla de la volatilidad "normal" diaria. Barth y So (2014) documentan que parte de esa varianza lleva prima de riesgo cuando el EA implica riesgo de volatilidad no diversificable.

El **variance swap** es el instrumento natural para aislar esta magnitud: su *fair strike* $K_{var}$ es la esperanza neutral al riesgo de la varianza realizada y, bajo difusión pura, se replica con una tira estática de opciones (Neuberger, 1994; Dupire, 1993; Demeterfi–Derman–Kamal–Zou, 1999). Con saltos, la equivalencia entre "lo que cuesta la tira de opciones" y "la varianza realizada esperada" **deja de ser exacta** (§3.5). Ese sesgo es precisamente lo que un modelo de saltos como Merton o Bates permite cuantificar.

### 1.2 Diseño en una línea

```
Option chain near-expiry ──► Calibración (Merton, Bates) ──► K_var teórico (Fourier, forma cerrada)
        │                                                          │
        └──► K_var model-free de mercado (strip de opciones) ◄─────┘   ← Comparación / error de pricing
                                   │
                    Monte Carlo (SDE con salto de earnings) ──► RV simulada, vol crush, P&L del swap
```

### 1.3 Preguntas de investigación

| # | Pregunta |
|---|----------|
| RQ1 | ¿Qué modelo (Merton o Bates) reproduce mejor el *fair strike* implícito en las opciones reales antes del anuncio? |
| RQ2 | ¿Qué fracción de $K_{var}$ se debe al riesgo de salto puro (Merton) y cuánto cambia al añadir volatilidad estocástica (Bates)? |
| RQ3 | ¿Cuán grande es el sesgo de replicación por saltos, $K^{\text{repl}}-K^{QV}$, en ventanas de earnings? |
| RQ4 | ¿Cómo reacciona la estructura temporal de varianza de Bates ante el colapso post-anuncio (*vol crush*)? |
| RQ5 | ¿Es $K_{var}$ un predictor insesgado de la varianza realizada en la ventana (prima de riesgo de varianza de evento)? |

### 1.4 Hipótesis contrastables

* **H1.** A horizontes de 5–10 días, la mayor parte de la curvatura/skew del smile la explica el componente de salto; Bates mejora de forma marginal a Merton *intra-vencimiento*, pero su ventaja aumenta al calibrar conjuntamente varios vencimientos (pre y post EA).
* **H2.** Con saltos de media negativa, $K^{\text{repl}}<K^{QV}$: el strike model-free de mercado **subestima** la varianza cuadrática esperada; el sesgo crece con $\lambda\,\mathbb{E}[J^3]$.
* **H3.** La varianza de evento representa la mayor parte de $K_{var}$ en vencimientos que cruzan el anuncio.
* **H4.** Un Bates homogéneo en el tiempo **no** puede reproducir el escalón de varianza a término asociado al EA sin parámetros dependientes del tiempo ($\lambda(t)$, $\theta(t)$).
* **H5.** $K_{var}>\mathbb{E}^{\mathbb{P}}[RV]$ en promedio, con heterogeneidad transversal (Barth y So, 2014).

> Las hipótesis se formulan para ser contrastadas, no como resultados establecidos.

---

## 2. Preliminares de cálculo estocástico

### 2.1 Espacio de probabilidad y procesos base

Sea $(\Omega,\mathcal{F},(\mathcal{F}_t),\mathbb{Q})$ un espacio filtrado que soporta:

* $W_t$: movimiento browniano estándar.
* $N_t$: proceso de Poisson de intensidad $\lambda$, independiente de $W$.
* $(J_i)_{i\ge1}$: saltos iid en logaritmo, $J\sim\mathcal{N}(\mu_J,\sigma_J^2)$, independientes de $W$ y $N$.

El **proceso de Poisson compuesto** $\sum_{i=1}^{N_t}J_i$ tiene

$$
\mathbb{E}\Big[\sum_{i=1}^{N_t}J_i\Big]=\lambda t\,\mu_J,\qquad
\operatorname{Var}\Big[\sum_{i=1}^{N_t}J_i\Big]=\lambda t\,(\mu_J^2+\sigma_J^2).
$$

### 2.2 Fórmula de Itô con saltos

Si $dX_t=a_t\,dt+b_t\,dW_t+\Delta X_t\,dN_t$ y $f\in C^{1,2}$,

$$
df(t,X_t)=\Big(f_t+a_tf_x+\tfrac12 b_t^2 f_{xx}\Big)dt+b_tf_x\,dW_t+\big[f(t,X_{t^-}+\Delta X_t)-f(t,X_{t^-})\big]dN_t .
$$

### 2.3 Variación cuadrática

Para $X$ con parte continua $X^c$ y saltos $\Delta X_s$:

$$
[X]_T=\langle X^c\rangle_T+\sum_{0<s\le T}(\Delta X_s)^2 .
$$

Para el log-precio, $[X]_T=\int_0^T\sigma_t^2\,dt+\sum_{s\le T}(\Delta X_s)^2$. Además, para particiones $0=t_0<\dots<t_n=T$ con malla $\to0$,

$$
\sum_{i=1}^n\big(X_{t_i}-X_{t_{i-1}}\big)^2\ \xrightarrow{\ \mathbb{P}\ }\ [X]_T,
$$

que es el fundamento de la **varianza realizada** como estimador de la variación cuadrática (Barndorff-Nielsen y Shephard, 2004).

### 2.4 Medida neutral al riesgo y compensación de saltos

Bajo $\mathbb{Q}$, el precio descontado con dividendos, $e^{-(r-q)t}S_t$, es martingala. Con saltos $e^{J}-1$, la deriva debe **compensar** el salto esperado:

$$
\frac{dS_t}{S_{t^-}}=(r-q-\lambda k)\,dt+\sigma\,dW_t+(e^{J}-1)\,dN_t,\qquad
k:=\mathbb{E}[e^{J}-1]=e^{\mu_J+\frac12\sigma_J^2}-1 .
$$

Notas:

* Los parámetros de salto $(\lambda,\mu_J,\sigma_J)$ **calibrados a opciones son parámetros bajo $\mathbb{Q}$**; incorporan la prima de riesgo de salto (Pan, 2002; Bates, 2000). Merton (1976) asumió riesgo de salto diversificable ($\mathbb{Q}=\mathbb{P}$ para el salto); la evidencia empírica no lo respalda en general.
* Bajo la medida física, $\lambda^{\mathbb{P}},\mu_J^{\mathbb{P}},\sigma_J^{\mathbb{P}}$ pueden diferir; la diferencia contiene información sobre la prima de riesgo de varianza y de salto.

### 2.5 Función característica

Para $X_T=\ln S_T$, $\varphi(u)=\mathbb{E}^{\mathbb{Q}}[e^{iuX_T}]$. Los modelos afines (Duffie, Pan y Singleton, 2000) tienen $\varphi$ en forma exponencial-afín, lo que permite valorar opciones por inversión de Fourier (§7). Propiedades usadas en este proyecto:

* **Martingala:** $\varphi(-i)=\mathbb{E}[S_T]=F_0$ (prueba numérica de consistencia).
* **Momentos:** $\mathbb{E}[X^n]=i^{-n}\varphi^{(n)}(0)$ (útil para la varianza realizada *discreta*, §3.4).
* **Independencia:** si el proceso es suma de componentes independientes, $\varphi$ es el producto de las funciones características (Bates = Heston × salto).

---

## 3. El variance swap: definición, replicación y sesgos

### 3.1 Contrato

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

### 3.2 Del contrato a la esperanza

Bajo $\mathbb{Q}$, con tasas deterministas, $K_{var}=\mathbb{E}^{\mathbb{Q}}[RV]$ (ignorando el descuento del pago único). La cuestión es **qué versión de $RV$**: (i) la variación cuadrática continua $[X]_T/T$, (ii) la suma discreta de cuadrados de log-retornos, (iii) la cantidad replicable por opciones. Bajo difusión pura las tres coinciden; con saltos no.

### 3.3 Replicación estática (Neuberger, Dupire, Carr–Madan)

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

### 3.4 Discretización y truncamiento con strikes de mercado

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

### 3.5 Resultado clave: el sesgo de replicación por saltos

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
* En un modelo de saltos **el swap (que paga $RV$ discreta) vale $\approx K^{QV}$, no $K^{\text{repl}}$**. Por eso el análisis debe reportar ambos y no mezclar definiciones al comparar modelo vs. mercado (RQ1 y RQ3).

Este punto es tratado sistemáticamente por Broadie y Jain (2008) y Carr y Lee (2009); Jarrow, Kchia, Larsson y Protter (2013) estudian la aproximación de la versión discreta por su límite continuo.

### 3.6 Variance swap vs. volatility swap (convexidad)

Como $\sqrt{\cdot}$ es cóncava, por Jensen $\mathbb{E}[\sqrt{RV}]<\sqrt{\mathbb{E}[RV]}$ y, a segundo orden,

$$
K_{vol}\approx\sqrt{K_{var}}-\frac{\operatorname{Var}(RV)}{8\,K_{var}^{3/2}} .
$$

Con saltos, $\operatorname{Var}(RV)$ es muy grande y esta aproximación puede fallar (Broadie y Jain, 2008); el análisis del proyecto se centra en varianza, que es replicable, y en $\sqrt{K_{var}}$ solo como unidad de reporte en puntos de volatilidad.

### 3.7 Varianza *forward* (calendar) e implicaciones para earnings

Por aditividad, entre dos vencimientos $T_1<T_2$:

$$
K_{T_1,T_2}=\frac{T_2K_{T_2}-T_1K_{T_1}}{T_2-T_1}.
$$

Es la herramienta model-free para leer la **estructura temporal de varianza** y localizar el salto de varianza que introduce el EA (§6).

---

## 4. Modelo de Merton (jump-diffusion)

### 4.1 SDE

Bajo $\mathbb{Q}$ (§2.4):

$$
\frac{dS_t}{S_{t^-}}=(r-q-\lambda k)\,dt+\sigma\,dW_t+(e^{J}-1)\,dN_t,\qquad J\sim\mathcal{N}(\mu_J,\sigma_J^2).
$$

**Parámetros (4):** $\sigma$ (volatilidad de difusión), $\lambda$ (intensidad anual de saltos), $\mu_J$ (media del salto en log), $\sigma_J$ (desviación del salto).

### 4.2 Solución y momentos

$$
\ln\frac{S_T}{S_0}=\Big(r-q-\tfrac12\sigma^2-\lambda k\Big)T+\sigma W_T+\sum_{i=1}^{N_T}J_i .
$$

$$
\mathbb{E}[R]=\Big(r-q-\tfrac12\sigma^2-\lambda k+\lambda\mu_J\Big)T,\qquad
\operatorname{Var}[R]=\big(\sigma^2+\lambda(\mu_J^2+\sigma_J^2)\big)T,
$$

y cumulantes de orden $n\ge3$: $\kappa_n=\lambda T\,\mathbb{E}[J^n]$. La **curtosis en exceso** es $\lambda T\,\mathbb{E}[J^4]/\operatorname{Var}[R]^2$ y decae como $1/(\lambda T)$: para ventanas de earnings con $\lambda T=O(1)$, las colas son muy gruesas.

### 4.3 Función característica

$$
\varphi_M(u)=\exp\Big\{iu\big[\ln S_0+(r-q-\tfrac12\sigma^2-\lambda k)T\big]-\tfrac12\sigma^2u^2T+\lambda T\big(e^{iu\mu_J-\frac12\sigma_J^2u^2}-1\big)\Big\}.
$$

### 4.4 Precio en serie de Merton

Condicionando en $N_T=n$ el log-precio es normal, de modo que

$$
C_M=\sum_{n=0}^{\infty}\frac{e^{-\lambda'T}(\lambda'T)^n}{n!}\;C_{BS}\big(S_0,K,T,\sigma_n,r_n,q\big),
$$

$$
\lambda'=\lambda(1+k),\qquad \sigma_n^2=\sigma^2+\frac{n\sigma_J^2}{T},\qquad r_n=r-\lambda k+\frac{n\ln(1+k)}{T}.
$$

Es una **mezcla de Poisson de modelos Black–Scholes**: con $\lambda T\approx1$ pocas iteraciones ($n\le 15$) bastan y da una validación independiente del método de Fourier.

### 4.5 Fair strike bajo Merton

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

**Descomposición del riesgo de salto puro (RQ2).** Definiendo $K^{BS}=\sigma^2$:

$$
\underbrace{K^{\text{repl}}_M-\sigma^2}_{\text{aporte del salto en el strike de mercado}}=2\lambda\big(k-\mu_J\big),\qquad
\text{Sesgo}=K^{\text{repl}}_M-K^{QV}_M=2\lambda\,\mathbb{E}\big[e^J-1-J-\tfrac12J^2\big].
$$

---

## 5. Modelo de Heston y modelo de Bates (SV + saltos)

### 5.1 Heston (1993)

$$
\begin{aligned}
\frac{dS_t}{S_t}&=(r-q)\,dt+\sqrt{v_t}\,dW_t^{S},\\
dv_t&=\kappa(\theta-v_t)\,dt+\xi\sqrt{v_t}\,dW_t^{v},\qquad d\langle W^S,W^v\rangle_t=\rho\,dt .
\end{aligned}
$$

**Parámetros (5):** $v_0,\kappa,\theta,\xi,\rho$. (Bajo $\mathbb{Q}$, $\kappa=\kappa^{\mathbb{P}}+\lambda_v$ y $\theta=\kappa^{\mathbb{P}}\theta^{\mathbb{P}}/\kappa$ si la prima de riesgo de varianza es $\lambda_v v$.)

* **Condición de Feller:** $2\kappa\theta\ge\xi^2$ garantiza $v_t>0$. Las calibraciones a opciones de corto plazo frecuentemente la violan; es aceptable si el esquema de simulación maneja $v=0$ (QE, §9).
* **Esperanza de la varianza instantánea:** $\mathbb{E}[v_t]=\theta+(v_0-\theta)e^{-\kappa t}$.
* **Strike bajo Heston:**

$$
K^{QV}_H=K^{\text{repl}}_H=\frac1T\,\mathbb{E}\!\int_0^Tv_t\,dt=\theta+\frac{(v_0-\theta)\big(1-e^{-\kappa T}\big)}{\kappa T}.
$$

### 5.2 Bates (1996)

$$
\begin{aligned}
\frac{dS_t}{S_{t^-}}&=(r-q-\lambda k)\,dt+\sqrt{v_t}\,dW_t^{S}+(e^{J}-1)\,dN_t,\\
dv_t&=\kappa(\theta-v_t)\,dt+\xi\sqrt{v_t}\,dW_t^{v},\qquad d\langle W^S,W^v\rangle_t=\rho\,dt,
\end{aligned}
$$

con $N_t$ y $J\sim\mathcal{N}(\mu_J,\sigma_J^2)$ independientes de $(W^S,W^v)$. **Parámetros (8):** $v_0,\kappa,\theta,\xi,\rho,\lambda,\mu_J,\sigma_J$.

### 5.3 Función característica (formulación estable, "Little Heston Trap")

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

### 5.4 Fair strike bajo Bates

Como la parte de difusión y la parte de salto son independientes,

$$
K^{QV}_B=\underbrace{\theta+\frac{(v_0-\theta)(1-e^{-\kappa T})}{\kappa T}}_{\text{volatilidad estocástica}}+\underbrace{\lambda\,(\mu_J^2+\sigma_J^2)}_{\text{saltos (QV)}},
$$

$$
K^{\text{repl}}_B=\theta+\frac{(v_0-\theta)(1-e^{-\kappa T})}{\kappa T}+2\lambda\big(k-\mu_J\big).
$$

(Demostración de $K^{\text{repl}}_B$: $\mathbb{E}\int dS/S=(r-q)T$ y $\mathbb{E}[\ln S_T/S_0]=(r-q-\lambda k+\lambda\mu_J)T-\tfrac12\mathbb{E}\int v\,dt$, de donde $\tfrac2T\mathbb{E}[\int dS/S-\ln S_T/S_0]=\tfrac1T\mathbb{E}\int v\,dt+2\lambda(k-\mu_J)$.)

**Prueba de consistencia obligatoria:** la integral de replicación de §3.3 evaluada con precios de opciones del modelo (Fourier) debe reproducir $K^{\text{repl}}$ analítico. Es un test unitario de todo el pipeline (verificado en §12 con error relativo $<10^{-6}$).

**Descomposición para RQ2:** dado que $K^{\text{repl}}_B-K^{\text{repl}}_M$ depende de que $\sigma^2$ sea reemplazado por el término de volatilidad estocástica, el experimento anidado recomendado es

$$
\text{BS}\ \to\ \text{Merton (+saltos)}\ \to\ \text{Heston (+SV)}\ \to\ \text{Bates (+SV+saltos)},
$$

reportando, para cada evento, la contribución de cada componente al strike y su cambio marginal.

### 5.5 Qué aporta cada componente en near-expiry

* **Saltos:** generan skew/curvatura que *no desaparece* cuando $T\to0$ (el smile en wings se hace más pronunciado con vencimientos cortos). Son el motor del smile en 5–10 días.
* **Volatilidad estocástica con $\rho<0$:** aporta skew de nivel casi constante en $T$ y, sobre todo, la **dinámica** de la varianza (mean-reversion) que gobierna la estructura a plazo. En un solo vencimiento corto, $(\kappa,\theta,v_0)$ están débilmente identificados (§8).
* **Limitación conocida:** Bates homogéneo produce curvas de varianza esperada *monótonas* en $T$ (regresan suavemente a $\theta$), por lo que no reproduce un escalón por evento (H4).

---

## 6. Modelar el evento de earnings y el vol crush

### 6.1 Descomposición de varianza total

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

### 6.2 Modelo de evento

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

**Consecuencia para $K_{var}$:** la contribución del evento al strike es $\varepsilon/T$ en QV o $2\Lambda_e(k_e-\mu_e)/T$ en replicación, con el mismo sesgo de tercer orden de §3.5.

### 6.3 Parámetros dependientes del tiempo

Para permitir un escalón en la estructura a plazo, se usa intensidad y variancia **por tramos**:

$$
\lambda(t)=\lambda_{pre}\,\mathbf{1}_{t<t_e^-}+\lambda_e\,\mathbf{1}_{|t-t_e|\le\epsilon}+\lambda_{post}\,\mathbf{1}_{t>t_e^+},\qquad
\theta(t),\ \kappa(t)\ \text{análogos}.
$$

Con parámetros constantes por tramos, el sistema de Riccati (§5.3) se **propaga tramo a tramo** en $\tau=T-t$ con la condición terminal de cada tramo (fórmula cerrada con $D(0)\ne0$ o integración RK4). La integral de intensidad $\Lambda(T)=\int_0^T\lambda(t)\,dt$ sustituye a $\lambda T$ en el término de salto (si $\mu_J,\sigma_J$ son constantes).

### 6.4 Vol crush: definición y lectura en el modelo

Para la opción del mismo vencimiento antes y después del anuncio, y con volatilidad base constante:

$$
\sigma^2_{pre}=\sigma_b^2+\frac{\varepsilon}{T},\qquad \sigma^2_{post}\approx\sigma_b^2,\qquad
CR=\frac{\sigma_{post}}{\sigma_{pre}}\approx\sqrt{\frac{\sigma_b^2}{\sigma_b^2+\varepsilon/T}} .
$$

Mecanismos en los modelos:

1. **Eliminación del salto de evento** ($\lambda(t)$ cae a $\lambda_{post}$): es el mecanismo dominante y natural en Merton/Bates por tramos.
2. **Reversión de $v_t$** hacia $\theta$ (Bates homogéneo): solo produce una caída gradual, no instantánea.
3. **Reset de varianza en $t_e$:** $v_{t_e^+}=c\cdot v_{t_e^-}$ con $c\in(0,1]$, o $\theta(t)$ por tramos. Modelos con saltos *positivos* en varianza (SVJJ; Duffie et al., 2000; Eraker, Johannes y Polson, 2003) generan alzas de varianza, no colapsos, y no son adecuados para el crush sin modificaciones.

**Impacto sobre el swap.** Por §3.1, la posición larga comprada antes del EA gana $R_e^2$ (aporta a $RV_{0,t}$) y pierde la diferencia entre $K_{t,T}$ pre y post: el vol crush ocurre en el término forward $K_{t,T}$, mientras el salto realizado ocurre en $RV_{0,t}$. La **participación del evento** $\eta=\frac{R_e^2}{\sum_iR_i^2}$ es un estadístico clave de cada trayectoria simulada.

---

## 7. Valoración de opciones por transformadas de Fourier

### 7.1 Gil-Pelaez / Heston (1993)

$$
C(K)=S_0e^{-qT}P_1-Ke^{-rT}P_2,
$$

$$
P_2=\frac12+\frac1\pi\int_0^\infty\operatorname{Re}\Big[\frac{e^{-iu\ln K}\varphi(u)}{iu}\Big]du,\qquad
P_1=\frac12+\frac1\pi\int_0^\infty\operatorname{Re}\Big[\frac{e^{-iu\ln K}\varphi(u-i)}{iu\,\varphi(-i)}\Big]du .
$$

### 7.2 Carr–Madan (1999): FFT con amortiguamiento

Con factor $\alpha>0$ de amortiguamiento y $k=\ln K$:

$$
C(K)=\frac{e^{-\alpha k}}{\pi}\int_0^\infty e^{-i\nu k}\,\psi(\nu)\,d\nu,\qquad
\psi(\nu)=\frac{e^{-rT}\,\varphi\big(\nu-(\alpha+1)i\big)}{\alpha^2+\alpha-\nu^2+i(2\alpha+1)\nu}.
$$

Se evalúa con FFT para toda una malla de strikes por vencimiento. Existe además la formulación de Lewis (2001) y el análisis de error de Lee (2004).

### 7.3 Aspectos numéricos para near-expiry

* Para $T\approx1/50$ años el decaimiento de $\varphi(u)$ en $u$ puede ser lento (más aún en Bates por $\xi$); use límites de integración amplios y **verifique convergencia** doblando el límite superior.
* Para una calibración rápida, **cuadratura de Gauss–Legendre fija** sobre la integral de Lewis/Gil-Pelaez, común a todos los strikes de un vencimiento (Lee, 2004; Cui, del Baño Rollin y Germano, 2017, para Heston).
* Precios de puts OTM muy pequeños: evitar cancelación numérica al usar paridad put-call en strikes muy bajos; integrar $P(K)$ vía Fourier directamente o usar Black–Scholes sobre la vol interpolada.

---

## 8. Calibración

### 8.1 Datos (Paso 1)

Por evento y fecha de observación:

* Cadena completa de calls y puts a 2–3 vencimientos: uno **pre-EA** (si existe), los **near-expiry que cruzan** el anuncio (5, 7, 10 días) y uno posterior (para anclar $\kappa,\theta$).
* Bid/ask, open interest, volumen; spot y ex-dividendos dentro de la ventana; tasa libre de riesgo y costo de préstamo.
* Fecha/hora exacta del anuncio (*before market open* vs. *after market close*), pues determina **qué retorno diario** contiene el salto.

### 8.2 Limpieza y preparación

1. **Forward implícito** por paridad put-call en strikes ATM: $F_0=K+e^{rT}(C-P)$, promediado en varios strikes; deduce $q$ efectivo (dividendos + *borrow*).
2. **Opciones americanas.** Las opciones sobre acciones individuales son americanas. Convertirlas a equivalentes europeas (binomial/Bjerksund–Stensland con dividendos discretos) o usar solo **OTM** (donde la prima de ejercicio anticipado es pequeña) para inferir la vol implícita europea.
3. **Filtros:** descartar bid$=0$, spreads relativos excesivos, precios que violen límites de no-arbitraje, y strikes fuera de un rango de delta razonable.
4. **No arbitraje:** comprobar convexidad en $K$ (mariposas) y monotonía calendario; corregir con SVI/Fengler antes de calibrar.
5. Trabajar en **tiempo de negociación** (§11).

### 8.3 Función objetivo

Sean $\Theta$ los parámetros y $i$ el índice de opción:

$$
\hat\Theta=\arg\min_{\Theta}\ \sum_{i=1}^{M}w_i\Big(\sigma^{model}_i(\Theta)-\sigma^{mkt}_i\Big)^2+\alpha\,\big\|\Theta-\Theta_{prior}\big\|^2 .
$$

* Equivalente en precios: $\sum_iw_i'\big(C_i^{model}-C_i^{mkt}\big)^2$ con $w_i'=1/\text{vega}_i^2$ (o $1/(\text{ask}-\text{bid})^2$), porque $\Delta C\approx\text{vega}\cdot\Delta\sigma$. **El RMSE de precios sin ponderar** sobreajusta las opciones caras (ITM/largas) y falla en wings baratos, que son los que determinan el sesgo por saltos.
* $w_i$: inverso del spread bid–ask o de la varianza del error de medición.
* El término $\alpha\|\Theta-\Theta_{prior}\|^2$ (regularización de Tikhonov) estabiliza parámetros débilmente identificados ($\kappa,\theta$).

### 8.4 Restricciones y optimizador

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

### 8.5 Identificabilidad

* Con **un solo vencimiento**, $(v_0,\theta,\kappa)$ son casi colineales: en $T\to0$ el modelo solo depende de $v_0$ y de la combinación $(\rho,\xi)$; $\kappa,\theta$ influyen apenas. Estrategias: calibrar **conjuntamente** varios vencimientos, fijar $\kappa$ y $\theta$ con vencimientos largos, o regularizar hacia priors históricos.
* Bates tiene 8 parámetros y Merton 4: comparar con criterios penalizados (AIC/BIC) y **validación fuera de muestra** (calibrar en strikes impares y testear en pares; o calibrar en $t-1$ y evaluar en $t$).
* Los saltos y la vol estocástica compiten para explicar el mismo smile: monitorear la correlación entre estimadores.

---

## 9. Simulación Monte Carlo

### 9.1 Merton: esquema exacto

Para paso $\Delta$ (típicamente 1 día = $1/252$):

$$
\ln S_{t+\Delta}=\ln S_t+\Big(r-q-\tfrac12\sigma^2-\lambda k\Big)\Delta+\sigma\sqrt\Delta\,Z_1+N\mu_J+\sigma_J\sqrt N\,Z_2,\qquad N\sim\text{Poisson}(\lambda\Delta),
$$

con $Z_1,Z_2\sim\mathcal{N}(0,1)$ independientes. La suma de $N$ saltos normales es $\mathcal{N}(N\mu_J,N\sigma_J^2)$, por lo que **no hay error de discretización** (el esquema es exacto en distribución).

### 9.2 Bates: esquema QE de Andersen (2008)

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

### 9.3 Ventana de earnings

Se simula sobre la malla diaria hasta el vencimiento ($n$ pasos) con parámetros por tramos:

1. **Pre-EA:** $\lambda_{pre}$ (bajo) y $v_0$ calibrado.
2. **Día del EA** (índice de la fila que contiene el salto según BMO/AMC): $\Lambda_e=\lambda_e\Delta_e\approx1$ (p. ej. $\lambda_e=252$ en un paso diario) **o** salto único con distribución elegida en §6.2.
3. **Post-EA:** $\lambda_{post}$ y, opcionalmente, reset de varianza $v\leftarrow c\,v$ (vol crush).

El compensador del salto en el paso de evento es $-\Lambda_ek_e$ (no $-\lambda k\Delta$ con la intensidad base).

### 9.4 Estadísticos por trayectoria y estimadores

* $RV^{(m)}=\frac{252}{n}\sum_{i=1}^{n}\big(\ln S^{(m)}_{t_i}/S^{(m)}_{t_{i-1}}\big)^2$ (misma convención que el contrato).
* $\hat K=\frac1M\sum_mRV^{(m)}$ con error estándar $\hat s/\sqrt M$. Contrastar $\hat K$ con: $K^{disc}$ (Merton, exacto), $K^{QV}$ y $K^{\text{repl}}$ analíticos.
* $P\&L^{(m)}=N_{var}\big(RV^{(m)}-K_{var}\big)$: media, desviación, asimetría, cuantiles al 1/5%.
* $\eta^{(m)}=R_{e}^2/\sum_iR_i^2$ (participación del día del evento) y $RV^{(m)}_{ex}$ (RV excluyendo el día del EA).
* Chequeo de **martingala**: $\mathbb{E}[S_T]e^{-(r-q)T}\approx1$ (en el ejemplo de validación, error $\sim10^{-5}$ con $4\times10^5$ trayectorias).

### 9.5 Reducción de varianza y semillas

* Variables antitéticas en $Z$; *control variate* usando una variable con esperanza conocida (por ejemplo, $\int v_t\,dt$ aproximada con la regla trapezoidal, cuya esperanza es la de §5.1).
* $M\ge10^5$–$10^6$ trayectorias; fijar semilla (`numpy.random.default_rng(seed)`) y reportar el intervalo de confianza.

---

## 10. Diseño empírico y métricas

### 10.1 Cantidades de comparación

| Cantidad | Definición |
|----------|-----------|
| $K^{mkt}_{repl}$ | Strike model-free del mercado (§3.3–3.4, tira de opciones con extrapolación) |
| $K^{model}_{repl}$ | Integral de replicación con precios del modelo calibrado; coincide con la fórmula analítica de §4.5 / §5.4 |
| $K^{model}_{QV}$ | Variación cuadrática esperada bajo el modelo |
| $\widehat{K}^{MC}$ | Estimador Monte Carlo de $\mathbb{E}[RV]$ discreta |
| $RV^{\mathbb{P}}$ | Varianza realizada observada tras el evento (RV histórica "sintética" por ticker/ventana) |

### 10.2 Error de pricing (RQ1)

* **Ajuste del smile:** $\text{RMSE}_{IV}=\sqrt{\frac1M\sum_i(\sigma_i^{model}-\sigma_i^{mkt})^2}$; también ponderado por vega.
* **Error en el strike:** $\epsilon_K=\dfrac{K^{model}_{repl}-K^{mkt}_{repl}}{K^{mkt}_{repl}}$ y en puntos de vol, $\sqrt{K^{model}}-\sqrt{K^{mkt}}$.
* **Comparación penalizada:** $\text{BIC}=M\ln(\text{RSS}/M)+p\ln M$ con $p=4$ (Merton) y $p=8$ (Bates).
* **Significancia:** contraste de **Diebold–Mariano** (1995) sobre diferencias de pérdidas absolutas/cuadráticas entre modelos, con errores estándar **clusterizados por fecha** (los EA se concentran en la misma semana de resultados) o *block bootstrap*.
* Desagregar por moneyness, vencimiento y tamaño esperado del movimiento.

### 10.3 Sesgo por saltos (RQ3)

Reportar, por evento: $K^{\text{repl}}-K^{QV}$ en absoluto y en % de $K^{QV}$, junto con $\hat\lambda\hat{\mathbb{E}}[J^3]$ del modelo calibrado, y contrastar con la relación teórica de §3.5.

### 10.4 Vol crush y estructura temporal (RQ4)

* **Crush observado:** $CR^{mkt}=\sigma^{ATM}_{post}/\sigma^{ATM}_{pre}$ del mismo vencimiento (o del vencimiento constante interpolado).
* **Crush del modelo:** recalcular $\sigma^{ATM}$ tras un paso de evento con parámetros post (§6.4) y comparar contra $CR^{mkt}$.
* **Varianza forward:** $K_{T_1,T_2}$ de mercado (§3.7) vs. la del modelo, antes y después del anuncio. Bajo H4, Bates homogéneo no reproduce el escalón; Bates por tramos sí, a costa de parámetros adicionales.
* **Varianza de evento:** $\hat\varepsilon^{mkt}$ (§6.1) vs. $\hat\varepsilon^{model}=\Lambda_e(\mu_e^2+\sigma_e^2)$ y vs. $\mathbb{E}^{\mathbb{P}}[R_e^2]$ histórico.

### 10.5 Prima de riesgo de varianza y predictibilidad (RQ5)

$$
VRP_i=K_{var,i}-RV_i,\qquad RV_i=\alpha+\beta\,K_{var,i}+u_i\ \ (\text{Mincer–Zarnowitz}),\quad H_0:(\alpha,\beta)=(0,1).
$$

Errores de Newey–West (1987) o clusterizados. La literatura muestra prima de riesgo de varianza en general (Carr y Wu, 2009; Bollerslev, Tauchen y Zhou, 2009) y, en EA, primas heterogéneas ligadas al riesgo no diversificable (Barth y So, 2014) y rendimientos de straddles previos a earnings (Gao, Xing y Zhang, 2018).

---

## 11. Checklist de consideraciones prácticas

**Datos y convenciones**

- [ ] **Tiempo de negociación vs. calendario.** Usar $\tau=(\text{días hábiles})/252$. En calendario, fines de semana diluyen la varianza y distorsionan $\sigma_b$ y $\varepsilon$.
- [ ] **BMO/AMC.** El salto de un anuncio *after close* aparece en el retorno del día siguiente; alinear el índice del salto con $RV$ close-to-close.
- [ ] **Dividendos y ex-date** en la ventana (ajustan $F_0$ y el spot ex-div).
- [ ] **Opciones americanas** (§8.2).
- [ ] **Liquidez de wings** en semanales; *stale quotes*.
- [ ] **Look-ahead bias:** calibrar solo con cotizaciones anteriores al anuncio.
- [ ] **Múltiples contrastes:** ajustar (Bonferroni/BH) al comparar muchos tickers y ventanas.

**Modelo y numérica**

- [ ] Condición de Feller y positividad de $v_t$ (QE).
- [ ] Martingala: $\varphi(-i)=F_0$.
- [ ] Integrales de Fourier convergidas (doblar $u_{max}$).
- [ ] Test unitario: integral de replicación vs. fórmula cerrada de §4.5/§5.4.
- [ ] Test unitario: Merton por serie vs. Fourier.
- [ ] Test unitario: MC vs. $K^{disc}$ (Merton).
- [ ] Rango de strikes y extrapolación: reportar sensibilidad de $K^{mkt}_{repl}$.
- [ ] Multiplicidad de mínimos y estabilidad de parámetros entre días.

**Interpretación**

- [ ] No mezclar $K^{\text{repl}}$ con $K^{QV}$ (§3.5).
- [ ] Los parámetros calibrados son bajo $\mathbb{Q}$; no interpretarlos como frecuencias históricas de saltos.
- [ ] Modelos de un solo salto gaussiano no explican concavidad/bimodalidad (§6.2 (v)).

---

## 12. Ejemplo numérico sintético (reproducible)

Parámetros **ilustrativos** (no calibrados a datos reales): $S_0=100$, $r=4\%$, $q=0$, $T=7/252$ (7 días hábiles).

| | Merton | Bates |
|---|---|---|
| Difusión / SV | $\sigma=0.30$ | $v_0=\theta=0.09$, $\kappa=3$, $\xi=1.0$, $\rho=-0.6$ |
| Saltos | $\lambda=36$, $\mu_J=-0.005$, $\sigma_J=0.07$ | $\lambda=36$, $\mu_J=-0.01$, $\sigma_J=0.06$ |
| $\lambda T$ (saltos esperados) | 1.0 | 1.0 |
| $K^{\text{repl}}$ analítico | 0.26663 (vol 51.6%) | 0.22201 (vol 47.1%) |
| $K^{\text{repl}}$ por integral de Fourier | 0.26663 | 0.22201 |
| $K^{QV}$ | 0.26730 (vol 51.7%) | 0.22320 (vol 47.2%) |
| Aporte del salto a $K^{\text{repl}}$ | 0.1766 (66%) | 0.1320 (59%) |
| IV en $K$=90 / 95 / 100 / 105 / 110 | 53.6 / 50.4 / 48.8 / 49.6 / 51.9 % | 51.2 / 47.9 / 45.1 / 43.9 / 45.2 % |

Lecturas:

* El strike de replicación **coincide con el analítico a ~6 dígitos**, lo que valida el pipeline de Fourier (§7) y las fórmulas (§4.5, §5.4).
* Con $\mu_J<0$, $K^{\text{repl}}<K^{QV}$ como predice §3.5, y la brecha es de orden $\lambda\mathbb{E}[J^3]/3$.
* Con un solo salto esperado en la ventana, **la mayor parte (59–66%) del strike proviene de riesgo de salto**, coherente con H3.
* Validación Monte Carlo (paso diario, $4\times10^5$ trayectorias, otros parámetros de prueba): $\widehat K^{MC}_{Merton}=0.1401\pm0.0003$ vs. $K^{disc}=0.1405$; $\widehat K^{MC}_{Bates}=0.1374\pm0.0003$ vs. $K^{QV}_B=0.1376$; martingala: $1.00002$.

---

## 13. Referencia de implementación en Python

Snippets validados numéricamente (Fourier, replicación y simulación). Dependencias: `numpy`, `scipy`.

### 13.1 Funciones características

```python
import numpy as np
from scipy.integrate import quad

def cf_merton(u, T, S0, r, q, sigma, lam, muJ, sJ):
    k = np.exp(muJ + 0.5*sJ**2) - 1
    drift = np.log(S0) + (r - q - 0.5*sigma**2 - lam*k)*T
    return np.exp(1j*u*drift - 0.5*sigma**2*u**2*T
                  + lam*T*(np.exp(1j*u*muJ - 0.5*sJ**2*u**2) - 1))

def cf_bates(u, T, S0, r, q, v0, kappa, theta, xi, rho, lam, muJ, sJ):
    k = np.exp(muJ + 0.5*sJ**2) - 1
    b = kappa - rho*xi*1j*u
    d = np.sqrt(b**2 + xi**2*(1j*u + u**2))
    g = (b - d)/(b + d)                       # "Little Heston Trap"
    e = np.exp(-d*T)
    D = (b - d)/xi**2 * (1 - e)/(1 - g*e)
    C = kappa*theta/xi**2*((b - d)*T - 2*np.log((1 - g*e)/(1 - g)))
    jump = lam*T*(np.exp(1j*u*muJ - 0.5*sJ**2*u**2) - 1 - 1j*u*k)
    return np.exp(1j*u*(np.log(S0) + (r - q)*T) + C + D*v0 + jump)
```

### 13.2 Precio de call (Gil-Pelaez) y strike por replicación

```python
def call_price(cf, S0, K, T, r, q, upper=400.0):
    """cf: callable u -> varphi(u).  Se evalua tambien en u-1j."""
    k = np.log(K)
    phi_mi = cf(-1j)                          # = E[S_T] = F0
    f2 = lambda u: np.real(np.exp(-1j*u*k)*cf(u)/(1j*u))
    f1 = lambda u: np.real(np.exp(-1j*u*k)*cf(u - 1j)/(1j*u*phi_mi))
    P2 = 0.5 + quad(f2, 1e-10, upper, limit=500)[0]/np.pi
    P1 = 0.5 + quad(f1, 1e-10, upper, limit=500)[0]/np.pi
    return S0*np.exp(-q*T)*P1 - K*np.exp(-r*T)*P2

def kvar_replication(cf, S0, T, r, q, kmin=0.05, kmax=8.0):
    """K_var = (2/T) e^{rT} [ int_0^F P/K^2 + int_F^inf C/K^2 ]"""
    F = S0*np.exp((r - q)*T)
    call = lambda K: call_price(cf, S0, K, T, r, q)
    put  = lambda K: call(K) - S0*np.exp(-q*T) + K*np.exp(-r*T)   # paridad
    Ip = quad(lambda K: put(K)/K**2,  kmin*S0, F,        limit=200)[0]
    Ic = quad(lambda K: call(K)/K**2, F,       kmax*S0,  limit=200)[0]
    return 2/T*np.exp(r*T)*(Ip + Ic)
```

### 13.3 Formas cerradas para validar

```python
def kvar_merton_closed(sigma, lam, muJ, sJ):
    k = np.exp(muJ + 0.5*sJ**2) - 1
    return dict(repl=sigma**2 + 2*lam*(k - muJ),
                qv=sigma**2 + lam*(muJ**2 + sJ**2))

def kvar_bates_closed(v0, kappa, theta, lam, muJ, sJ, T):
    k = np.exp(muJ + 0.5*sJ**2) - 1
    sv = theta + (v0 - theta)*(1 - np.exp(-kappa*T))/(kappa*T)
    return dict(repl=sv + 2*lam*(k - muJ),
                qv=sv + lam*(muJ**2 + sJ**2))
```

### 13.4 Simulación exacta de Merton con salto de evento

```python
def simulate_merton_earnings(M, n, t_e, S0, r, q, sigma, lam_base, muJ, sJ,
                             Lambda_e, mu_e, s_e, seed=0):
    """t_e: indice (0..n-1) del retorno diario que contiene el salto de EA.
    En ese paso el salto base se reemplaza por el salto de evento (Lambda_e ~ 1)."""
    rng = np.random.default_rng(seed)
    dt = 1/252
    k  = np.exp(muJ + 0.5*sJ**2) - 1
    ke = np.exp(mu_e + 0.5*s_e**2) - 1
    R = np.empty((M, n))
    for i in range(n):
        mu_i, s_i, k_i, Lam = muJ, sJ, k, lam_base*dt
        if i == t_e:
            mu_i, s_i, k_i, Lam = mu_e, s_e, ke, Lambda_e
        N  = rng.poisson(Lam, M)
        Z1 = rng.standard_normal(M); Z2 = rng.standard_normal(M)
        R[:, i] = ((r - q - 0.5*sigma**2)*dt - Lam*k_i
                   + sigma*np.sqrt(dt)*Z1 + N*mu_i + s_i*np.sqrt(N)*Z2)
    RV  = 252*np.mean(R**2, axis=1)
    eta = R[:, t_e]**2/np.sum(R**2, axis=1)      # participacion del evento
    return RV, eta
```

### 13.5 Esqueleto de calibración

```python
from scipy.optimize import differential_evolution, least_squares

def residuals(theta, model, mkt):
    """mkt: arrays T, K, iv_mkt, w (pesos). Devuelve sqrt(w)*(iv_model-iv_mkt)."""
    iv_model = model_implied_vols(theta, mkt)      # via Fourier + inversion BS
    return np.sqrt(mkt["w"])*(iv_model - mkt["iv"])

# 1) Busqueda global (acotada)      2) Refinamiento local
# res_g = differential_evolution(lambda th: np.sum(residuals(th, model, mkt)**2),
#                                bounds=bounds, seed=1, maxiter=200, tol=1e-8)
# res_l = least_squares(residuals, res_g.x, args=(model, mkt), bounds=(lb, ub))
```

---

## 14. Estructura sugerida del repositorio

```
.
├── README.md                     # Este documento (marco teórico)
├── data/
│   ├── raw/                      # Cadenas de opciones y earnings dates
│   └── processed/                # Cadenas limpias, forwards, vols interpoladas
├── src/
│   ├── models/
│   │   ├── merton.py             # cf, serie de Merton, K_var cerrado
│   │   ├── heston_bates.py       # cf (Little Trap), K_var cerrado
│   │   └── event.py              # factor de evento, parámetros por tramos
│   ├── pricing/
│   │   ├── fourier.py            # Gil-Pelaez, Carr-Madan, Lewis/Gauss-Legendre
│   │   └── variance_swap.py      # replicación, discretización, extrapolación
│   ├── calibration/
│   │   ├── prep.py               # forward por paridad, americanas, filtros
│   │   └── fit.py                # DE + least_squares, regularización
│   ├── simulation/
│   │   ├── merton_mc.py
│   │   └── bates_qe.py
│   └── analysis/
│       ├── metrics.py            # RMSE_IV, BIC, Diebold-Mariano, Mincer-Zarnowitz
│       └── volcrush.py           # CR, varianza forward, varianza de evento
├── tests/                        # ver checklist §11 (tests unitarios)
├── notebooks/
└── requirements.txt              # numpy, scipy, pandas, matplotlib
```

---

## 15. Limitaciones y extensiones

**Limitaciones del diseño**

* **Modelos afines de una sola componente de salto gaussiana** no reproducen bimodalidad ni concavidad del smile pre-EA; Merton/Bates son *benchmarks* estándar, no el estado del arte para eventos.
* **Identificación** débil de $(\kappa,\theta,v_0)$ con un único vencimiento corto.
* **Prima de riesgo de salto** absorbida en los parámetros $\mathbb{Q}$: no se puede separar de la intensidad histórica sin datos $\mathbb{P}$.
* **Opciones americanas**: el proceso de conversión introduce ruido de modelo.
* **Independencia** entre salto de evento y proceso base; en la práctica el tamaño del salto puede correlacionarse con el nivel de $v_t$.

**Extensiones naturales**

* Saltos con **ley asimétrica o bimodal** (Kou; mezcla de dos gaussianas), y **SVCJ/SVJJ** (Duffie et al., 2000; Eraker et al., 2003) con saltos en $v$ correlacionados.
* **Superficie cross-section** de EA: regresiones de $\varepsilon$ en características de la empresa (Dubinsky et al., 2019; Barth y So, 2014).
* Comparar contra **modelos de varianza forward / rough volatility** (Bergomi, 2016).
* Contratos relacionados: **corridor variance swaps** y **gamma swaps** (Carr y Lewis, 2004) para aislar el tramo del salto.

---

## 16. Bibliografía

> **Nota:** las notas de research de bancos (Goldman Sachs, J.P. Morgan) son documentos de *sell-side*; su disponibilidad en línea varía. Verifique DOI/páginas y accesibilidad antes de someter el paper.

### 16.1 Research de bancos (variance swaps y volatilidad)

* Demeterfi, K., Derman, E., Kamal, M. y Zou, J. (1999). *More Than You Ever Wanted to Know About Volatility Swaps*. Goldman Sachs Quantitative Strategies Research Notes, marzo de 1999. Versión publicada: *A Guide to Volatility and Variance Swaps*, **Journal of Derivatives**, 6(4), 9–32.
* Derman, E., Kamal, M., Kani, I., McClure, J., Pirasteh, C. y Zou, J. (1998). *Investing in Volatility*. Futures & Options World. (Goldman Sachs).
* Bossu, S., Strasser, E. y Guichard, R. (2005). *Just What You Need to Know About Variance Swaps*. J.P. Morgan, Equity Derivatives Investor Marketing / Quantitative Research and Development, Londres, mayo de 2005.
* Allen, P., Einchcomb, S. y Granger, N. (2006). *Variance Swaps*. J.P. Morgan Securities Ltd., European Equity Derivatives Research.
* Bergomi, L. (2016). *Stochastic Volatility Modeling*. Chapman & Hall/CRC. (Práctica cuantitativa de Société Générale).

### 16.2 Modelos y replicación

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

### 16.3 Variance/volatility swaps con saltos y discretización

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

### 16.4 Earnings announcements y opciones

* Dubinsky, A., Johannes, M., Kaeck, A. y Seeger, N. J. (2019). Option Pricing of Earnings Announcement Risks. **Review of Financial Studies**, 32(2), 646–687.
* Barth, M. E. y So, E. C. (2014). Non-Diversifiable Volatility Risk and Risk Premiums at Earnings Announcements. **The Accounting Review**, 89(5), 1579–1607.
* Gao, C., Xing, Y. y Zhang, X. (2018). Anticipating Uncertainty: Straddles Around Earnings Announcements. **Journal of Financial and Quantitative Analysis**, 53(6), 2587–2617.
* Alexiou, L., Goyal, A., Kostakis, A. y Rompolis, L. (2025). Pricing Event Risk: Evidence from Concave Implied Volatility Curves. **Review of Finance**, 29(4), 963–1007.
* Kachhara, D., Markin, J. K. E. y Singh, A. (2023). *Option Smile Volatility and Implied Probabilities: Implications of Concavity in IV Curves*. arXiv:2307.15718.
* Patell, J. y Wolfson, M. (1979). Anticipated Information Releases Reflected in Call Option Prices. **Journal of Accounting and Economics**, 1(2), 117–140.
* Beaver, W. (1968). The Information Content of Annual Earnings Announcements. **Journal of Accounting Research**, 6 (Supplement), 67–92.
* Isakov, D. y Périgon, C. (2001). Evolution of Market Uncertainty Around Earnings Announcements. **Journal of Banking & Finance**, 25(9), 1769–1788.

### 16.5 Métodos numéricos, simulación y calibración

* Carr, P. y Madan, D. (1999). Option Valuation Using the Fast Fourier Transform. **Journal of Computational Finance**, 2(4), 61–73.
* Lewis, A. L. (2001). *A Simple Option Formula for General Jump-Diffusion and Other Exponential Lévy Processes*. SSRN.
* Lee, R. W. (2004). Option Pricing by Transform Methods: Extensions, Unification and Error Control. **Journal of Computational Finance**, 7(3), 51–86.
* Albrecher, H., Mayer, P., Schoutens, W. y Tistaert, J. (2007). The Little Heston Trap. **Wilmott Magazine**, enero de 2007, 83–92.
* Andersen, L. (2008). Simple and Efficient Simulation of the Heston Stochastic Volatility Model. **Journal of Computational Finance**, 11(3), 1–42.
* Lord, R., Koekkoek, R. y van Dijk, D. (2010). A Comparison of Biased Simulation Schemes for Stochastic Volatility Models. **Quantitative Finance**, 10(2), 177–194.
* Cui, Y., del Baño Rollin, S. y Germano, G. (2017). Full and Fast Calibration of the Heston Stochastic Volatility Model. **European Journal of Operational Research**, 263(2), 625–638.
* Gatheral, J. y Jacquier, A. (2014). Arbitrage-Free SVI Volatility Surfaces. **Quantitative Finance**, 14(1), 59–71.
* Fengler, M. R. (2009). Arbitrage-Free Smoothing of the Implied Volatility Surface. **Quantitative Finance**, 9(4), 417–428.

### 16.6 Inferencia estadística

* Diebold, F. X. y Mariano, R. S. (1995). Comparing Predictive Accuracy. **Journal of Business & Economic Statistics**, 13(3), 253–263.
* Newey, W. K. y West, K. D. (1987). A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix. **Econometrica**, 55(3), 703–708.
* Mincer, J. y Zarnowitz, V. (1969). The Evaluation of Economic Forecasts. En *Economic Forecasts and Expectations*, NBER.

---

**Licencia y uso:** documento de investigación con fines académicos; no constituye asesoramiento financiero.
