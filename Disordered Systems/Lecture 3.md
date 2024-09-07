16.04.2024

## 2.3 Curie-Weiss model

### 2.3.1 Setup: 

$\sigma_i = \pm 1\, ;\, i = \{1, \dots ,N\}$
$H(\mathbf{\sigma}) = -\frac{1}{2}\sum_{i, j}J_{ij}\sigma_i\sigma_j - \sum_i h_i\sigma_i$
where $J_{ij} = \frac{J_0}{N}, \; i \neq j \;\& \; h_i=h$

Also, $m(\mathbf{\sigma} = \frac{1}{N}\Sigma_i\sigma_i)$

$$
\therefore H(\mathbf{\sigma}) = -\frac{NJ_0}{2}\left[m^2(\mathbf{\sigma}) - Nhm(\mathbf{\sigma})\right] + \mathcal{O}(1)
$$

Now, in
$$
Z = \sum_{\mathbf{\sigma}}e^{-\beta H(\mathbf{\sigma})} = \sum_{\mathbf{\sigma}}e^{\frac{N\beta J_0}{2}m^2(\mathbf{\sigma}) + N\beta hm(\mathbf{\sigma})}
$$
use $$
e^{\frac{N\beta J_0}{2} m^2(\sigma)} = \int\frac{dz}{\sqrt{2\pi \beta\frac{J_0}{N}}}e^{-\frac{1}{2}\frac{z^2}{\beta \frac{J_0}{N}} + Nzm(\mathbf{\sigma})}
$$

to get
$$
\begin{align} \\
Z &= \sum_{\mathbf{\sigma}}\int \frac{dz}{\sqrt{\frac{2\pi\beta J_0}{N}}}e^{-\frac{N}{2}\frac{z^2}{\beta J_0}}e^{N(\beta h + z)m(\mathbf{\sigma})} \\
&= \int \frac{dz}{\sqrt{\frac{2\pi\beta J_0}{N}}}e^{-\frac{N}{2}\frac{z^2}{\beta J_0}}\sum_{\mathbf{\sigma}}e^{N(\beta h + z)m(\mathbf{\sigma})} \\
&=\int \frac{dz}{\sqrt{\frac{2\pi\beta J_0}{N}}}e^{-\frac{N}{2}\frac{z^2}{\beta J_0}}\cdot\left[2\cosh(\beta h + z)\right]^N \\
&= \int \frac{dz}{\sqrt{\frac{2\pi\beta J_0}{N}}}e^{-N\beta f(z)}
\end{align}
$$
where $f(z) = \frac{z^2}{2\beta^2 J_0}  - T\ln\left[2\cosh(z+\beta h)\right]$

Expanding the function $f(z)$ around its minimum:
$f(z) := f(z^*) + \frac{1}{2}(z-z^*)^2f^{\prime\prime}(z^*) + \dots$ as $f^\prime(z) = 0$

$$
\begin{align}
\therefore \int e^{-N\beta f(z)} dz &\approx e^{-N\beta f(z^*)}\sqrt{\frac{2\pi}{N\beta f^{\prime\prime}(z^*)}}\cdot e^{-\frac{1}{2}\frac{\left(z-z^*\right)^2}{\sigma^2}} \\
&= e^{-N\beta f(z^*)}\sqrt{\frac{2\pi}{\sigma^2}}\cdot e^{-\frac{1}{2}\frac{\left(z-z^*\right)^2}{\sigma^2}}
\end{align}
$$
where $\sigma^2 \ = frac{1}{N\beta f^{\prime\prime}(z^*)}$

Hence, the partition function becomes 
$$
Z = \frac{1}{\sqrt{2\pi \frac{\beta J_0}{N}}}\sqrt{\frac{2\pi}{N\beta f^\prime(z)}}e^{-N\beta f(z^*)}
$$

and the *Free energy per spin* becomes
$$
f = \frac{F}{N} = - \frac{T}{N}\ln Z = f(z^*) - \frac{T}{2N}\ln\left[\frac{1}{\beta^2 J_0f^{\prime\prime}(z^*)}\right]
$$

**NOTE**: The only thing that matters in the end is how high the peak is and NOT how wide it is!

$\therefore \text{as } N \xrightarrow[]{}\infty:\;f\equiv f(z^*) =\min_z f(z)$
and this happens when 
$$
\begin{align}
f(z^*) &\overset{!}{=} 0 \\
\implies \frac{z}{\beta^2 J_0} &= T\frac{\sinh{(\beta h + z)}}{\cosh (\beta h + z)}\\ \\
\therefore \;z^*&= \beta^2J_0T\tanh(z^* + \beta h)
\end{align}
$$

This is a **self-consistent equation** in $z^*$.

i.e. $$
\begin{align}
f &\equiv \tilde{f}(z^*(h), h) \\
\implies \frac{\partial f}{\partial h}&= \cancelto{0}{\frac{\partial \tilde{f}}{\partial z^*}} \frac{\partial z^*}{\partial h} + \frac{\partial \tilde{f}}{\partial h}  (\text{as } f \text{ attains a minimum @ }z^*)\\
&=\frac{\partial \tilde{f}}{\partial h} \; 
\end{align}
$$

Therefore,
$$
m = -\left.\frac{\partial f}{\partial h}\right|_{z^*} = \tanh(\beta h + z^*)
$$

One can, thus form a relation between $z^*$ and m using $m = \frac{z^*}{\beta J_0}$ and write a **self-consistent equation** for m as:
$$
m = \tanh(\beta J_0 m + \beta h)
$$
this can be interpreted as the energy change by flipping  a single spin having a field $h \xrightarrow[]{} h + J_0 m$. OR the interactions $J_0$ between the spins give an additional field $J_0 m$  acting on each spin resulting in a *mean-field* model of ferromagnetism.

**NOTE**: The temperature in the original happiness problem of  putting people in different rooms represents the strength of how strongly one wants to optimize the problem.

~~< Have to put numerical results for phase transitions here... but until the day arrives: >~~
![[curieWeiss.png]]

Now, we have to find the energy $E = \langle H(\mathbf{\sigma)}\rangle$

$$
\begin{align}
\frac{E}{N} &= -\frac{1}{N}\frac{\partial}{\partial \beta}\ln Z = \frac{\partial}{\partial \beta}(\beta f) \; \left(\because f = - \frac{T}{N}\ln Z\right) \\
&= f + \beta\frac{\partial f}{\partial \beta} \\
&= \frac{z^2}{2\beta^2 J_0} - T \ln (2\cosh (\beta h + z)) + \beta \left(\frac{-2z^2}{2\beta^2 J_0}\right) + \frac{\beta}{\beta^2}\ln 2\cosh (\beta h + z) \\
&= -\frac{z^2}{2\beta^2 J_0} - h \tanh(\beta h+z^*) \\
&= -\frac{J_0}{2}m^2 - hm
\end{align}
$$

This energy per spin has the same form as we started with, but is represented in terms of magnetization instead of summations over individual spins.

For $T\xrightarrow[]{}0: \;m\xrightarrow[]{}sgn(h)$, i.e. every one is in the same room in the party problem. Essentially $\frac{E}{N} = \frac{-J_0}{2} - \left|h\right|$

## 2.4 Random Field Ising Model

### 2.4.1 Setup and Calculation:

Now, all $h_i$'s are different.

$$
\begin{align}
Z &= \sum_{\mathbf{\sigma}}\int \frac{dz}{\sqrt{\frac{2\pi\beta J_0}{N}}}e^{-\frac{N}{2}\frac{z^2}{\beta J_0}}e^{N(\beta h_i + z)m(\mathbf{\sigma})} \\
&= \int \frac{dz}{\sqrt{\frac{2\pi\beta J_0}{N}}}e^{-\frac{N}{2}\frac{z^2}{\beta J_0}}\sum_{\mathbf{\sigma}}e^{N\beta h_i + z)m(\mathbf{\sigma})} \\ 
&= \int \frac{dz}{\sqrt{\frac{2\pi\beta J_0}{N}}}e^{-\frac{N}{2}\frac{z^2}{\beta J_0}}\sum_{\mathbf{\sigma}}e^{\sum_{\sigma{i}}(\beta h_i + z)\sigma_i} \; \left(\because m = \sum_i\sigma_i/N\right) \\
&=\int \frac{dz}{\sqrt{\frac{2\pi\beta J_0}{N}}}e^{-\frac{N}{2}\frac{z^2}{\beta J_0}}\cdot\Pi_{i}\left[2\cosh(\beta h + z)\right]^N \\
&= \int \frac{dz}{\sqrt{\frac{2\pi\beta J_0}{N}}}e^{-N\beta f(z)}
\end{align}
$$
where $f(z) = \frac{z^2}{2\beta^2 J_0}  - T\frac{1}{N}\sum_i\ln\left[2\cosh(z+\beta h)\right]$

And, by a similar calculation as before, the **self-consistent equation** for m is:

$$
\begin{align}
m &= \frac{1}{N}\sum_i\tanh(\beta J_0 m + \beta h_i)\\
&=\int dh \hat{P}(h)\tanh(\beta J_0 m + \beta h)\\ \
\end{align}
$$
where $\hat{P}(h) = \frac{1}{N}\sum_i\delta(h-h_i)$ is the *empirical distribution* of $h_i$ OR a *finely binned histogram* of $h_i$.
If one knows from which distribution $h_i$ have been derived, then for large systems $\hat{P}(h)$ takes the shape of the distribution of $h_i$