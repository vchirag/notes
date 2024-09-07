10.04.2024

# 2. Fully connected disordered systems

## 2.1 A problem class with many flavours

$$
\text{Given: }H(\mathbf{\sigma}) =-\frac{1}{2}\sum_{i,j}J_{ij}\sigma_i\sigma_j -\sum_ih_i\sigma_i
$$

where $\sigma_i = \pm 1$, $J_{ij} = J_{ji};\; J_{ii} = 0$

One can study the following problems with this form of Hamiltonian:
1. **Curie-Weiss mean-field ferromagnet**: $J_{ij} = \frac{J_0}{N};\; h_i = h$
	- It is a long-range, mean-field model. 
1. **Ising model**: $J_{ij} = J_0$ if $i,j$ are neighbours on a lattice
	- It is a short-range model
	- **Random Field** Ising model isused to study *avalanche dynamics* 
1. **Sherrington-Kirkpatric model**: $J_{ij} = J_0 +\text{randomness};\; h_i = h$
2. **Edward-Anderson model**: Like SK but $J_{ij} \neq 0$ only for neighbours on a lattice (the same spatial structure as the *Ising model*)
3. **Viana-Bray model**: Like SK but $J_{ij}\neq0$ only for some *friendship* network
4. **Restricted Boltzmann machine**: hidden spins as generating bases on the visible ones.
5. **Hopfield model**: $J_{ij} = \frac{1}{N}\sum_{\mu=1}^{P}\xi_i^{\mu}\xi_j^{\mu}$, has an *associative memory* of storing $P$ patterns $\xi^{\mu}$

## 2.2 Recap of Statistical Physics

1. Boltzmann distribution: $P(\mathbf{\sigma}) = \frac{1}{Z}\exp\{-\beta H(\mathbf{\sigma})\}$
2. *Partition function* $Z=\sum_{\mathbf{\sigma}}\exp\{-\beta H(\mathbf{\sigma})\}$ where $\beta = \frac{1}{T}$ and $k_B = 1$
3. Average Energy $E = \langle H\rangle = \sum_{\mathbf{\sigma}}P(\mathbf{\sigma})H(\mathbf{\sigma}) = -\frac{\partial}{\partial\beta}\ln Z$
	- Proof: 
	$$\begin{align} \\
-\frac{\partial}{\partial \beta}\ln Z &= -\frac{1}{Z}\frac{\partial Z}{\partial \beta}\\ \\
&= \frac{1}{Z}\sum_{\mathbf{\sigma}}H(\mathbf{\sigma})\exp\{-\beta H(\mathbf{\sigma})\}\\ \\
&= \sum_{\mathbf{\sigma}}H(\mathbf{\sigma})P(\mathbf{\sigma})\\ \\
&= E
	 \end{align}$$
4. Limit $\beta \rightarrow \infty$ or $T \rightarrow 0$ gives the ground state limit.
5. *Entropy* $S = -\sum_{\mathbf{\sigma}}P(\mathbf{\sigma})\ln P(\mathbf{\sigma}) = -\sum_{\mathbf{\sigma}}P(\mathbf{\sigma})\left[-\beta H(\mathbf{\sigma}) - \ln Z\right] = \beta E + \ln Z$ 
	- For uniform $P(\mathbf{\sigma}) = \frac{1}{2^N}$, S measures the spread of the distribution: $0 \;(@\,T\rightarrow 0)\leq S\leq N\ln2\;(@\,T\rightarrow\infty)$
6. *Free Energy* $F = -T\ln Z$. Hence, $S = -\frac{\partial F}{\partial T} = \beta E - \beta F$ 
7. S is related to the number of configurations @ a given energy with the aid of *Density of States* $\Omega(E) = \sum_{\mathbf{\sigma}}\delta(H(\mathbf{\sigma})-E)$
	- The number of states $\mathbf{\sigma}$ in a small interval $E < H(\mathbf{\sigma}) < E + dE$ are then $\Omega(E)dE$
	- For large N, $\Omega(E) \sim e^{S(E)}$ upto some sub-exponential corrections in N

## 2.3 Curie-Weiss Model (Mean-Field Ferromagnet)

### 2.3.1 Setup

Choose all the interactions to be equal: $J_{ij} = \frac{J_0}{N}$ where $i\neq j$ and the biases/preferences t be the same $h_i = h$.
NOTE: $N$ in the denominator makes the interaction finite and $i \neq j$ definition avoids any *self-interaction*.
$$\begin{align}
H(\mathbf{\sigma}) &=-\frac{1}{2}\sum_{i,j}J_{ij}\sigma_i\sigma_j -\sum_ih_i\sigma_i\\
&=-\frac{1}{2}\sum_{i\neq j}\frac{J_0}{N}\sigma_i\sigma_j -h\sum_i\sigma_i\\
&= -\frac{1}{2}\frac{J_0}{N}\sum_{i, j}\sigma_i\sigma_j + \frac{1}{2}\frac{J_0}{N}\sum_{i}\sigma_i\sigma_i -h\sum_i\sigma_i\\ \\
&= -\frac{1}{2}\frac{J_0}{N}\sum_{i, j}\sigma_i\sigma_j + \frac{1}{2}\frac{J_0}{N}\left(1+\cdots+1\right) -h\sum_i\sigma_i\\
&= -\frac{1}{2}\frac{J_0}{N}\sum_{i, j}\sigma_i\sigma_j + \frac{J_0}{2} -h\sum_i\sigma_i
\end{align}
$$
**Define**: *magnetization* $m(\mathbf{\sigma}) := \frac{1}{N}\sum_i\sigma_i$

Then 
$$
\begin{align}
H(\mathbf{\sigma})&= -\frac{1}{2}\frac{J_0}{N}\sum_{i, j}\sigma_i\sigma_j + \frac{J_0}{2} -h\sum_i\sigma_i\\ \\
&= -\frac{J_0}{2N}\left[Nm(\mathbf{\sigma})\right]^2 - Nhm(\mathbf{\sigma}) +\frac{J_0}{2} \\
&\xrightarrow[]{\text{large } N} -\frac{J_0}{2N}\left[Nm(\mathbf{\sigma})\right]^2 - Nhm(\mathbf{\sigma})
\end{align}
$$

where the last term was ignored in the large $N$ limit as it is a constant $\mathcal{O}(N)$.

$$
\therefore H(\mathbf{\sigma}) = -\frac{J_0}{2N}\left[Nm(\mathbf{\sigma})\right]^2 - Nhm(\mathbf{\sigma})
$$

### 2.3.2 Hubbard-Stratonovic Saddle Point

**Aim**:  to calculate the Partition function
$$
\begin{align}
Z &=\sum_{\mathbf{\sigma}}\exp\{-\beta H(\mathbf{\sigma})\} \\
&= \sum_{\mathbf{\sigma}}\exp\left\{\frac{\beta J_0}{2N}\left[Nm(\mathbf{\sigma})\right]^2 + N\beta hm(\mathbf{\sigma})\right\}\\
&=\sum_{\mathbf{\sigma}}\exp\left\{\frac{\beta J_0}{2N}\sum_{i, j}\sigma_i\sigma _j + \beta h \sum_i\sigma_i\right\}
\end{align}
$$

This sum would be easier to calculate if it is possible to factorize the interaction term like
$$
\sum_{\mathbf{\sigma}}g_1(\sigma_1)g_2(\sigma_2)\dots g_N(\sigma_N) = \left[\sum_{\sigma_1 = \pm1}g_1(\sigma_1)\right]\left[\sum_{\sigma_2 = \pm1}g_2(\sigma_2)\right]\dots\left[\sum_{\sigma_N = \pm1}g_N(\sigma_N)\right]
$$

To ease the process, we can use the Hubbard-Stratonovic integral or the Gaussian linearization integral:
$$
\int\frac{d z}{\sqrt{2\pi v}}e^{-\frac{z^2}{2v}}e^{xz} = e^{\frac{1}{2}vx^2}
$$

The Gaussian has one non-zero cumulant which is quadratic.

Now if we substitute $x = Nm(\mathbf{\sigma});\; v = \frac{\beta J_0}{N}$, then the intersection term in $e^{-\beta H(\mathbf{\sigma})}$ can be represented as $e^{\frac{1}{2}vx^2}=e^{\frac{N\beta J_0}{2}m^2(\mathbf{\sigma})}$

$$\begin{align}
\therefore \sum_{\mathbf{\sigma}}e^{-\beta H(\mathbf{\sigma})} &= \sum_{\mathbf{\sigma}}\exp\left\{\frac{\beta J_0}{2}Nm^2(\mathbf{\sigma})+ N\beta h m (\mathbf{\sigma})\right\}\\ \\
&=\sum_{\mathbf{\sigma}}\int \frac{dz}{\sqrt{2\pi\frac{\beta J_0}{N}}}\exp\left\{-\frac{1}{2}\frac{Nz^2}{\beta J_0} + Nzm(\mathbf{\sigma})+ N\beta hm(\mathbf{\sigma})\right\} \\
&= \sum_{\mathbf{\sigma}}\int \frac{dz}{\sqrt{2\pi\frac{\beta J_0}{N}}}\exp\left\{-\frac{1}{2}\frac{Nz^2}{\beta J_0} + Nm(\mathbf{\sigma})\left[z+ \beta h\right]\right\}\\
&= \sum_{\mathbf{\sigma}}\int \frac{dz}{\sqrt{2\pi\frac{\beta J_0}{N}}}\exp\left\{-\frac{1}{2}\frac{Nz^2}{\beta J_0} + \left[z+ \beta h\right]\sum_{\sigma_i}\sigma_i\right\} \\
&= \int \frac{dz}{\sqrt{2\pi\frac{\beta J_0}{N}}}\exp\left\{-\frac{1}{2}\frac{Nz^2}{\beta J_0}\right\} \cdot\sum_{\mathbf{\sigma}}\exp\left\{\left[z+ \beta h\right]\sum_{\sigma_i}\sigma_i\right\} \\
&= \int \frac{dz}{\sqrt{2\pi\frac{\beta J_0}{N}}}\exp\left\{-\frac{1}{2}\frac{Nz^2}{\beta J_0}\right\} \cdot\sum_{\{\sigma_1, \sigma_2\dots\sigma_N\}}\left(e^{\left[\beta h+z\right]\sigma_1}e^{\left[\beta h+z\right]\sigma_2}\cdots e^{\left[\beta h+z\right]\sigma_N}\right) \\
&= \int \frac{dz}{\sqrt{2\pi\frac{\beta J_0}{N}}}\exp\left\{-\frac{1}{2}\frac{Nz^2}{\beta J_0}\right\} \cdot\left(\sum_{\sigma_1 = \{-1, +1\}}e^{\left[\beta h+z\right]\sigma_1}\right)\cdots\left(\sum_{\sigma_N = \{-1, +1\}}e^{\left[\beta h+z\right]\sigma_N}\right) \\
&= \int \frac{dz}{\sqrt{2\pi\frac{\beta J_0}{N}}}\exp\left\{-\frac{1}{2}\frac{Nz^2}{\beta J_0}\right\} \cdot\left(\sum_{\sigma_1 = \{-1, +1\}}e^{\left[\beta h+z\right]\sigma_1}\right)^N \\
&= \int \frac{dz}{\sqrt{2\pi\frac{\beta J_0}{N}}}\exp\left\{-\frac{1}{2}\frac{Nz^2}{\beta J_0}\right\} \cdot\left(e^{\left[\beta h+z\right]} + e^{-\left[\beta h+z\right]}\right)^N\\
&= \int \frac{dz}{\sqrt{2\pi\frac{\beta J_0}{N}}}\exp\left\{-\frac{1}{2}\frac{Nz^2}{\beta J_0}\right\} \cdot\left[2\cosh \left(\beta h+z\right)\right]^N \\
&=\int\frac{dz}{\sqrt{2\pi\beta \frac{J_)}{N}}}e^{-Nf(z)}
\end{align}
$$
where $f(z) = \frac{z^2}{2\beta J_0}  - \ln\left[2\cosh(z+\beta h)\right]$
![[h-s integral.png]]

(*Laplace method/ saddle-point integration*) If now we expand near $f(z^*) := f(z^*) + \frac{1}{2}(z-z^*)^2f^{\prime\prime}(z^*) + \dots$,
then the integral becomes:
$$
\begin{align}
\int dze^{-Nf(z)} &= \int e^{-Nf(z^*)-\frac{N}{2}(z-z^*)^2f^{\prime\prime}(z^*) + \dots}\\ \\
&= e^{-Nf(z^*)}\sqrt{\frac{2\pi}{Nf^{\prime\prime}(z^*)}}
\end{align}
$$

i.e. larger the $N$ , sharper the peak is. And the partition function is described majorly by the value of the function at the minimum point $z^*$.