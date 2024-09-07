09.04.2024

**Course**: Statistical Physics of Disordered Systems- From Glassy Energy Landscapes to Machine Learning and Optimization

# 1. Introduction & Motivation:

## 1.1 Optimization problems in Social Networks

**Scenario**: Let there be a party with $N$ friends s.t. not all of the friends go together. Let there exist 2 rooms- one for each group of friends that go well with each other.
**Aim**: To partition people according to their liking in these 2 rooms optimally

*Define*:  Matrix elements:
$$
\begin{align} \\
J_{ij}:
\begin{cases} >0 \;\; \text{if } i, j\; \text{like each other} \\ <0 \;\; \text{if } i, j\; \text{dislike each other} \\=0 \;\; \text{if } i, j\; \text{do NOT know each other OR have no preference}
\end{cases}
\end{align}
$$
Label the rooms  $: \pm 1$

Let assignment of guest $i$ in either room be $\sigma_i : \{+1, -1\}$

Therefore, we have to optimize

$$\text{happiness} = \sum_{i<j}J_{ij}\left(\frac{1+\sigma_i\sigma_j}{2}\right)$$
where $\frac{1+\sigma_i\sigma_j}{2} =\begin{cases}1,\; \; \text{if}\; \sigma_i = \sigma_j\\ 0,\;\; \text{if}\; \sigma_i \neq \sigma_j\end{cases}$ and the sum is over distinct pairs of friends.

i.e. maximizing $\text{happiness}$ will return an optimum vector $\mathbf{\sigma} = (\sigma_1, ..., \sigma_N)$ 

Now, one can try asking some questions:
1. What is the optimal value and how does it depend on the statistics of $[\mathbf{J}]$?
2. How many optimal values exist (if they exist)? Are there many good optima and corresponding $\mathbf{\sigma}$ very different?
3. Is there an efficient way to find this optimal value?
	- Not so easy as there are $2^N$ possible options
4. Is this a hard problem?
	- No, not always! If there are 2 groups s.t. $J_{ij} > 0$ and $J_{ij} = 0$ within one of these groups respectively. Then the solution is trivial.
	- However, for large N ($\approx 10^{23}$) , one can typically get subgroups with all $J_{ij} < 0$.
	- e.g. if $J_{12} < 0, \;\; J_{23} < 0,\;\; J_{13} < 0$, then everyone can NOT be happy simultaneously. This is an example of a **FRUSTRATED** system where one has to find compromises.
## 1.2 Connection to Statistical Physics

The party problem is equivalent to finding the configuration $\mathbf{\sigma}$ with the lowest possible energy. Essentially, it is a minimization problem where one has to minimize:
$$
H(\mathbf{\sigma}) = -\sum_{i<j}J_{ij}\sigma_i\sigma_j$$
where  $J_{ij} = J_{ji}$ and $J_{ii} = 0$

However, if the individual $i$ has a preference $h_i$, then one can modify $H$ as:
$$
\begin{align}
H(\mathbf{\sigma}) &= -\sum_{i<j}J_{ij}\sigma_i\sigma_j -\sum_ih_i\sigma_i$\\
&= -\frac{1}{2}\sum_{i,j}J_{ij}\sigma_i\sigma_j -\sum_ih_i\sigma_i
\end{align}
$$
This problem can be made even simpler by looking at the distribution of the states and NOT at a particular lowest energy state. It's always easy to consider a *Boltzmann distribution* of states
$$P(\mathbf{\sigma})\propto \exp\{-\beta H(\mathbf{\sigma})\}$$

Average energy for $\beta \rightarrow \infty$ will give the lowest energy, but the energy landscape can be very **ROUGH**. essentially, there could be a lot of possible solutions which might not lie close to each other.
e.g.  1. *Ferromagnets*: $J_{ij} > 0 \; \forall \; i,j$ are magnetized in 1 direction.
	2. *Random Field Ising Model*: Preferential attachment, leads to avalanches.
	3. *Spin glass*: Solid state material with randomly placed magner impurities. RKKY interaction.
	4. *Hopfield Model*: Associative memory arising from the structure of $J_{ij}$.
## 1.3 Inference and Inverse problems

**Given Scenario**: No information about $J_{ij}$. The guests decide where they go. *We* observe people during $M$ parties and let the number of observations be $\sigma^m$ where $m\in\{1, \dots, M\}$.
*Define*: Let $P(\mathbf{\sigma}|\mathbf{J})$ be the Boltzmann distribution for a given $\mathbf{J}$ and let $D = \{\mathbf{\sigma^m}\}$ be the data after $m$ parties.

**Aim**: To choose $\mathbf{J}$ that maximizes the likelihood of data $$ P(D|\mathbf{J}) = \Pi_{m=1}^M P(\mathbf{\sigma}^m|\mathbf{J})$$
or minimizes 
$$-\ln P(D|\mathbf{J}) = -\sum_{m=1}^M\ln P(\mathbf{\sigma}^m|\mathbf{J})$$

i.e. we are given $\sigma$ and we have to find $J_{ij}$. The roles have been reversed where the randomness is given by the data $D = \{\mathbf{\sigma}^m\}$ and the optimization is over $\mathbf{J}$.

This exact thing can  be done for proteins by observing tons of amino acids which have preferential attachment to the other amino acids. Then reverse calculate the interaction between such amino acids.

In essence, we can look at the distribution of $\mathbf{J}$ and no the *best* J.

We can also define a **prior** distribution of $\mathbf{J}$ as $P(\mathbf{J})$ and later define the joint probability as:
$$
P(D, \mathbf{J}) = P(D|\mathbf{J})\cdot P{(\mathbf{J})}
$$
where $P(D|\mathbf{J})$ is the **conditional distribution** and $P(D, \mathbf{J})$ is a **generative model** for the ensembles of guests.

Then, $$P(\mathbf{J}|D)P(D) = P(D, \mathbf{J})=   P(D|\mathbf{J})P(\mathbf{J})$$
which gives the **Bayesian estimate/inference** as
$$
\begin{align}
P(\mathbf{J}|D) &\propto P(D|\mathbf{J})P(\mathbf{J})\\\\
\left(\begin{matrix}\text{posterior of }\mathbf{J}\\ \text{after we have seen}\\ \text{the data } D\end{matrix}\right) &=(\text{likelihood})\cdot(\text{prior})
\end{align}
$$

Some questions which arise:
1. How many samples in $D$ do we require to arrive at a good estimate of $\mathbf{J}$?
2. What if the assumption about the distribution $\mathbf{\sigma}^m$ is wrong?

## 1.4 Connection to Machine Learning

So far, outlined cases of **Unsupervised learning**. However, **Supervised learning** infers relation between inputs $\mathbf{x}$ and outputs $\mathbf{y}$ via some *prediction function* $f(\mathbf{x}|\mathbf{y})$.

**Define**: For each example a **LOSS** function $\epsilon(\mathbf{\sigma}, \mathbf{y})$ where $\mathbf{\sigma} = (\mathbf{x}, \mathbf{y})$.
This loss function could be, for example, $\epsilon(\mathbf{\sigma}, \mathbf{J}) = \frac{1}{2}\left[y-f(\mathbf{x},\mathbf{J})\right]^2$

**Learning**/inference can be performed by minimizing the *total error*
$$E(\mathbf{J}|D) = \sum_{m=1}^M\epsilon(\mathbf{\sigma^m}, \mathbf{J})$$

 This can be generalized to a distribution by setting:
 $$
 \begin{align}
 P(\mathbf{J}| D) &\propto \exp\{-\beta E(\mathbf{J}|D)\}P(\mathbf{J})\\ \\
\left(\begin{matrix}\text{Boltzmann-like}\\\text{distribution with disorder}\end{matrix}\right) &= \left(\text{Boltzmann factor}\right)\cdot\left(prior\right)
 \end{align}
 $$
i.e. 2 similar inputs might have completely different outputs.

This can be reduced to the earlier case if 
$$-\beta \epsilon(\mathbf{J}|D) = lnP(D|\mathbf{J})$$
i.e. if the energy function is the *logarithm of **LIKELIHOOD***.

Hence, one can say, all problems, viz. supervised, unsupervised, Bayesian, soft-optimization, have similar structures where the data $D$ usually comes from some distribution.