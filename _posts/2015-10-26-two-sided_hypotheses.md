---
layout: post
title: UMP tests for two-sided hypotheses
---

The (non-) existence of uniformly most powerful (or UMP) tests for two-sided hypotheses is an interesting phenomenon. 

### Example of existence

First let's look at an example when such a test does exist. This is Problem 3.2 in TSH.

For $i = 1,\dots, n$ let $X\subscript{i}$ be i.i.d. $\mathrm{Uniform}(0,\theta)$ random variables, denote their realizations by lower case $x\subscript{i}$s, and let $X$ denote the vector of the $X\subscript{i}$s. Consider the hypothesis $H : \theta = \theta\subscript{0}$ against the alternative $K : \theta \neq \theta\subscript{0}$.

Denote $x\subscript{(n)} := \max\\{x\subscript{1}, \dots, x\subscript{n}\\}$. Let $\phi$ be a hypothesis test which rejects $H : \theta = \theta\subscript{0}$ in favor of a two-sided alternative, if either $x\subscript{(n)} \geq \theta\subscript{0}$ or $x\subscript{(n)} < \theta\subscript{0} \sqrt[n]{\alpha}$.

#### Proof

Using the fundamental lemma of Neyman and Pearson, it is straightforward to prove that $\phi$ is UMP. Namely, $\phi$ is a UMP test at level $\alpha$ by Neyman-Pearson, if for any fixed $\theta\subscript{1} \neq \theta\subscript{0}$, the test $\phi$ can be written as 

$$\phi(x) = \begin{cases} 
1, \quad &\mathrm{if}\, p\subscript{\theta\subscript{1}}(x) > k p\subscript{\theta\subscript{0}}(x),\\\\\\
0, \quad &\mathrm{if}\, p\subscript{\theta\subscript{1}}(x) < k p\subscript{\theta\subscript{0}}(x), 
\end{cases}$$

with a suitable $k$, and if it satisfies

$$\mathrm{E}\subscript{\theta\subscript{0}} \phi(X) = \alpha.$$

We have that

$$\begin{eqnarray}
\nonumber
\mathrm{E}\subscript{\theta\subscript{0}} \phi(X) &=& P\subscript{\theta\subscript{0}}\left(X\subscript{(n)} > \theta\subscript{0}\right) + P\subscript{\theta\subscript{0}}\left(X\subscript{(n)} < \theta\subscript{0}\sqrt[n]{\alpha}\right)\\\\\\
&=& 0 + \left(\frac{\theta\subscript{0} \sqrt[n]{\alpha}}{\theta\subscript{0}}\right)^n = \alpha.
\nonumber
\end{eqnarray}$$

As for the other Neyman-Pearson condition, we have to consider multiple cases:

* If $\theta\subscript{1} > \theta\subscript{0}$, then $k = \left(\frac{\theta\subscript{0}}{\theta\subscript{1}}\right)^n$ yields the desired result.
* If $\theta\subscript{0}\sqrt[n]{\alpha} < \theta\subscript{1} < \theta\subscript{0}$, then $k = \left(\frac{\theta\subscript{0}}{\theta\subscript{1}}\right)^n$ can be used as well.
* If $\theta\subscript{1} < \theta\subscript{0}\sqrt[n]{\alpha} < \theta\subscript{0}$, then $k = 0$.

<div align="right">
$\blacksquare$
</div>

### Example of non-existence

Thus, we saw an example of a UMP test for a two-sided hypothesis. 

However, when the underlying distribution comes from an exponential family, then a UMP test does not exist for $H : \theta = \theta\subscript{0}$ vs. $K : \theta \neq \theta\subscript{0}$ (Problem 3.54 in TSH). This follows quite easily from the consideration of UMP tests for the one-sided hypotheses $H\subscript{1} : \theta \leq \theta\subscript{0}$ vs. $K\subscript{1} : \theta > \theta\subscript{0}$, and $H\subscript{2} : \theta \geq \theta\subscript{0}$ vs. $K\subscript{2} : \theta < \theta\subscript{0}$.
A detailed proof follows.

#### Proof

According to Theorem 3.4.1 in TSH, a UMP test of $H\subscript{1}$ exists and can be written as

$$\phi\subscript{1}(x) = \begin{cases} 
1, \quad &\mathrm{if}\, T(x) > C\subscript{1},\\\\\\
0, \quad &\mathrm{if}\, T(x) < C\subscript{1}.
\end{cases}$$

Similarly, a UMP test of $H\subscript{2}$ exists and can be written as

$$\phi\subscript{2}(x) = \begin{cases} 
1, \quad &\mathrm{if}\, T(x) < C\subscript{2},\\\\\\
0, \quad &\mathrm{if}\, T(x) > C\subscript{2}.
\end{cases}$$

Clearly, $\phi\subscript{1}$ and $\phi\subscript{2}$ are level-$\alpha$ tests for $H$ vs. $K$ as well.

Let $\phi\subscript{0}$ be a level-$\alpha$ test of $H$ vs. $K$. Fix a $\theta\subscript{1} > \theta\subscript{0}$ and a $\theta\subscript{2} < \theta\subscript{0}$. Assume that 

$$\mathrm{E}\subscript{\theta\subscript{i}} \phi\subscript{0}(X) \geq \mathrm{E}\subscript{\theta\subscript{i}} \phi\subscript{i}(X)$$

for $i = 1,2$. Then $\phi\subscript{0}$ is most powerful for testing $\theta\subscript{0}$ vs. $\theta\subscript{1}$ and for testing $\theta\subscript{0}$ vs. $\theta\subscript{2}$. Thus, by the fundamental lemma of Neyman and Pearson the UMP test can be rewritten as

$$
\begin{equation}
\phi\subscript{0}(x) = \begin{cases} 
1, \quad &\mathrm{if}\, p\subscript{\theta\subscript{1}}(x) > k\subscript{1} p\subscript{\theta\subscript{0}}(x),\\\\\\
0, \quad &\mathrm{if}\, p\subscript{\theta\subscript{1}}(x) < k\subscript{1} p\subscript{\theta\subscript{0}}(x), 
\end{cases}
\label{eq1}
\end{equation}
$$

$$
\begin{equation}
\phi\subscript{0}(x) = \begin{cases} 
1, \quad &\mathrm{if}\, p\subscript{\theta\subscript{2}}(x) > k\subscript{2} p\subscript{\theta\subscript{0}}(x),\\\\\\
0, \quad &\mathrm{if}\, p\subscript{\theta\subscript{2}}(x) < k\subscript{2} p\subscript{\theta\subscript{0}}(x).
\end{cases}
\label{eq2}
\end{equation}
$$

Let $x$ be such that $\phi\subscript{0}(x) = 1$. Now, from the monotonicity of the likelihood ratio, it follows that

* if $T(y) > T(x)$ then $\phi\subscript{0}(y) = 1$ (by equation $\eqref{eq1}$),
* if $T(y) < T(x)$ then $\phi\subscript{0}(y) = 1$ (by equation $\eqref{eq2}$).

That is, either $\phi\subscript{0}(y) = 1$ for all $y$ or $\phi\subscript{0}(x) \neq 1$ for all $x$. A contradiction. It follows that $\phi\subscript{0}$ can not be more powerful than $\phi\subscript{1}$ for testing $\theta\subscript{0}$ vs. $\theta\subscript{1}$ and than $\phi\subscript{2}$ for testing $\theta\subscript{0}$ vs. $\theta\subscript{2}$. Thus, a UMP test for $H$ vs. $K$ does not exist.

<div align="right">
$\blacksquare$
</div>

Even though a UMP test for the two-sided hypothesis considered above does not exist, there exist a UMP unbiased test (i.e. a test that is uniformly most powerful among all unbiased tests). For detail see Section 4.2 in TSH.
