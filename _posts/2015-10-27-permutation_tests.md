---
layout: post
title: Permutation tests
---

When a parametric probabilistic model cannot be assumed, one can still construct exact level-$\alpha$ hypotheses tests as permutation tests. Here, based on sections 5.8 and 5.9 of TSH, I discuss the concept by considering as an example a permutation test for the difference of two means. 

Assume that each of the random variables $X\subscript{1}, \dots, X\subscript{m}$ has mean $\eta$ and that each of $Y\subscript{1}, ..., Y\subscript{n}$ has mean $\xi$. Additionally assume that the distributions of all those variables differ only with respect to the mean, for example, $X\subscript{i} \sim \mathrm{i.i.d.}\, f(x\subscript{i})$ and $Y\subscript{i} \sim \mathrm{i.i.d.}\, f(y\subscript{i} - \Delta)$ with $\Delta = \eta - \xi$. The density function $f$ is not known apart from the fact that it is continuous a.e. We want to test the hypothesis $H : \Delta = 0$.

Let $N:=n+m$, denote the random vector containing all $X$s and $Y$s as $Z := (X^T, Y^T)^T$, and let $S(z)$ be the set of all permutations of the entries of a realization $z$ of the random vector $Z$.
Then a level-$\alpha$ test $\phi$ has to satisfy

$$\int \phi(z) \prod\subscript{i=1}^N f(z\subscript{i}) dz = \alpha.$$

Interestingly, it turns out that this equality holds if and only if

$$\frac{1}{N!} \sum\subscript{w\in S(z)} \phi(w) = \alpha.$$

A more general result that accounts for population stratification is given by theorem 5.8.1 in TSH.

The power of $\phi$  against an alternative $h(z)$ is given by

$$\int \phi(z) h(z) dz = \int \mathrm{E}\left(\phi(Z) \middle| T=t\right) dP^T(t).$$

Let $T(Z)$ be the order statistic. It holds that $S(z) = S(T(z)) = S(t)$, and from the expression of the conditional expectation $\mathrm{E}\left(\phi(Z) \middle| T=t\right)$ (see Example 2.4.1 and Problem 2.6), it can be further derived that the most powerful test $\phi$ maximizes

$$\sum\subscript{z\in S(t)} \phi(z) \frac{h(z)}{\sum\subscript{w\in S(z)} h(w)}$$

subject to

$$\frac{1}{N!} \sum\subscript{z\in S(t)} \phi(z) = \alpha.$$

Now, the Neyman-Pearson fundamental lemma implies that the hypothesis should be rejected whenever $\frac{h(z)N!}{\sum\subscript{w\in S(z)} h(w)}$ is too large.
This leads to a most powerful test $\phi$ given by

$$\phi(z) = \begin{cases}
1, \quad\mathrm{if}\, h(z) > C(T(z)), \\\\\\
\gamma, \quad\mathrm{if}\, h(z) = C(T(z)), \\\\\\
0, \quad\mathrm{if}\, h(z) < C(T(z)).
\end{cases}$$

Thus the test is carried out by... 

1. ordering the points in $S(z)$ in a decreasing order according to $h$,
2. rejecting if $h(z)$ is one of the $k$ largest values and rejecting with probability $\gamma$ if $h(z)$ is $(k+1)$st largest, where $k$ and $\gamma$ are determined by

   $$k+\gamma = \alpha \cdot N!$$

More general versions of this approach, which incorporate population stratification and randomization, are given in section 5.8-5.13 in TSH.

The above test is not UMP because it depends on $h$. However, it can be shown that if under the null hypothesis each $Z\subscript{i}$ follows the same normal distribution $\mathcal{N}(\xi, \sigma^2)$, then the derived test is most powerful among all unbiased tests of level $\alpha$ against all normal alternatives under consideration (see Lemma 5.9.1 in TSH for an even more general result).
Such an approach is appropriate when the data is assumed to be approximately normal but the assumption is not considered reliable. The permutation test is maximizing the power against all normal alternatives, while still being unbiased against all other alternatives.

