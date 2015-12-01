---
layout: post
title: UMP invariant tests for linear models 
---

Consider a linear model $y = X\beta + \varepsilon$, where $y\in\mathbb{R}^n$, $X\in\mathbb{R}^{n\times p}$ with $p < n$, $\mathrm{rank}(X) = p$, $\beta\in\mathbb{R}^p$, and $\varepsilon\subscript{i} \sim \mathrm{i.i.d.}\, \mathcal{N}(0,\sigma^2)$ for all $i\in\\{1,\dots,n\\}$.

For convenience, denote $\xi := \mathrm{E}(y) = X\beta$. It holds that $\xi\in\Pi$, where $\Pi$ denotes a $p$-dimensional subspace of $\mathbb{R}^n$ (spanned by the columns of $X$).

Assume that we want to test a hypothesis of the form $\mathrm{H} : C^T \beta = 0$, where $C \in \mathbb{R}^{p\times r}$ with $r < p$ and $\mathrm{rank}(C) = r$. It holds that $C^T = B^T X$ where $B^T = C^T(X^T X)^{-1} X^T \in \mathbb{R}^{r\times n}$. Thus, we can rewrite the hypothesis as $\mathrm{H} : B^T \xi = 0$. That is, under the null hypothesis $\xi$ lies in a $(p-r)$-dimensional subspace of $\Pi$. If we denote this $(p-r)$-dimensional subspace by $\omega$, then the testing problem becomes

$$
\begin{equation}
\label{test}
\mathrm{H} : \xi \in \omega \quad\mathrm{vs.}\quad \mathrm{K} : \xi\in\Pi.
\end{equation}
$$

In the following we will construct the uniformly most powerful test among all invariant tests for testing problems of this general form.

## Orthogonal coordinate transformation

Let $Q\in\mathbb{R}^{n\times n}$ be an orthogonal matrix such that its first $p$ rows $q\subscript{1}, \dots, q\subscript{p}$ span $\Pi$, and such that $\mathrm{span}(q\subscript{r+1}, \dots, q\subscript{p}) = \omega$.

Let $z := Qy$, and denote $\eta := \mathrm{E}(z) = Q\xi$. It follows that 

* $z\subscript{p+1} = \dots = z\subscript{n} = 0$ if and only if $y\in\Pi$,
* $z\subscript{1} = \dots = z\subscript{r} = 0$ and $z\subscript{p+1} = \dots = z\subscript{n} = 0$ if and only if $y\in\omega$.

Then the hypothesis H can be given can be given in a simpler form:

$$
\begin{equation}
\label{orthotest}
\mathrm{H} : \eta\subscript{1} = \dots = \eta\subscript{r} = 0.
\end{equation}
$$

## Reduction of the problem via the principle of invariance

In general, we call a testing problem $\mathrm{H} : \xi \in \Omega\subscript{H}$ vs. $\mathrm{K} : \xi\in\Omega\subscript{K}$ *invariant* under a transformation $g$ of the sample space, if the induced transformation $\bar{g}$ on the parameter space preserves both $\Omega\subscript{H}$ and $\Omega\subscript{K}$, that is, if $\bar{g}\Omega = \Omega$, $\bar{g}\Omega\subscript{H} = \Omega\subscript{H}$ and $\bar{g}\Omega\subscript{K} = \Omega\subscript{K}$. See Section 6.1 in TSH for a more detailed definition. In the following, we denote by $G$ the group of all such transformations.

In the present case, the testing problem ($\ref{orthotest}$) remains invariant under the groups of transformations:

* All transformations of the form $z\subscript{i}^\prime = z\subscript{i} + c\subscript{i}$ with any $c\subscript{i}\in\mathbb{R}^{p-r}$ for $i = \\{r+1, \ldots, p\\}$.
* All orthogonal transformations of $z\subscript{1}, \ldots, z\subscript{r}$.
* All scale changes $z^\prime = cz$ with any $c\in\mathbb{R}$.

A function which is constant on each [orbit](https://en.wikipedia.org/wiki/Group_action#Orbits_and_stabilizers) of $G$, but takes on a different value for each orbit is called *maximal invariant*. See Section 6.2 in TSH for a more detailed definition.

In the present case, it is easy to derive (see Section 7.1 in TSH for a step-by-step derivation) that a maximal invariant under the above transformations is given by

$$
\begin{equation}
\label{orthoteststat}
W = \frac{\sum\subscript{i=1}^r z\subscript{i}^2 / r}{\sum\subscript{i=p+1}^n z\subscript{i}^2 / (n-p)},
\end{equation}
$$ 

and a maximal invariant in the parameter space is given by

$$
\begin{equation}
\label{psi}
\psi^2 = \frac{\sum\subscript{i=1}^r \eta\subscript{i}^2}{\sigma^2}.
\end{equation}
$$ 

Thus, the principle of invariance reduces the problem ($\ref{test}$) to

$$
\begin{equation}
\label{simplifiedtest}
\mathrm{H} : \psi^2 = 0 \quad\mathrm{vs.}\quad \mathrm{K} : \psi^2 > 0.
\end{equation}
$$

## Uniformly most powerful invariant test

Since $z \sim \mathcal{N}(\eta, \sigma^2 I)$, the expression in ($\ref{orthoteststat}$) makes it clear that $W \sim F\subscript{n-p}^r$ under the null hypothesis.

Now, Theorem 6.3.2 in TSH implies that the distribution of $W$ only depends on $\psi^2$.
Additionally, it can be shown (see Problems 7.2 and 7.3 in TSH) that the likelihood ratio $\frac{p\subscript{\psi\subscript{1}}(w)}{p\subscript{0}(w)}$ is increasing in $w$ for any $\psi\subscript{1}$. Therefore, by the variant of the Neyman-Pearson fundamental lemma given in Theorem 3.4.1 in TSH, it follows that the uniformly most powerful invariant test of ($\ref{simplifiedtest}$) rejects H if and only if $W > C$, where $C$ is determined by

$$\int\subscript{C}^\infty F\subscript{n-p}^r(w) \mathrm{d}w = \alpha.$$

Moreover, it is worth pointing out that in the case that $r=1$, the test reduces to a two-sided $t$-test

$$
\begin{equation}
\label{orthottest}
t = \frac{|z\subscript{1}|}{\sqrt{\sum\subscript{i=p+1}^n z\subscript{i}^2 / (n-p)}} > C\subscript{0}.
\end{equation}
$$

This $t$-test is not only UMP invariant but also UMP unbiased (see Problem 5.5 in TSH).

### Test statistic in terms of the original variables

Of course, in practice it is rather inconvenient having to find the orthogonal transform $z = Qy$. Therefore, we re-express the test statistic in terms of the original variables $y$.

Let $P = X(X^T X)^{-1}X^T$ and $\hat{\xi} = Py$. That is, $P$ is the orthogonal projection onto $\Pi$, and $\hat{\xi}$ is the least squares estimator of $\mathrm{E}(y) = \xi$.

Then $y - \hat{\xi}$ is orthogonal to $\Pi$ and $\hat{\xi}$ is orthogonal to $\Pi^c$. Therefore and by the orthogonality of $Q$ we have that

$$\sum\subscript{i=p+1}^n z\subscript{i}^2 = \\|Q(y - \hat{\xi})\\|\subscript{2}^2 = \\|y - \hat{\xi}\\|\subscript{2}^2.$$

Similarly we conclude that

$$\sum\subscript{i=1}^r z\subscript{i}^2 + \sum\subscript{i=p+1}^n z\subscript{i}^2 = \\|y - \hat{\hat{\xi}}\\|\subscript{2}^2,$$

where $\hat{\hat{\xi}}$ is the projection of $y$ onto $\omega$.

Thus, we can write the test statistic ($\ref{orthoteststat}$) as

$$
\begin{equation}
\label{teststat}
W = \frac{\left[\\|y - \hat{\hat{\xi}}\\|\subscript{2}^2 - \\|y - \hat{\xi}\\|\subscript{2}^2 \right] / r}{\\|y - \hat{\xi}\\|\subscript{2}^2 / (n-p)} = \frac{\\|\hat{\xi} - \hat{\hat{\xi}}\\|\subscript{2}^2 / r}{\\|y - \hat{\xi}\\|\subscript{2}^2 / (n-p)}.
\end{equation}
$$ 

### Example: Derivation of the well-known t-test

Assume for some $c\in\mathbb{R}^p$ (i.e. $r=1$) we want to test $\mathrm{H} : c^T \beta = 0$. As discussed in the very beginning this test is equivalent to testing $\mathrm{H} : b^T \xi = 0$ where $b^T = c^T(X^T X)^{-1}X^T$. It follows that 

$$
\begin{eqnarray}
\hat{\xi} &=& X(X^T X)^{-1} X^T y, \nonumber \\\\\\
\hat{\hat{\xi}} &=& \hat{\xi} - \frac{b^T y}{\\|b\\|\subscript{2}^2} b = \hat{\xi} - \frac{c^T \hat{\beta}}{\\|b\\|\subscript{2}^2} b, \nonumber 
\end{eqnarray}
$$

where $\hat{\beta} = (X^T X)^{-1} X^T y$ is the usual least squares solution. Then it holds that

$$
\\|\hat{\xi} - \hat{\hat{\xi}}\\|\subscript{2}^2 = \left\\|\frac{c^T \hat{\beta}}{\\|b\\|\subscript{2}^2} b \right\\|\subscript{2}^2 = \frac{(c^T \hat{\beta})^2}{\\|b\\|\subscript{2}^2},
$$

and consequently ($\ref{teststat}$) becomes

$$W = \frac{(c^T \hat{\beta})^2}{s^2 c^T (X^T X)^{-1} c},$$

where $s^2 = \\|y - \hat{\xi}\\|\subscript{2}^2 / (n-p)$. As shown above, $W$ has the $F$ distribution with 1 and $(n-p)$ degrees of freedom. Thus, the statistic

$$t = \frac{c^T \hat{\beta}}{s \sqrt{c^T (X^T X)^{-1} c}}$$

has the $t$ distribution with $(n-p)$ degrees of freedom.

In particular, the hypothesis $\mathrm{H} : \beta\subscript{i} = 0$ is rejected if $|t| > t\subscript{\alpha/2, n-p}$, where the test statistic is given by 

$$t = \frac{\hat{\beta}\subscript{i}}{s \sqrt{\left[(X^T X)^{-1}\right]\subscript{i,i}}}.$$
