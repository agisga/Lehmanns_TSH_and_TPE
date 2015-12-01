---
layout: post
title: Robustness of hypotheses tests in linear models
---

Consider a linear model $y = X\beta + \varepsilon$, where $y\in\mathbb{R}^n$, $X\in\mathbb{R}^{n\times p}$ with $p < n$, $\mathrm{rank}(X) = p$, and $\beta\in\mathbb{R}^p$. Also, assume that the $y\subscript{i}$ are independent with $\mathrm{Var}(y\subscript{i}) = \sigma^2$. The well-know least squares estimator of $\beta$ is $\hat{\beta} = (X^T X)^{-1} X^T y$. It is the unbiased linear estimator with the minimum variance, and fulfills some other optimality conditions as well (e.g. see [this](/Lehmanns_TSH_and_TPE/LSE_so_nice/) and [that](/Lehmanns_TSH_and_TPE/LSE_so_nice_2/)).

Consider testing the linear hypothesis $\mathrm{H} : c^T \beta = 0$ against a two-sided alternative for some $c\in\mathbb{R}^p$.
In [a previous write-up](/Lehmanns_TSH_and_TPE/small_sample_tests_for_LM/) I have derived a $t$-test of level $\alpha$, which rejects H if $|t| > t\subscript{\alpha/2, n-p}$, where the test statistic is given by

$$
\begin{equation}
\label{ttest}
t = \frac{c^T \hat{\beta}}{\sqrt{\frac{\\|y - X\hat{\beta}\\|\subscript{2}^2}{n-p} c^T (X^T X)^{-1} c}}.
\end{equation}
$$

If it holds that $\varepsilon\subscript{i} \sim \mathrm{i.i.d.}\, \mathcal{N}(0,\sigma^2)$ for all $i\in\\{1,\dots,n\\}$, then this test is UMP invariant (see [my write-up](/Lehmanns_TSH_and_TPE/small_sample_tests_for_LM/)) and UMP unbiased (see Problem 5.5 in TSH).

A reasonable question to ask is whether this test remains a valid level $\alpha$ test (at least in an asymptotic sense) without assuming normality. This property is referred to as *robustness of validity* in the literature (e.g. see p. 421 in TSH).

Without loss of generality we can assume that $\\|c\\|\subscript{2}^2 = 1$. Also, without loss of generality we can assume that the matrix $X$ is orthogonal. That is because we can decompose $X = QR$, where $Q$ is an $n\times p$ matrix with orthogonal columns and $R$ is an upper triangular $p\times p$ matrix. Then defining $z := R\beta$ the linear model is equivalent to $y = Qz + \varepsilon$ and the hypothesis H is equivalent to $\mathrm{H}^\prime : d^T z = 0$, where $d$ is the unique solution to the triangular linear system $R^T d = c$.

Let $d := Xc$. By the orthogonality assumption the numerator of ($\ref{ttest}$) becomes

$$c^T \hat{\beta} = c^T (X^T X)^{-1} X^T y = (Xc)^T y = d^T y.$$

Additionally, under the null hypothesis $\mathrm{H} : c^T \beta = 0$ we have that

$$d^T y = d^T (y - X\beta).$$

The entries of $(y-X\beta)$ are identically and independently distributed with mean $0$ and variance $\sigma^2$. Thus, by Lemma 11.3.3 in TSH (which is a consequence of the Lindeberg Central Limit theorem) it follows that 

$$d^T y \overset{\mathcal{L}}{\longrightarrow} \mathcal{N}(0, \sigma^2),$$

if $\max d\subscript{i}^2 \to 0$ as $n \to \infty$ (see also [my writeup on asymptotic normality of $\hat{\beta}$](/Lehmanns_TSH_and_TPE/LSE_so_nice_2/)). That is, the numerator of the test statistic ($\ref{ttest}$) is asymptotically normal.

Notice that the condition $\max d\subscript{i}^2 \to 0$ is satisfied if and only if $\max\subscript{i} \\|x\subscript{i}\\|\subscript{2}^2\to 0$, because $\max\subscript{i} d\subscript{i}^2 = \max\subscript{i} (x\subscript{i}^T c)^2$ where $x\subscript{i}$ denotes the $i$th row of $X$.

Now, by orthogonality of the columns of $X$ and since $c$ has length 1, the denominator of the test statistic ($\ref{ttest}$) reduces to $\sqrt{\frac{\\|y - X\hat{\beta}\\|\subscript{2}^2}{n-p}}$. This term tends in probability to $\sigma$, which is essentially a consequence of the Weak Law of Large Numbers (see p. 454 in TSH for a step-by-step derivation).

Thus, by Slutsky's theorem it follows that the test statistic $t$ of ($\ref{ttest}$) converges in distribution to $\mathcal{N}(0,1)$. Because the critical value $t\subscript{\alpha/2, n-p}$ converges to the $(1-\alpha/2)$ quantile of the standard normal distribution as well, we conclude that the $t$-test in question is asymptotically robust against non-normality (provided $\max\subscript{i} \\|x\subscript{i}\\|\subscript{2}^2\to 0$).

#### Remark 

Theorem 11.3.1 in TSH establishes the robustness of validity for the more general testing problem $\mathrm{H} : C\beta = 0$ where $C$ is a $q\times p$ matrix (I have derived such a test in [a different writeup](/Lehmanns_TSH_and_TPE/small_sample_tests_for_LM/)).
