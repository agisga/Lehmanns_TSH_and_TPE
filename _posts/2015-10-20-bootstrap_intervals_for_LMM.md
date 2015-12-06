---
layout: post
title: Bootstrap confidence intervals and hypotheses tests for linear mixed models
---

Bootstrap confidence intervals and hypotheses tests are introduced in Chapter 15 of TSH. Bootstrap methods are typically used when the small sample distribution as well as the asymptotic sampling distribution of an estimator is not tractable. Also, even when asymptotic tests or confidence regions (such as Wald's) can be obtained, their performance is often provably inferior to the bootstrap methods.

In the following we will consider parametric bootstrap confidence intervals for the fixed effects coefficients in linear mixed models. But first I would like to review the general concept of bootstrap intervals. Later, bootstrap hypotheses tests (likelihood ratio test and goodness of fit) for linear mixed models are introduced as well.

## Bootstrap confidence intervals

Assume that the data vector $X^n  = (x\subscript{1}, \ldots, x\subscript{n})^T \in \mathbb{R}^n$ contains i.i.d. observations from an unknown distribution $P\in\mathcal{P}$. Denote the parameter of interest $\theta(P)$. In order to construct a confidence interval for $\theta(P)$, we consider a real-valued functional $R\subscript{n}(X^n, \theta(P))$ depending on $X^n$ and $\theta(P)$, the so-called *root*. Denote the distribution of the root by $J\subscript{n}(P)$, that is,

$$J\subscript{n}(x, P) = \mathrm{P}(R\subscript{n}(X^n, \theta(P)) \leq x).$$

If the distribution of $R\subscript{n}(X^n, \theta(P))$ is known, then we can construct a (one-sided) $(1-\alpha)$ confidence interval for $\theta(P)$ as

$$\\{\theta\in\Theta : R\subscript{n}(X^n, \theta(P)) \leq J\subscript{n}^{-1}(1-\alpha, P) \\},$$

and a two-sided interval is obtained analogously.

The bootstrap method replaces the unknown distribution $P$ by the empirical distribution $\hat{P}\subscript{n}$. In nonparametric problems the empirical distribution can be defined as $\hat{P}\subscript{n}(E) = \frac{1}{n} \sum\subscript{i=1}^n I\\{x\subscript{i} \in E\\}$. In parametric problems where $\mathcal{P} = \\{P\subscript{\psi} : \psi\in\Psi\\}$ we can use $\hat{P}\subscript{n} = P\subscript{\hat{\psi}\subscript{n}}$. The sample, which is used to compute the empirical distribution, is usually obtained by either resampling observations from $X^n$ with replacement, or in parametric settings with help of random number generators.

In the case that $J\subscript{n}(x, \hat{P})$ is not continuous and strictly increasing, we define

$$J\subscript{n}^{-1}(1-\alpha, \hat{P}) = \inf\\{x : J\subscript{n}(x, \hat{P}) \geq 1 - \alpha\\}.$$

The (two-sided) $(1-\alpha)$ bootstrap confidence interval is consequently given by

$$
\begin{equation}
\label{interval}
B\subscript{n}(1 - \alpha, X^n) = \\{\theta \in \Theta : J\subscript{n}^{-1}(\alpha/2, \hat{P}\subscript{n}) \leq R\subscript{n}(X^n, \theta) \leq J\subscript{n}^{-1}(1-\alpha/2, \hat{P}\subscript{n})\\}.
\end{equation}
$$

### Consistency of the bootstrap estimator

Chapter 15 in TSH presents some important results about the consistency of the bootstrap estimator $J\subscript{n}(\hat{P}\subscript{n})$.

Under the assumption that a continuous limiting distribution exists (i.e. $J\subscript{n}(P\subscript{n}) \rightarrow J(P)$), Theorem 15.4.1 in TSH proves that

1. $\sup\subscript{x} \left|J\subscript{n}(x, P) - J\subscript{n}(x, \hat{P}\subscript{n}) \right| \rightarrow 0$ with probability 1,
2. $J\subscript{n}^{-1}(1-\alpha, \hat{P}\subscript{n}) \rightarrow J^{-1}(1-\alpha)$ with probability 1,
3. $\mathrm{P}\left( \theta(P) \in B\subscript{n}(1-\alpha, X^n) \right) \rightarrow 1-\alpha$, i.e. $B\subscript{n}$ is pointwise consistent in level.

When 1-2 hold, the bootstrap is said to be *strongly consistent*. When 1-2 only hold in probability, then the bootstrap is said to be *weakly consistent*. Property 3 holds for a weakly consistent estimator as well.

Now, consider the root $R\subscript{n}(X^n, \theta(P)) = \sqrt{n}\left(g(\hat{\theta}\subscript{n}) - g(\theta(P)) \right)$, where $\hat{\theta}\subscript{n}$ is an efficient likelihood estimator, and $P$ is a distribution from a quadratic mean differentiable family. In particular, it follows that $R\subscript{n}(X^n, \theta(P))$ is asymptotically normal with mean zero and asymptotic variance $\frac{\partial g}{\partial \theta} \mathcal{I}^{-1}(\theta) \left(\frac{\partial g}{\partial \theta}\right)^T$. Theorem 15.4.2 in TSH proves that bootstrap confidence intervals based on this root are weakly consistent.

## Application to a linear mixed model

Consider a linear mixed model of the general form

$$
\begin{eqnarray}
\label{LMM}
(Y^n | b = \tilde{b}) &\sim& \mathcal{N}(X\beta + Z\tilde{b}, \sigma^2 W^{-1}), \\\\\\
b &\sim& \mathcal{N}(0, \Sigma\subscript{\theta}), \nonumber
\end{eqnarray}
$$

where $Y^n = (y\subscript{1}, \ldots, y\subscript{n})^T \in\mathbb{R}^n$ and $b\in\mathbb{R}^q$ are random vectors (response and random effects), $\beta\in\mathbb{R}^p$ is the vector of fixed effects, and $W\in\mathbb{R}^{n\times n}$ is a diagonal matrix of known prior weights. The random effects covariance matrix $\Sigma\subscript{\theta}\in\mathbb{R}^{q\times q}$ depends on the variance component parameter vector $\theta\in\mathbb{R}^l$.

Let $\hat{\beta}$, $\hat{\theta}$ and $\hat{\sigma}^2$ denote the maximum likelihood estimators of $\beta$, $\theta$ and $\sigma^2$ respectively. I have discussed their derivation in my writeup on [asymptotic tests for linear mixed models](/Lehmanns_TSH_and_TPE/Wald_and_LRT_in_LMM/).

### Bootstrap sample

Parametric bootstrap samples of the parameter estimates in the linear mixed model ($\ref{LMM}$) are obtained by the following procedure (as outlines in [this paper](http://personal.bgsu.edu/~jshang/AICb_assumption.pdf)).

1. Fit a linear mixed model to obtain the estimated fixed effects $\hat{\beta}$, the estimated random effects covariance matrix $\hat{\Sigma} = \Sigma\subscript{\hat{\theta}}$, and the estimated scaling factor (or residual variance) $\hat{\sigma}^2$. 

2. Generate a bootstrap sample as $Y^{\ast} = X\hat{\beta} + Zb^{\ast} + \varepsilon^{\ast}$, where we randomly sample $b^{\ast} \sim N(0, \hat{\Sigma})$ and $\varepsilon^{\ast} \sim N(0, \hat{\sigma}^2 W^{-1})$.

3. Re-fit the linear mixed model to the bootstrap data to obtain bootstrap parameter estimates.

4. Repeat steps 2-3 $N$ times.

We denote the bootstrap parameter estimates obtained in step 3 by $\beta^\ast$, $\theta^\ast$ and $(\sigma^\ast)^2$. We denote the $\alpha$th quantile of the bootstrap sample for a parameter $\omega:=g(\beta, \theta, \sigma^2)$ by $\omega\subscript{\alpha}^\ast$. That is, if $\hat{P}$ is the empirical distribution of $\omega$ obtained from the bootstrap sample, then it holds that $\omega\subscript{\alpha}^\ast = \inf\left\\{x : \hat{P}(x) \geq \alpha\right\\}$.
 
### Basic bootstrap interval

Let $\gamma = \beta\subscript{i}$ denote the $i$th entry of the vector $\beta$, and let $\hat{\gamma} = \hat{\beta}\subscript{i}$ respectively denote the $i$th entry of $\hat{\beta}$. We aim to construct a $(1-\alpha)$ confidence interval for $\gamma$.

Let the root be $R\subscript{n}(Y^n, \gamma) = \sqrt{n}(\hat{\gamma} - \gamma)$. Then the bootstrap interval ($\ref{interval}$) becomes

$$B\subscript{n}(1-\alpha, Y^n) = \left\\{\gamma\in\mathbb{R}^p : \left[\sqrt{n}(\gamma^{\ast} - \hat{\gamma})\right]\subscript{\alpha/2} \leq \sqrt{n}(\hat{\gamma} - \gamma) \leq \left[\sqrt{n}(\gamma^{\ast} - \hat{\gamma})\right]\subscript{1-\alpha/2} \right\\},$$

which simplifies to

$$B\subscript{n}(1-\alpha, Y^n) = \left\\{\gamma\in\mathbb{R}^p : \gamma\subscript{\alpha/2}^{\ast} - \hat{\gamma} \leq \hat{\gamma} - \gamma \leq \gamma\subscript{1-\alpha/2}^{\ast} - \hat{\gamma} \right\\}.$$

The resulting $(1-\alpha)$ confidence interval for $\gamma$ is given by

$$(2\hat{\gamma} -\gamma\subscript{(1-\alpha/2)}^{\ast}, 2\hat{\gamma} -\gamma\subscript{(\alpha/2)}^{\ast}),$$

and is called *basic bootstrap interval* (see also (5.6) in Chapter 5 of A. C. Davison and D. V. Hinkley, Bootstrap Methods and their Application).

### Studentized bootstrap interval

Again, we denote by $\gamma = \beta\subscript{i}$ the $i$th entry of the vector $\beta$, and by $\hat{\gamma} = \hat{\beta}\subscript{i}$ respectively the $i$th entry of $\hat{\beta}$, and again we aim to construct a $(1-\alpha)$ confidence interval for $\gamma$.

We obtain a parametric *studentized bootstrap interval*, which is also known as *bootstrap-t* interval, by considering the root 

$$R\subscript{n}(Y^n, \gamma) = \sqrt{n}\frac{(\gamma^\ast - \gamma)}{\delta^\ast},$$

where $(\delta^\ast)^2$ is a consistent estimator of the variance of the estimator $\gamma^\ast$ (see also (5.7) in Chapter 5 of Davison & Hinkley, Bootstrap Methods and their Application). If likewise, $\hat{\delta}^2$ is a consistent estimator of the variance of the variance of the MLE $\hat{\gamma}$, then the $(1-\alpha)$ bootstrap interval ($\ref{interval}$) becomes

$$B\subscript{n}(1-\alpha, Y^n) = \left\\{\gamma\in\mathbb{R}^p : \left[\sqrt{n}\frac{\gamma^{\ast} - \hat{\gamma}}{\delta^\ast}\right]\subscript{\alpha/2} \leq \sqrt{n}\frac{\hat{\gamma} - \gamma}{\hat{\delta}} \leq \left[\sqrt{n}\frac{\gamma^{\ast} - \hat{\gamma}}{\delta^\ast}\right]\subscript{1-\alpha/2} \right\\}.$$

Denoting $z^\ast = \frac{(\gamma^\ast - \hat{\gamma})}{\delta^\ast}$ it follows that the bootstrap-t interval is given by

$$\hat{\gamma} - \hat{\delta} z\subscript{(1-\alpha/2)}^{\ast} \leq \gamma \leq \hat{\gamma} - \hat{\delta} z\subscript{(\alpha/2)}^{\ast}.$$

### Comparison of methods

Theoretical results given in Chapter 5 of A. C. Davison and D. V. Hinkley, *Bootstrap Methods and their Application* guarantee that for statistics which are approximately normal, the studentized bootstrap confidence intervals are second order accurate, meaning that a confidence interval with confidence level of $(1-\alpha)$ contains the true value with a probability of $(1-\alpha) + \mathcal{O}(n^{-1})$, where $n$ is the sample size. The basic and percentile bootstrap methods however are only first order accurate in general. That is, the interval coverage is correct only up to an order of $n^{-1/2}$. Nevertheless, for equi-tailed confidence intervals (as are all intervals considered above), the basic and percentile methods are second order accurate as well. The Wald Z confidence intervals, which I have introduced in my write up on [asymptotic tests for LMM](/Lehmanns_TSH_and_TPE/Wald_and_LRT_in_LMM/), are first order even when they are equi-tailed. Also note that all theoretical results here assume that the bootstrap sample is sufficiently large.

In general, it appears that basic and studentized bootstrap intervals are superior in accuracy compared to the Wald intervals under any circumstances. Additionally, the studentized bootstrap method adjusts for nonconstant variance and skewness as well as bias.

Of course, the Wald Z method has the advantage of being computationally efficient and convenient. All bootstrap intervals are computationally very heavy, especially for big data sets. Thus, it is probably best to use the Wald Z intervals in the data exploration phase, and compare different kinds of bootstrap intervals once it is more clear what to look for.

Also, see my [blog post](/bootstap_confidence_intervals/) on the implementation of bootstrap confidence intervals in Ruby for the [`mixed_models`](https://github.com/agisga/mixed_models) software package.

## Bootstrap hypotheses tests

Again, denote by $P\in\mathcal{P}$ the unknown distribution of the data. Using the duality of tests and confidence regions, the confidence intervals introduced above can be used to test hypotheses about a parameter $\theta(P)$. That is, given a consistent in level bootstrap interval, a bootstrap hypothesis test of $\mathrm{H} : \theta(P) = \theta\subscript{0}$, which rejects $\mathrm{H}$ if and only if $\theta\subscript{0}$ is outside the confidence interval, is consistent in level. However, not all testing problems can be reduced to testing parameters. For example, consider the goodness of fit problem, where the objective is to test whether the data come from a parametric subfamily of a non-parametric family of distribution.

In general, if we divide $\mathcal{P}$ into two disjoint subsets $\mathcal{P} = \mathcal{P}\subscript{0} \cup \mathcal{P}\subscript{1}$, then we can test

$$\mathrm{H} : P\in\mathcal{P}\subscript{0} \quad\mathrm{vs.}\quad P\in\mathcal{P}\subscript{1}.$$

Given a test statistic $T\subscript{n}$, we want to find a critical value $c\subscript{n, 1-\alpha}$ such that

$$
\begin{eqnarray}
P(T\subscript{n} > c\subscript{n, 1-\alpha}) &\rightarrow& \alpha, \,\mathrm{if}\, P\in\mathcal{P}\subscript{0}, \nonumber \\\\\\
P(T\subscript{n} > c\subscript{n, 1-\alpha}) &\rightarrow& 1, \,\mathrm{if}\, P\in\mathcal{P}\subscript{1}, \nonumber 
\end{eqnarray}
$$

as $n\rightarrow \infty$.

Denote the distribution of $T\subscript{n}$ under $P$ by $G\subscript{n}(P)$ and its quantiles by $g\subscript{n}(1-\alpha, P)$, that is,

$$
\begin{eqnarray}
G\subscript{n}(t, P) &=& P(T\subscript{n} \leq t), \nonumber \\\\\\
g\subscript{n}(1-\alpha, P) &=& \inf\\{t : G\subscript{n}(t, P) \geq 1-\alpha \\}. \nonumber
\end{eqnarray}
$$

Denote by $\hat{Q}\subscript{n}$ the empirical estimate of $P$ under constraint that $P\in\mathcal{P}\subscript{0}$. The bootstrap approach is to estimate the null sampling distribution of $T\subscript{n}$ by $G\subscript{n}(\hat{Q}\subscript{n})$. Then the critical value of the bootstrap test is $g\subscript{n}(1-\alpha, \hat{Q}\subscript{n})$. Theorem 15.6.1 in TSH establishes that this bootstrap test is asymptotically level $\alpha$.

### Bootstrap likelihood ration test

Rather than using the fact that the twice the logarithm of the likelihood ratio has a chi squared distribution under the null (see [my previous writeup on asymptotic tests for LMM](/Lehmanns_TSH_and_TPE/Wald_and_LRT_in_LMM/)), one can bootstrap $T\subscript{n}$ using the parametric resampling strategy introduced before. Thus, if $\left(\hat{\hat{\beta}}, \hat{\hat{\theta}}, \hat{\hat{\sigma}}^2\right)$ are the parameter estimates of model ($\ref{LMM}$) under the null hypothesis, then the critical value of the test is given by $g\subscript{n}\left(1-\alpha, \left(\hat{\hat{\beta}}, \hat{\hat{\theta}}, \hat{\hat{\sigma}}^2\right)\right)$.

Example 15.6.2 in TSH cites a result that the bootstrap LRT has error $\mathcal{O}(n^{-2})$ in rejection probability, while the usual LRT based on the chi squared distribution has an error of $\mathcal{O}(n^{-1})$.

See Section 4.2.3 in Davison and Hinkley (1997) "Bootstrap Methods and their Application" for more detail on the construction and properties of the test.

### Goodness of fit

One might want to test whether the data really comes from a distribution such as assumed in ($\ref{LMM}$).
To this end, let $\hat{P}\subscript{n}$ denote the empirical distribution of $Y^n$ from ($\ref{LMM}$), let $(\hat{\beta}, \hat{\theta}, \hat{\sigma}^2)\subscript{n}$ be the usual maximum likelihood estimators of the parameters in ($\ref{LMM}$), and define the test statistic to be

$$T\subscript{n} = \sqrt{n} \delta\left(\hat{P}\subscript{n}, P\subscript{(\hat{\beta}, \hat{\theta}, \hat{\sigma}^2)\subscript{n}} \right),$$

where $\delta$ is some metric.

We can bootstrap $T\subscript{n}$ based on $P\subscript{(\hat{\beta}, \hat{\theta}, \hat{\sigma}^2)\subscript{n}}$ (same resampling approach as introduced previously), in order to find its sampling distribution under the hypothesis that model ($\ref{LMM}$) is valid. Then a critical value for the test can be readily obtained.
 
See Example 15.6.5 in TSH for references to papers on this method.
