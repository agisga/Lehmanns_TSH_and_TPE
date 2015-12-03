---
layout: post
title: Asymptotic tests and confidence regions for linear mixed models
---

Section 12.4 in TSH discusses the construction of tests and confidence regions, which are based on the asymptotic distribution of the maximum likelihood estimator (MLE). The following explores the application of these approaches to linear mixed models.

### Linear mixed model

A linear mixed model has the general form

$$
\begin{eqnarray}
\label{LMM}
(y | b = \tilde{b}) &\sim& \mathcal{N}(X\beta + Z\tilde{b}, \sigma^2 W^{-1}), \\\\\\
b &\sim& \mathcal{N}(0, \Sigma\subscript{\theta}), \nonumber
\end{eqnarray}
$$

where $y\in\mathbb{R}^n$ and $b\in\mathbb{R}^q$ are random vectors (response and random effects), $\beta\in\mathbb{R}^p$ is the vector of fixed effects, and $W\in\mathbb{R}^{n\times n}$ is a diagonal matrix of known prior weights. The random effects covariance matrix $\Sigma\subscript{\theta}\in\mathbb{R}^{q\times q}$ depends on the variance component parameter vector $\theta\in\mathbb{R}^l$.

The method of maximum likelihood can be used to find parameter estimates $(\hat{\beta}, \hat{\theta})$ and predictions $\hat{b}$ that fit the observed data best. Then, based on the MLE, one can perform hypotheses tests and construct confidence sets regarding the true model parameters.

In the following, we will mostly restrict our attention to the estimation of the fixed effect coefficients $\beta$, and the corresponding significance tests, as well as confidence regions.

### MLE

Without loss of generality we assume that $W=I$ in model ($\ref{LMM}$).
Define $V := \mathrm{Var}(y)$, and notice that

$$V = \mathrm{Var}(y) = \mathrm{Var}(\mathrm{E}(y | b)) + \mathrm{E}(\mathrm{Var}(y | b)) = Z\Sigma\subscript{\theta}Z^T + \sigma^2 I.$$

Thus, model ($\ref{LMM}$) can be rewritten as

$$y \sim \mathcal{N}(X\beta, V),$$

so that the log-likelihood is

$$l = -\frac{N}{2} \log(2\pi) - \frac{1}{2}\log |V| - \frac{1}{2} (y - X\beta)^T V^{-1} (y - X\beta).$$

Differentiation with respect to $\beta$ yields

$$
\begin{equation}
\label{firstderivative}
\frac{\partial l}{\partial \beta} = X^T V^{-1} y - X^T V^{-1} X \beta,
\end{equation}
$$

and consequently the MLE of $\beta$ is

$$
\begin{equation}
\label{MLE}
\hat{\beta} = \left[X^T \hat{V}^{-1} X\right]^{-1} X^T \hat{V} y,
\end{equation}
$$

with $\hat{V} = Z\Sigma\subscript{\hat{\theta}}Z^T + \hat{\sigma}^2 I$ where $\hat{\theta}$ and $\hat{\sigma}^2$ are the MLE of $\theta$ and $\sigma^2$ respectively.

The estimation of $\theta$ and $\sigma^2$ is more involved. It is shown in [Bates et. al. "Fitting linear mixed-effects models using lme4"](http://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf) that the maximum likelihood estimator of $\sigma^2$ is given by

$$\hat{\sigma}^2 = \frac{\\|y - X\beta - Z\mathrm{E}\subscript{\hat{\theta}}(b|y)\\|\subscript{2}^2 + \\|\mathrm{E}\subscript{\hat{\theta}}(u|y)\\|\subscript{2}^2}{n},$$

with $u \sim \mathcal{N}(0, \sigma^2 I)$ such that $b = \Lambda\subscript{\theta} u$ where $\Sigma\subscript{\theta} = \sigma^2 \Lambda\subscript{\theta} \Lambda\subscript{\theta}^T$ (via the Cholesky decomposition). However, in general there is no closed form expression for the MLE of $\theta$, which therefore has to be estimated numerically (except in special cases, such as balanced ANOVA designs).

### Asymptotic normality of the MLE

Differentiating ($\ref{firstderivative}$) with respect to $\beta$ yields $\frac{\partial^2 l}{\partial \beta \partial \beta} = -X^T V^{-1} X$. It follows that the information matrix with respect to $\beta$ is

$$\mathcal{I} = X^T V^{-1} X.$$

Unfortunately, the standard results on asymptotics of MLE, such as described in Chapter 6 of TPE, cannot be used in this case because the observations $y\subscript{i}$ are neither identically distributed nor independent. Nevertheless, it holds that

$$
\begin{equation}
\label{asymptotics}
\hat{\beta} \overset{\mathcal{L}}{\longrightarrow} \mathcal{N}(\beta, \mathcal{I}^{-1}),
\end{equation}
$$

under normality assumptions on the random effects and the error terms (as given by ($\ref{LMM}$)), as well as further regularity assumptions. A proof can be found in [the PhD thesis of J. C. Pinheiro](http://www.math.ku.dk/~erhansen/web/stat1/pinheiro.pdf).

### Wald tests and confidence regions 

Assume that for a fixed $j \in\\{1,\ldots,p\\}$, we want to test the statistical significance of the fixed effect $\beta\subscript{j}$, that is, the hypothesis $\mathrm{H} : \beta\subscript{j} = 0$ against a two-sided alternative. Clearly, H can be written as $g(\beta) = 0$ for a suitable function $g : \mathbb{R}^p \to \mathbb{R}$. 

Thus, we can restrict our attention to the general class of hypotheses tests of the form

$$
\begin{equation}
\label{hypotheses}
\mathrm{H} : g(\beta) = 0 \quad\mathrm{vs.}\quad \mathrm{K} : g(\beta) > 0,
\end{equation}
$$

where the alternative $\mathrm{K}$ may also be two-sided. In the following, we use Wald's approach, as it is illustrated in Section 12.4.2 in TSH, in order to construct a hypotheses test and confidence regions.

An applicatoin of the Delta method to ($\ref{asymptotics}$) yields that

$$g(\hat{\beta}) \overset{\mathcal{L}}{\longrightarrow} \mathcal{N}(g(\beta), \delta^2),$$

where $\delta^2 = \frac{\partial g}{\partial \beta} \mathcal{I}^{-1} \frac{\partial g}{\partial \beta}^T$.

Thus, we obtain a level $\alpha$ test (in an asymptotic sense) for the testing problem ($\ref{hypotheses}$) by rejecting the null $\mathrm{H}$ in...

* ... a one-sided test if and only if $g\left(\hat{\beta}^{(n)}\right) > \hat{\delta}^{(n)} z\subscript{1-\alpha}$, 
* ... a two-sided test if and only if $\left| g\left(\hat{\beta}^{(n)}\right) \right| > \hat{\delta}^{(n)} z\subscript{1-\alpha / 2}$, 

where $\left(\hat{\delta}^{(n)}\right)^2 = \frac{\partial g}{\partial \beta} \left(\hat{\beta}^{(n)}\right) \left(\hat{\mathcal{I}}^{(n)}\right)^{-1} \left( \frac{\partial g}{\partial \beta}\left(\hat{\beta}^{(n)}\right) \right)^T$ with the superscript $n$ denoting the sample size underlying the computations, and where $\hat{\mathcal{I}}^{(n)} = \left(X^{(n)}\right)^T \left(\hat{V}^{(n)}\right)^{-1} X^{(n)}$ is an estimate of $\mathcal{I}$, and $z\subscript{1-\alpha}$ denotes the $(1-\alpha)$ percentile of the standard normal distribution.

By the equivalence of confidence regions and hypotheses tests (see [another of my writeups about TSH](http://0.0.0.0:4000/Lehmanns_TSH_and_TPE/unbiasedness/)), it follows that $g(\beta)$ is contained in the $(1-\alpha)$ confidence region if and only if $g(\beta) \geq g\left( \hat{\beta}^{(n)} \right) - z\subscript{1-\alpha/2}\delta$ and $g(\beta) \leq g\left( \hat{\beta}^{(n)} \right) + z\subscript{1-\alpha/2}\delta$, which yields a Wald confidence interval.

In order to get a Wald confidence ellipsoid for the vector $\beta$, we observe from ($\ref{asymptotics}$) that

$$\mathcal{I}^{1/2} (\hat{\beta} - \beta) \overset{\mathcal{L}}{\longrightarrow} \mathcal{N}(0, I),$$

and consequently

$$(\hat{\beta} - \beta)^T \mathcal{I} (\hat{\beta} - \beta) \overset{\mathcal{L}}{\longrightarrow} \chi\subscript{p}^2.$$

Hence, the Wald confidence ellipsoid is given by

$$\left\\{\beta : (\hat{\beta} - \beta)^T \hat{\mathcal{I}}^{(n)} (\hat{\beta} - \beta) \leq \chi\subscript{p, 1-\alpha}^2 \right\\}.$$

Analogously, a test for $\mathrm{H} : \beta = \beta^\ast$ for a $\beta^\ast \in \mathbb{R}^p$ can be constructed.

### Likelihood ratio test

Theorem 12.4.2 in TSH (also known as Wilk's Theorem) essentially states that twice the logarithm of the likelihood ratio of two nested models converges in law to a Chi squared distribution with the degrees of freedom being the difference in the dimensions of the two (nested) parameter spaces.

However, this Theorem cannot be applied to the linear mixed model ($\ref{LMM}$) directly, because the observations are not identically and independntly distributed. Nevertheless, when testing a linear hypothesis $\mathrm{H} : K^T \beta = 0$  with $K\in\mathbb{R}^{p\times q}$, $\mathrm{rank}(K) = q$ and $q < p$, it can be proven that under normality assumptions on the random effects and the error terms, under the null hypothesis $\mathrm{H}$ it holds that

$$2 \log R\subscript{n} = 2 \log \frac{L}{L\subscript{0}} \overset{\mathcal{L}}{\longrightarrow} \chi\subscript{q}^2,$$

where $L$ denotes the likelihood of the full model ($\ref{LMM}$) and $L\subscript{0}$ the likelihood of the model restricted subject to $K^T \beta = 0$ (it looks like this can be proven by following the approach taking in the proof of Theorem 12.4.2; alternatively it should be possible to derive this result directly from the asymptotic distribution of $\hat{\beta}$ given above; though I haven't carried out either approach completely; however, I have seen this result stated without proof in multiple linear mixed models books, so I trust it (for now) :smile:).

This leads to a simple level $\alpha$ test that rejects $\mathrm{H}$ if and only if

$$2 \log R\subscript{n} > \chi\subscript{q, 1-\alpha}^2,$$

where $\chi\subscript{q, 1-\alpha}^2$ is the $(1-\alpha)$ percentile of $\chi\subscript{q}^2$.
