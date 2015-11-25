---
layout: post
title: On the notion of unbiasedness of estimators, hypotheses tests, and confidence intervals
---

The following discusses various well-known definitions of unbiasedness, their generalizations and relationships with each other, as well as some of the underlying intuition (such as the relationship between hypotheses tests and confidence intervals).

## Unbiased estimators

The well-known and widely used definition of an unbiased estimator $\hat{\theta}$ of a parameter $\theta$ is

$$\mathrm{E}\subscript{\theta}(\hat{\theta}) = \theta.$$

However it can be generalized as follows. Assume that there is a loss function $L(\theta, \hat{\theta})$, which only depends on the correct parameter $\theta$ and the estimate $\hat{\theta}$ (i.e. it measures how far off the estimator is from the parameter that it aims to estimate).
Then $\hat{\theta}$ is said to be unbiased for $\theta$ with respect to $L$, if for all $\theta^\prime$ it holds that

$$\mathrm{E}\subscript{\theta}(L(\theta^\prime, \hat{\theta})) \geq \mathrm{E}\subscript{\theta}(L(\theta, \hat{\theta})).$$

That is, if $\hat{\theta}$ is on average closer to the correct parameter $\theta$ than to any wrong parameter $\theta^\prime$ in the parameter space.

When estimating a real valued $\theta$ with the square of the error as loss, the above condition becomes

$$\mathrm{E}\subscript{\theta}\left(\left| \theta^\prime - \hat{\theta} \right|^2\right) \geq \mathrm{E}\subscript{\theta}\left(\left| \theta - \hat{\theta}\right|^2\right).$$

If $\mathrm{E}\subscript{\theta}\hat{\theta}$ is one of the possible values of $\theta$, then by adding and subtracting $\mathrm{E}\subscript{\theta}\hat{\theta}$ inside the parentheses on both sides of the equation it follows that the above unbiasedness condition is satisfied if and only if

$$\mathrm{E}\subscript{\theta}(\hat{\theta}) = \theta.$$

This equivalence also holds under somewhat more general assumptions, see exercise 1.2 in TSH.

## Unbiased tests

Consider a level $\alpha$ test $\phi$ of the hypothesis $H : \theta \in \Omega\subscript{H}$ against an alternative $K : \theta \in \Omega\subscript{K}$.
Denote the power function of $\phi$ by $\beta\subscript{\phi}(\theta) = \mathrm{E}\subscript{\theta} \phi(X)$.
Then it is natural to define unbiasedness of $\phi$ by the criterion

$$
\begin{eqnarray}
\nonumber
\beta\subscript{\phi}(\theta) &\leq& \alpha \quad \mathrm{if}\, H : \theta \in \Omega\subscript{H}, \\\\\\
\beta\subscript{\phi}(\theta) &\geq& \alpha \quad \mathrm{if}\,  K : \theta \in \Omega\subscript{K}. 
\nonumber
\end{eqnarray}
$$

In particular, it follows that $\beta\subscript{\phi}(\theta) = \alpha$ on the common boundary of $\Omega\subscript{H}$ and $\Omega\subscript{K}$. In fact, a test that is the most powerful among all such tests, is UMP unbiased (Lemma 4.1.1 in TSH). 

However, the definition of an unbiased test can be generalized in the same way as that of an unbiased estimator shown above.
Assume that there is a loss function $L(\theta, \phi(x))$, which only depends on the true value of $\theta$ and the decision $\phi(x)$ takes by the test $\phi$. Then the hypothesis test is unbiased with respect to $L$, if for all $\theta^\prime$ it holds that

$$\mathrm{E}\subscript{\theta}(L(\theta^\prime, \phi(X))) \geq \mathrm{E}\subscript{\theta}(L(\theta, \phi(X))).$$

For the test $\phi$ of $H$ vs. $K$ let the loss function be equal to $\alpha$ if a Type II error is committed and equal $(1-\alpha)$ if a Type I error is committed. Then 

$$
\mathrm{E}\subscript{\theta}(L(\theta^\prime, \phi(X))) = 
\begin{cases}
\alpha (1 - \beta\subscript{\phi}(\theta)) \quad &\mathrm{if}&\, \theta^\prime \in \Omega\subscript{K}\\\\\\ 
(1-\alpha) \beta\subscript{\phi}(\theta) \quad &\mathrm{if}&\, \theta^\prime \in \Omega\subscript{H},
\end{cases}
$$

It follows that if $\theta \in \Omega\subscript{H}$ then $\alpha (1 - \beta\subscript{\phi}(\theta)) \geq (1-\alpha) \beta\subscript{\phi}(\theta)$, and consequently

$$\beta\subscript{\phi}(\theta) \leq \alpha.$$

Similarly, by considering $\theta\in\Omega\subscript{K}$, we get $\beta\subscript{\phi}(\theta) \geq \alpha$. Thus the usual definition is a special case of the more general loss-function-based definition.

## Unbiased confidence sets

As is well-known, the defining condition for a confidence interval $\left(\underline{\theta}, \overline{\theta}\right)$ is

$$P\subscript{\theta}\left(\underline{\theta}(X) \leq \theta \leq \overline{\theta}(X)\right) \geq 1-\alpha,$$

for all $\theta$.

### Hypotheses tests vs. confidence intervals

It is well-known that hypotheses tests and confidence intervals generally do exactly the same thing.
However, to describe with mathematical rigour in what sense it is true requires a little thinking.

Consider a level $\alpha$ test of a two-sided hypothesis test $H : \theta = \theta\subscript{0}$ vs. $K : \theta \neq \theta\subscript{0}$, and denote its acceptance region by $A(\theta\subscript{0})$.
Define the inclusion region of the confidence set to be

$$S(x) := \\{ \theta : x\in A(\theta) \\},$$

that is, $\theta \in S(x)$ if and only if $x\in A(\theta)$. Then $S(x)$ defines a $(1-\alpha) \cdot 100\\%$ confidence set, because for all $\theta$ we have

$$P\subscript{\theta}(\theta \in S(x)) = P\subscript{\theta}(x\in A(\theta)) \geq 1 - \alpha.$$

Conversely, if we start out with a family of confidence sets $\\{S(x) : x\in\mathcal{X}\\}$, and define $A(\theta) := \\{x : \theta\in S(x)\\}$, then for any $\theta$ it holds that

$$P\subscript{\theta}(x\in A(\theta)) = P\subscript{\theta}(\theta \in S(x)) \geq 1 - \alpha.$$

It follows that $P\subscript{\theta}(\mathrm{Type\,I\,error}) \leq \alpha$, that is, $A(\theta)$ is the acceptance region of a level $\alpha$ test.

### Unbiased and uniformly most accurate unbiased confidence sets

Now it suggests itself to define an unbiased confidence set as one that stems from an unbiased hypothesis test by the above procedure. 
In the two-sided case discussed above this condition reduces to

$$P\subscript{\theta}\left(\underline{\theta}(X) \leq \theta^\prime \leq \overline{\theta}(X)\right) \leq 1 - \alpha$$

for all $\theta^\prime$ and $\theta$ such that $\theta \neq \theta^\prime$. That is, the inclusion probability of the null hypothesis parameter $\theta^\prime$ in the confidence interval, when the alternative $\theta$ is true, is less than the confidence level. Lemma 5.5.1 in TSH shows that the confidence set derived from an unbiased level $\alpha$ hypothesis test has indeed the form of an interval.

Similarly, uniformly most accurate confidence intervals correspond to uniformly most powerful tests (see section 3.5 in TSH for more detail).
However, UMP tests usually do not exist, which is a reason to concentrate on unbiasedness instead. In particular, UMP unbiased tests correspond to uniformly most accurate unbiased confidence sets, i.e.  $S(x)$ such that for all $\theta^\prime$ and $\theta$ with $\theta\in K(\theta^\prime)$ the probability $P\subscript{\theta}(\theta^\prime\in S(x))$ is minimized.
