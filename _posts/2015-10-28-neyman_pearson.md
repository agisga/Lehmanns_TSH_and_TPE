---
layout: post
title: An informal summary of Neyman-Pearson and generalizations
---

The following offers an informal view on the fundamental lemma of Neyman and Pearson and generalizations thereof.
For a mathematically rigorous presentation see the corresponding results in TSH, which are cited in this article.

*Notation*: MP = "most powerful", UMP = "uniformly most powerful", $H$ denotes the null hypothesis, $K$ denotes the alternative hypothesis, $\alpha$ denotes the level of the hypothesis test, lower case Roman letters denote realizations of random variables (upper case).

1. *[Simple hypotheses]* What is actually called fundamental lemma of Neyman and Pearson in TSH (Theorem 3.2.1) is concerned with a test of two simple hypotheses. Under the null hypothesis the random variable $X$ is assumed to follow a probability distribution with density $p\subscript{0}$, while under the alternative hypothesis the density is $p\subscript{1}$. 

     Consider $H : p\subscript{0}$ vs. $K : p\subscript{1}$. MP test $\phi$ exists. It rejects the null if $\frac{p\subscript{1}}{p\subscript{0}} > k$, accepts the null if $\frac{p\subscript{1}}{p\subscript{0}} < k$, and rejects with probability $\gamma$ if $\frac{p\subscript{1}}{p\subscript{0}} = k$, where $\gamma$ and $k$ are chosen to satisfy $\mathrm{E}\subscript{p\subscript{0}} \phi(X) = \alpha$.

2. *[Monotone likelihood ratio, one-sided, one-param. exp. fam.]* For one-parameter families of distributions and one-sided hypotheses, the Neyman-Pearson lemma can be generalized to construct a UMP test if the distributions in question have monotone likelihood ratios. This is Theorem 3.4.1 in TSH.

     Consider $H : \theta \leq \theta\subscript{0}$ vs. $K : \theta > \theta\subscript{0}$ ($\theta \in \mathbb{R}$). If $\frac{p\subscript{\theta^\prime}(x)}{p\subscript{\theta}(x)}$ is nondecreasing in $T(x)$ for any $\theta < \theta^\prime$, then a UMP test $\phi$ exists. It rejects if $T(x) > C$, accepts if $T(x) < C$, and rejects with probability $\gamma$ if $T(x) = C$, where $C$ and $\gamma$ are determined by $\mathrm{E}\subscript{\theta\subscript{0}} \phi(X) = \alpha$.

     By interchanging the inequalities one obtains a UMP test for the dual problem $H : \theta \geq \theta\subscript{0}$ vs. $K : \theta < \theta\subscript{0}$.

     Additionally, this test minimizes the Type I error subject to $\mathrm{E}\subscript{\theta\subscript{0}} \phi(X) = \alpha$.

3. *[Two-sided null in one-param. exp. fam.]* An analogous UMP test exists for a two-sided null hypothesis $H : \theta \leq \theta\subscript{1} \,\mathrm{or}\, \theta \geq \theta\subscript{2}$ in one-parameter exponential families. It rejects if $C\subscript{1} < T(x) < C\subscript{2}$, accepts if $T(x) < C\subscript{1}$ or $T(x) > C\subscript{2}$, rejects with probability $\gamma$ if $T(x) = C\subscript{i}$ (for $i=1$ or $i=2$), and satisfies $\mathrm{E}\subscript{\theta\subscript{1}} \phi(X) = \alpha = \mathrm{E}\subscript{\theta\subscript{2}} \phi(X)$. Subject to the last condition, this test minimizes the Type I error. See Theorem 3.7.1 in TSH.

     A UMP test for a two-sided alternative hypothesis $K : \theta \leq \theta\subscript{1} \,\mathrm{or}\, \theta \geq \theta\subscript{2}$ does not exist (<a href="{{ site.baseurl }}/two-sided_hypotheses/">e.g. see my corresponding write-up</a>). However, a UMP *unbiased* test analogous to the above exists (see Section 4.2 in TSH).

4. *[UMP unbiased tests, multi-param. exp. fam.]* For multi-parameter exponential families the existence of a UMP test typically cannot be established. However, UMP *unbiased* tests can be constructed without great difficulties. Assume that $\theta\in\mathbb{R}$ is the parameter to be tested, and that $(U, T)$ is a sufficient statistic, where $U$ corresponds to $\theta$ and $T$ corresponds to all other parameters. Then UMP unbiased tests exist for most of the usual hypotheses, and can be written in the same way as in the one-parameter case, except that now all constants specifying the rejection region depend on $T$ (e.g. the rejection rule has the form $u > C(t)$, etc.). Also, the size of the test is measured conditional on $T$.

    See Theorem 4.4.1 in TSH.

5. *[UMP unbiased and independent of sufficient statistic]* UMP unbiased tests for multi-parameter exponential families, as discussed in the last point, are independent of $T$ if a number of additional conditions are satisfied. For example, assume that $V = h(U, T)$ is independent of $T$ (with $\theta = \theta\subscript{1}$ and $\theta = \theta\subscript{2}$) and that $h$ is increasing in $u$. Then a UMP unbiased test for a two-sided null hypothesis rejects if $C\subscript{1} < v < C\subscript{2}$, accepts if $v < C\subscript{1}$ or $v > C\subscript{2}$, etc.

    See Theorem 5.1.1 in TSH for more.

6. *[UMP invariant tests]* If the problem of testing $H : \Omega\subscript{0}$ vs. $K : \Omega\subscript{1}$ remains invariant under a finite group $G = \\{g\subscript{1}, g\subscript{2}, \dots, g\subscript{N} \\}$, then there exists a UMP invariant test that rejects when $\frac{\sum p\subscript{\overline{g}\subscript{i} \theta\subscript{1}} (x)}{\sum p\subscript{\overline{g}\subscript{i} \theta\subscript{0}} (x)} > C$ (for any $\theta\subscript{0} \in \Omega\subscript{0}$ and any $\theta\subscript{1}$ in $\Omega\subscript{1}$). See Theorem 6.3.1 in TSH.

7. *[Finite composite null]* When the null hypothesis specifies that $X$ is distributed according to one of finitely many densities $p\subscript{1}, p\subscript{2}, \dots, p\subscript{m}$, and the alternative hypothesis is $p\subscript{m+1}$, then there exists a test $\phi$ that maximizes $\int \phi p\subscript{m+1} d\mu$. For suitable constants $k\subscript{1}, k\subscript{2}, \dots, k\subscript{m}$, this test rejects the null if $p\subscript{m+1}(x) > \sum\subscript{i=1}^m k\subscript{i} p\subscript{i}(x)$, it accepts the null if $p\subscript{m+1}(x) < \sum\subscript{i=1}^m k\subscript{i} p\subscript{i}(x)$, and it satisfies $\int \phi p\subscript{i} d\mu \leq \alpha$ for $i = 1,2,\dots,m$.

     See Theorem 3.6.1 and Corollary 3.6.1 in TSH for more detail.

8. *[Least favorable distributions]*  Assume a setting similar to the one in the last point, except that the number of distributions under the null hypothesis does not need to be finite. That is, $H : f\subscript{\theta}, \theta \in \omega$ vs. $K : g$.
One can define a *least favorable* distribution $\Lambda$ over $\omega$ and assume that $\theta \sim \Lambda$. As $\Lambda$ is least favorable, one can expect that it leads to a hypothesis test that works best in the worst case (i.e. at values $\theta$ closest to $K$). Thus, $\Lambda$ will typically be a distribution of $H$ that is closest to $K$. In particular, one would have $\Lambda(\omega^\prime) = 1$ for some "boundary region" $\omega^\prime$ of $\omega$. Then a MP test $\phi$ exists. It rejects if $g(x) > k \int f\subscript{\theta}(x) d\Lambda(\theta)$, accepts if $g(x) < k \int f\subscript{\theta}(x) d\Lambda(\theta)$, and satisfies $\sup\subscript{\theta\in\omega} \mathrm{E}\subscript{\theta} \phi(X) = \alpha$.

    See Theorem 3.8.1 and Corollary 3.8.1 in TSH for rigour and detail.

9. *[Maximin tests]* The same approach can also be generalized to test $H : f\subscript{\theta}, \theta \in \omega$ vs. $K : f\subscript{\theta}, \theta \in \omega^\prime$. However, one has to ditch the UMP condition in favor of the condition that $\inf\subscript{\theta\in\omega^\prime} \mathrm{E}\subscript{\theta} \phi(X)$ is maximized under the constraint $\sup\subscript{\theta\in\omega} \mathrm{E}\subscript{\theta} \phi(X) \leq \alpha$. Such a test is called a *maximin* test, and is established in Theorem 8.1.1 and Corollary 8.1.1 in TSH. It rejects if $\int\subscript{\omega^\prime} f\subscript{\theta}(x) d\Lambda^\prime(\theta)(x) > C \int\subscript{\omega} f\subscript{\theta}(x) d\Lambda(\theta)$, for suitably chosen distributions $\Lambda^\prime$ and $\Lambda$.

    Further, Theorem 8.5.1 (Hunt-Stein) and Lemma 8.4.1 in TSH establish the existence of almost invariant tests satisfying the maximin property.
