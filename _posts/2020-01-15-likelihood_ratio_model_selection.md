---
layout: post
title:  "Asymptotics of the log likelihood ratio and a Bayesian model selection \"paradox\"."
date:   2020-01-15 02:00:00 -0800
categories: bayes
---

In an early draft of Jonathan Huggins' and Jeff Miller's [BayesBag
paper](https://arxiv.org/abs/1912.07104) I learned of a particular ``paradox''
in Bayesian model selection, in which different models with precisely equal
explanatory power are not asymptotically given equal posterior probability by a
Bayesian posterior.  Rather, the Bayes posterior concentrates entirely on one
model or the other, though on which model it concetrates is determined by a fair
coin flip.  In other words, the symmetry between the two equally good models is
maintained on resampling the data, but the Bayesian posterior uncertainty
doesn't seem to adequately represent this symmetry for any particular dataset.
Apparently the idea goes back to a 1966 paper by Berk (see ibid. for the full
citation), so I've been calling it Berk's paradox, though I'm not sure that's
standard. I found that the following simple decomposition of a log likelihood
ratio helped me to understand Berk's paradox, and some other issues besides.

## Notation

Let's introduce some notation.  Suppose we independent and identically
distributed data, $$y = (y_1, \ldots, y_N)$$, and that $$N$$ is large enough for
standard asymptotic approximations to hold.  Say we have two different
probabilistic models for $$y$$ models, parameterized by $$\theta_1$$ and
$$\theta_2$$.  In a standard abuse of notation, we'll write the probabilty of
$$y_n$$ under these two models as $$p(y_n | \theta_1)$$ and $$p(y_n |
\theta_2)$$, respectively.  Identities that hold for either model will be
written in terms of $$p(y_n | \theta)$$. Denote log probabilities by $$\ell(y_n |
\theta) := \log p(y_n | \theta)$$. Finally, suppose that the true distribution
of the data is given by $$y_n \sim q(y_n)$$ (which is not necessary either of
the models), and denote the KL divergence from $$q$$ to $$p$$ by $$KL(q(y_1) ||
p(y_1  | \theta))$$.

Write the true optimal parameters as

$$
\begin{align*}
\bar \theta :=& \mathrm{argmin}_\theta KL(q(y_1) || p(y_1 | \theta)) \\
\end{align*}
$$

and write the maximum likelihood estimators as

$$
\begin{align*}
\hat\theta :=& \mathrm{argmax}_\theta \sum_{n=1}^N \ell(y_n | \theta).\\
\end{align*}
$$

We will be interested in the posterior probability that $$p(y_n | \theta_1)$$
is the data generating process, which we denote $$p(M_1 | y)$$.  For simplicity,
take the equiprobable prior $$p(M_1) = 0.5$$.  By Bayes' rule,

$$
\begin{align*}
p(M_1 | y) =&
    \frac{\int p(y | \theta_1) p(\theta_1) d\theta_1}
         {\int p(y | \theta_1) p(\theta_1) d\theta_1 +
          \int p(y | \theta_2) p(\theta_2) d\theta_2}.
\end{align*}
$$


## Asymptotic decomposition

Let us now do a simple asymptotic analysis of $$p(M_1 | y)$$.
For reasonable priors $$p(\theta)$$, by the Laplace
approximation as $$N \rightarrow \infty$$,

$$
\begin{align*}
\log \int p(y | \theta) p(\theta) d\theta =& \ell(y | \hat\theta) + o(N).
\end{align*}
$$

Define

$$
\Delta(y | \theta_1, \theta_2) := \ell(y | \theta_1)  - \ell(y | \theta_2),
$$

so that

$$
p(M_1 | y) =
    \frac{\exp(\Delta(y | \hat\theta_1, \hat\theta_2))}
    {\exp(\Delta(y | \hat\theta_1, \hat\theta_2)) + 1} + o(N).
$$

We can decompose the key quantity
$$\Delta(y | \hat\theta_1, \hat\theta_2)$$
into three terms, each of a different asymptotic order:

$$
\begin{align*}
\Delta(y | \hat\theta_1, \hat\theta_2) =
    \sum_{n=1}^N \Delta(y_n | \hat\theta_1, \hat\theta_2)
=&
\sum_{n=1}^N \left(
    \Delta(y_n | \hat\theta_1, \hat\theta_2) -
    \Delta(y_n | \bar\theta_1, \bar\theta_2)
    \right) + \\
& \sum_{n=1}^N \left(
    \Delta(y_n | \bar\theta_1, \bar\theta_2) -
    \mathbb{E}_q  \left[ \Delta(y_1 | \bar\theta_1, \bar\theta_2) \right]
    \right) - \\
& N \mathbb{E}_q  \left[ \Delta(y_1 | \bar\theta_1, \bar\theta_2) \right].
\end{align*}
$$

Let us give these terms their own descriptive names:

$$
\begin{align*}
T_{opt} :=& \sum_{n=1}^N \left(
    \Delta(y_n | \hat\theta_1, \hat\theta_2) -
    \Delta(y_n | \bar\theta_1, \bar\theta_2)
    \right) \\
T_{clt} :=& \sum_{n=1}^N \left(
    \Delta(y_n | \bar\theta_1, \bar\theta_2) -
    \mathbb{E}_q  \left[ \Delta(y_1 | \bar\theta_1, \bar\theta_2) \right]
    \right)\\
T_{kl} :=& N
    \mathbb{E}_q  \left[ \Delta(y_1 | \bar\theta_1, \bar\theta_2) \right],
\end{align*}
$$

so that
$$\Delta(y | \hat\theta_1, \hat\theta_2) = T_{opt} + T_{clt} + T_{kl}$$.  As
I now show, in general we expect these terms to be $$O(1)$$, $$O(\sqrt{N})$$,
and $$O(N)$$, respectively.

The first term, $$T_{opt}$$, is the average difference between the log likelihood
at the empirical optima $$\hat\theta$$ and at  the population optima
$$\bar\theta$$.  By an analysis similar to the proof of Wilke's theorem, one can
show that $$T_{opt} = O(1)$$.

The second term, $$T_{clt}$$, is $$O(\sqrt{N})$$ by a standard central limit
theorem argument unless the variance $$V := \mathrm{Var}(\Delta(y_1 |
\bar\theta_1, \bar\theta_2))$$ is zero.

The final term, $$T_{kl}$$, is $$N$$ times the relative KL divergence of the two
models, evaluated at the "true" optimal parameters, since

$$
\begin{align*}
\mathbb{E}_q  \left[ \Delta(y_1 | \bar\theta_1, \bar\theta_2) \right] =&
\mathbb{E}_q  \left[ \log q(y_1) - \ell(y_1 | \bar\theta_2)) \right] -
\mathbb{E}_q  \left[ \log q(y_1) - \ell(y_1 | \bar\theta_1)) \right] \\
=& KL(q || p(y_1 | \bar\theta_2)) - KL(q || p(y_1 | \bar\theta_1)).
\end{align*}
$$

Consequently, $$T_{kl} = O(N)$$ unless the two models have equal $$KL$$
divergence to the truth.

## Discussion

From this decomposition, a number of interest facts follow immediately.

#### Consistency

A consistency result is clearest.  If the $$\theta_1$$ model is a better fit in
terms of KL divergence, then $$KL(q || p(y_1 | \bar\theta_2)) > KL(q || p(y_1 |
\bar\theta_1))$$, and $$T_{kl} \rightarrow \infty$$ at rate $$N$$, dominating
the other terms. Consequently, $$p(M_1 | y) \rightarrow 1$$.  Conversely, if the
$$\theta_2$$ model is a better fit, then $$T_{kl} \rightarrow - \infty$$ and
$$p(M_1 | y) \rightarrow 0$$.  This consistency result obtains despite the
randomness in the data, which is of smaller order.

#### Berk's paradox

The second consequence is Berk's paradox.  Suppose that the two models are
equally explanatory so that $$KL(q || p(y_1 | \bar\theta_1)) = KL(q || p(y_1 |
\bar\theta_2))$$ and, consequently, $$T_{kl} = 0$$.  In that case, the
$$O(\sqrt{N})$$ term $$T_{clt}$$ dominates.  Since $$\frac{1}{\sqrt{N}} T_{clt}
\rightarrow \mathcal{N}(0, V)$$, with equal prbability $$T_{clt}$$ goes to
$$\infty$$ or $$-\infty$$ according to whether $$\frac{1}{\sqrt{N}} T_{clt}$$ is
greater or less than zero, respectively.  Correspondingly, $$p(M_1 | y)$$ goes
to either $$1$$ or zero respectively with equal probability.

In finite samples, Berk's paradox can occur with different severity depending on
the variance $$V$$.  Specifically, if $$\ell(y_n | \theta_1)$$ and $$\ell(y_n |
\theta_2)$$ are highly correlated, then $$V$$ will be small, and $$p(M_1 | y)$$
will approach its extremes relatively slowly.  This occurs when the two models
are very similar to one another --- imagine, for example, comparing two very
similar specific parameter values of a single parametric model.

Note also that, in finite samples, Berk's paradox can occur even if the two
models are not exactly equivalent if the KL difference $$KL(q || p(y_1 |
\bar\theta_1)) - KL(q || p(y_1 | \bar\theta_2))$$ is of order smaller than
$$1 / \sqrt{N}$$.

#### Frequentist testing

Third, one can see clearly why the frequentsist likelihood ratio test comes with
the restriction that the tested models be nested.  Under the null hypothesis
that the two models are equivalent to one another, $$p(y_1 | \bar\theta_1) =
p(y_1 | \bar\theta_2)$$ almost surely, so $$V = 0$$, and $$T_{clt}$$ is
zero. Of course, the null equivalence of the two models implies that $$T_{kl}$$
is zero as well. When, additionally, the two models are equivalent to
$$q(y_1)$$, then the $$T_{opt}$$ term has a $$\chi^2$$ distribution by classical
results (Wilke's theorem), and it is this distribution that is used to test the
null.

Were one to apply the frequentist likelihood ratio test with equivalent but
non-nested models, then $$V$$ would not be zero, and one would also suffer from
the phenomenon of Berk's paradox.  In particular, one would asymptotically
reject the null exactly half the time for any finite cutoff.


## Practical Consequences of Berk's paradox

Is Berk's paradox simply an ivory tower peculiarity?  I do not think it is. On
the contrary, I believe it illustrates a fundamental _brittleness_ of Bayesian
pairwise model selection, particularly when we need Bayesian thinking
most---when there is real model uncertainty, i.e., when $$p(M_1 | y)$$ is not
near $$0$$ or $$1$$.

Suppose you're a workaday Bayesian, and you want to perform a pairwise model
comparison using Bayes' rule.  As you are instructed by your textbooks, you
calculate $$p(M_1 | y)$$.  Suppose your models are complicated so you have no
idea what the posterior's frequentist properties are.  Recall that you don't
actually observe the terms $$T_{opt}$$, $$T_{clt}$$, and $$T_{kl}$$ because you
don't know $$q$$, $$\bar\theta_1$$, or $$\bar\theta_2$$.

What if $$p(M_1 | y)$$ is close to one?  You don't know whether this is the case
because $$T_{kl}$$ or $$T_{clt}$$ has dominated.  In other words, you don't
really know whether the models are nearly equivalent and the high confidence
posterior is due to chance, or because the posterior is extremely confident due
to chance.  Still, you can take some solace from the consistency result that
even small differences in $$KL$$ divergence dominate asymptotically.

What if $$p(M_1 | y)$$ is at some value that is not close to $$0$$ nor to $$1$$?
Then you really don't know much except that the $$KL$$ divergence is probably
small.  Either $$V$$ is very small, and the models have nearly the same
likelihood (in which case comparing them is uninteresting), or $$V$$ is not
small and you would, with all likelihood, not come to the same conclusion with a
new draw from the same generating distribution.

#### Sparse regression example

Here's a concrete example of something I, at least, would have done wrong
before thinking about Berk's paradox in this way.  Suppose you have a sparse
linear model and can actually calculate or estimate the posterior probability of
different sparsity patterns.  Denote a sparsity pattern $M_s$, for $s$ ranging
over the power set of including or excluding regressors.  Before I was aware
of Berk's paradox, I would have happily taken all patterns $$s$$ such that
$$p(M_s | y)$$ was sufficintly close to $$1$$ or $$0$$ and considered those
regressors to be confidently included or excluded.   Then I would have taken
all the others, for which $$p(M_s | y)$$ was at some intermediate value, and
treated them as "uncertain" for further, more careful analysis.  Unfortunately,
as show above, this procedure will be extremely brittle to resampling the data.
One will not pick up the same pattern of "certain" and "uncertain" regressors
with new datasets drawn from the same generating distribution, even if the
data actually comes from one of the specified models.  (One can verify
with simulations that this phenomenon does in fact occur with sparse
linear regression.)


## Paths forward

As I write this, I'm not sure what the solution is.  A careful reading
of the proof of Theorem 4.1 in the BayesBag paper shows that they are using
the bootstrap to estimate $$V$$ and then rescale the test statistic.  I think
that is a good idea, though I do wonder whether the BayesBag posterior is more
than is needed to effect this solution.

I should emphasize that the above objections only apply to the attempt to choose
one model over another, or over a group of alternative models.  Bayesian model
averaging (BMA) does not really suffer from this problem, in that it does not
matter for the purpose of BMA which of several nearly equivalent models you
choose.

So I currently feel a bit lost.  If you are reading this and have ideas or
thoughts, please do reach out!
