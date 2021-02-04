---
layout: post
title:  "Introduction: A statistician's guide to functional calculus."
date:   2020-01-27 02:00:00 -0800
categories: functionals
---


With this post I would like to embark on one of the tasks that motivated me to
start this blog: a readable introduction, intended for workaday statisticians,
to differentiable calculus with functionals.  My goal is to write the sort of
thing that would have helped me a couple years ago as I started learning about
this useful and enlightening subject during my PhD. To that end, I hope to
sketch some of the many ways functionals show up in statistics, what new
concepts are needed to extend ordinary Euclidian calculus to spaces of
functions, and provide lots of examples of how infinite-dimensional calculus can
go right and wrong. I may sketch some proofs when they provide intuition, but I
will generally not reproduce rigorous results. I will try, however, to point to
where rigorous results can be found.  My only ambition is to provide context and
motivation for reading the considerable rigorous literature that already exists.

Finally, I should add that my understanding of this topic is quite modest and
far from complete. Writing these posts will be as much to refine my own
understanding as to provide understanding to the reader.  I will surely make
errors, and may start down some dead ends without realizing it.

## Functionals

A "functional" is simply a map from a function to a real number.  It is, in
other words, a function of a function, whose domain is the reals. We will at
times refer to maps from functions to other functions, or, in the extreme, to
functions between generic Banach spaces (which we will define later). Perhaps
confusingly, such exotic maps are still referred to as plain old functions.
As such, in my opinion, the term "functional" seems a bit pretentious and
could perhaps be discarded.

Here, the term "functional calculus" means taking derivatives and forming Taylor
series approximations to functionals.  Of course, such a calculus also requires
thinking carefully about the notions of smoothness and continuity.

For the rest of this post, let's just take a look at some examples of
functionals in statistics without examining their properties just yet.

## The evaluation functional.

The evaluation functional is perhaps the simplest example.  Let $$f: \mathbb{R}
\mapsto \mathbb{R}$$, by which I mean that $$f$$ is a function that takes in a
real and returns a real.  I will also write this as $$f \in \{\mathbb{R} \mapsto
\mathbb{R}\}$$, i.e., $$\{\mathbb{R} \mapsto \mathbb{R}\}$$ means the set of
functions from $$\mathbb{R}$$ to $$\mathbb{R}$$.

Fix some $$x\in \mathbb{R}$$.
The "evaluation functional" $$\phi_x: \{\mathbb{R} \mapsto \mathbb{R}\} \mapsto
\mathbb{R}$$ maps $$f$$ to its value at $$x$$: $$\phi_x(f) = f(x)$$. (If you
think this functional is so simple that nothing much can be said about it,
go read a little bit about reproducing kernel Hilbert spaces.)

Note that we could equally well have written the evaluation functional as
$$\phi: \mathbb{R} \times \{\mathbb{R} \mapsto \mathbb{R}\} \mapsto \mathbb{R}$$
and taken $$\phi(x, f) = f(x)$$.  Which range we choose depends on what we want
to do or emphasize.  Many operations in statistics have functions implicit in
them, and are written with these functions implicity held fixed.  If, instead,
we fix every other aspect of the operation but vary the function, we define a
functional.

## The sample mean.

The sample mean provides quite a rich set of functionals.  Let
us define some independent and identically distributed (IID) data,
$$x = (x_1,\ldots,x_N)$$, and consider the quantity

$$
\frac{1}{N} \sum_{n=1}^N g(x_n).
$$

Let us consider some ways to represent this as a functional by varying
different aspects.

#### The distribution function.

Let us fix $$g$$ and write $$F(\chi) := P(x_n \le \chi)$$
for the cumulative distribution function of $$x_n$$.  Thus, $$F:\mathbb{R}
\mapsto [0,1]$$. Let us write

$$
\phi_c(F) = \int g(x)dF(x).
$$

(I use the subscript $$c$$ for "CDF".) With this functional, $$\phi_c(F) =
\mathbb{E}[g(x)]$$. If we define the empirical distribution to be
$$\hat{F}(\chi) := \frac{1}{N} \sum_{n=1}^N 1(\chi \le x_n)$$, then
$$\phi_c(\hat{F}) = \frac{1}{N} \sum_{n=1}^N g(x_n)$$.  By considering different
arguments to $$\phi_c$$, we can quantify the effect of different distribution
functions on the sample mean.  Similarly, the quantity $$\phi_c(\hat{F}) -
\phi_c(F)$$ measures how close the sample mean is to the population mean.

Functionals that look like $$\phi_c$$ are particularly important in classical
frequentist statistics. One can see how a notion of "continuity" of $$\phi_d$$,
together with a notion of "closeness" of $$\hat{F}$$ and $$F$$, could provide
consistency results.  One can also see how the sensitivity of $$\phi_d(F)$$
to "small" changes of $$F$$ gives a notion of robustness.  Making such ideas
precise is one reason to study functional calculus.

#### The density.

Of course, the distribution function is not the only way to describe a
distribution.  Let $$f:\mathbb{R} \mapsto \mathbb{R}$$ denote the density of
the data $$x$$ with respect to the Lebesgue measure, and write

$$
\phi_d(f) = \int g(x) f(x) dx.
$$

I use the subscript $$d$$ for "density".  The population expectation can be
written as $$\phi_d(f) = \mathbb{E}[g(x)]$$.  Technically the empirical
distribution cannot be represented in this form, but you might imagine fudging
it with the diract delta function $$\delta(\cdot)$$ and $$\hat{f} = \frac{1}{N}
\sum_{n=1}^N \delta(x - x_n)$$.  (Whether such a fudge is legal is definitely
something that we will need to think about carefully.)

Functionals that look like $$\phi_d$$ are particularly important in Bayesian
statistics when considering sensitivity to the prior density, as posterior
expectations can typically be written as ratios of integrals over the prior
density.  Again, some notion of the "smoothness" of a functional like $$\phi_d$$
naturally informs the robustness of a Bayesian posterior, and functional
calculus is one tool to do so.


#### The target function.

If we fix the vector $$x$$, then we can define a functional

$$
\phi(g) = \frac{1}{N} \sum_{n=1}^N g(x_n).
$$

With this definition, the sample mean of $$x$$ is $$\bar{x} = \phi(x \mapsto
x)$$, and the sample variance is given by $$\hat\sigma^2 = \phi(x \mapsto x^2) -
\phi(x \mapsto x)^2$$.  One might imagine choosing some class of functions,
$$\mathcal{G}$$, and finding

$$
g^* := \sup_{g \in \mathcal{G}} \phi(g).
$$

For example, the Rademacher complexity of a function class takes the form of
$$\mathbb{E}\left[\frac{1}{N} \sum_{n=1}^N g^*(x_n)\right]$$.  Though I should
mention that (as of the time of wrigin at least), I don't know much about
applications of functional calculus to such set complexity measures.

## Optimization problems and maximum likelihood estimators.

Suppose we want to solve an optimization problem

$$
\hat\theta = \mathrm{argmin}_\theta \frac{1}{N}\sum_{n=1}^N \ell(\theta, x_n).
$$

As the optimization objective is a sample mean, this can also be written
as a distribution functional:

$$
\phi_d(F) = \mathrm{argmin}_\theta \int \ell(\theta, x) dF(x).
$$

With this definition, $$\phi_d(\hat{F}) = \hat\theta$$.

We might also consider the loss function itself:

$$
\phi_\ell(\ell) = \mathrm{argmin}_\theta \frac{1}{N}\sum_{n=1}^N \ell(\theta, x_n).
$$

We also might consider the effect of a regularizer:

$$
\phi_r(R) = \mathrm{argmin}_\theta \frac{1}{N}\sum_{n=1}^N \ell(\theta, x_n)
+ R(\theta),
$$

with $$\phi_r(\theta \mapsto 0) = \hat\theta$$.

By considering the smoothness of $$\phi_\ell$$ and $$\phi_r$$ we can quantify
the effect on the optimum of changing the optimization problem. For example,
suppose we intend $$\hat\theta$$ to be a maximum likelihood estimator, in which
case $$\ell$$ would be the negative log likelihood.  If we are worried that the
model is misspecified, we might ask how much using a different model class would
have changed our results.

As a special case, suppose we have a generalized linear model with regressors
$$x_n$$, responses $$y_n$$, and an inverse link function $$\gamma$$:

$$
\mathbb{E}[y_n] = \gamma(\theta^T x_n).
$$

A maximum likelihood estimator $$\hat\theta$$ would depend on the functional
form of $$\gamma$$, say through the log likelihood $$\ell_\gamma(\theta, x_n)$$.
We could then write

$$
\phi_\gamma(\gamma) =
\mathrm{argmin}_\theta \frac{1}{N}\sum_{n=1}^N \ell_\gamma(\theta, x_n).
$$

Here, we might imagine restricting $$\gamma$$ to have certain properties, such
as symmetry, unimodality, and so on.

## Generic functionals of the distribution function.

More general functionals of the distribution function $$F$$ provide an important
set of use cases for functional calculus in statistics.  For example, the median
functional can be defined as

$$
\phi_{med}(F) := \inf \left\{x: F(x) \ge \frac{1}{2}\right\}.
$$

One might analyze, for example, $$\sqrt{N}(\phi_{med}(\hat{F}) -
\phi_{med}(F) )$$ and hope that smoothness properties of $$\phi_{med}$$
give rise to a central limit theorem for the median.

## Bayesian statistics

Given a parameter $$\theta$$, data $$x$$, prior $$\pi(\theta)$$ (here with respect
to the Lebesgue measure), and model $$p(x|\theta)$$, the Bayesian posterior
expectation of a quantity $$g(\theta)$$ is

$$
\mathbb{E}[g(\theta) | x]
= \frac{\int g(\theta) p(\theta | x) \pi(\theta) d\theta}
       { \int p(\theta | x) \pi(\theta) d\theta}.
$$

Here, there are loads of functionals to be found.  Many of the
examples we have already considered apply.  For example, we could take the
expectation to be a functional of the model $$p(\cdot | x)$$ or $$p(\theta |
\cdot)$$, of the true distribution, $$F$$, of the data (which would enter if
$$x$$ contains IID data), or of the density of the data, or of the argument
$$g$$.

A special focus in the
Bayesian literature has been given to the map from the prior density to the
posterior expectation:

$$
\phi_{prior}(\pi) =
\frac{\int g(\theta) p(\theta | x) \pi(\theta) d\theta}
       { \int p(\theta | x) \pi(\theta) d\theta}.
$$

We might imagine the prior living within a parametric class or some constrained
non-parametric class.
We could also take other summary statistics of the posterior
which are themselves functionals, such as the median, and write them as
functionals of any of the above inputs, such as the prior:

$$
\phi_{med}(\pi) = \inf \left\{
    \theta:
    \frac{\int_{-\infty}^\theta p(\theta' | x) \pi(\theta') d\theta'}
           { \int p(\theta' | x) \pi(\theta') d\theta'}
            \ge \frac{1}{2}
    \right\}.
$$

## Where next?

The behavior of the functionals encompasses quite a few of the fundamental
concerns of statistics.  Take a functional, $$\phi$$ of the distribution
function, for example.  The asymptotic bias and variance of a statistic $$T =
\phi(\hat{F})$$ is given by the variability of $$\phi(\hat{F})$$ around
$$\phi(F)$$.  Predictive error can be estimated by considering the effect on
$$\hat{F}$$ of creating a held-out set.  Examining robustness to
misspecification of all sorts requires considering varying functional inputs to
statistical procedures.

Often in statistics, functionals can be quite complicated or expensive to
evaluate.  To make matters worse, the set of arguments is infinite dimensional!
Consequently, it can be very useful for both theoretical understanding and for
computation to form local approximations to functionals, as one forms Taylor
series approximations in ordinary calculus.  Since the underlying spaces are
infinite dimensional, it will turn out that we have to be a bit more careful
than in ordinary Euclidian calculus, and some intuition fails to pass over.  But
once one is aware of these issues, one can begin to operate on functionals with
the same mechanical confidence as with ordinary functions.

The first step will be formalzing a generic conception of spaces of functions.
In statistics, this will typically mean extending the domain of our functionals
to more capcious spaces which include some functions which do not represent
proper distributions.  Namely, we will be doing calculus on Banach
spaces---which we will define and describe next time.
