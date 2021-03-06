

--------------
--------------

# Old stuff


Of course, for any particular $$\theta$$, by the law of large numbers (LLN),

$$
\frac{1}{N}\sum_{n=1}^N \ell(x_n | \theta) \rightarrow
    \mathbb{E}_{x_1}[\ell(x_1 | \theta)].
$$

The posterior involves integrals over $$\theta$$, however, so we'll need a
stronger result---we need the above limit, and limits like it, to hold uniformly
over all $$\theta$$.  We'll talk about uniform LLNs later as needed.

The second key quantity is the limiting Fisher information,

$$
\mathcal{I} :=
    -\left.\frac{\partial^2 \mathbb{E}_x[\ell(x_1, \theta)]}
               {\partial\theta \partial \theta} \right|_{\theta_0}.
$$

We'll assume that we can exchange expectation and differentiation, and that a
uniform LLN holds for the sample derivatives.  For convenience let's also define

$$
\sigma := \mathcal{I}^{-1/2}.
$$

Finally, we'll need to center and scale the parameter $$\theta$$, so we define

$$
\begin{align*}
\tau &:= \sqrt{N} \frac{\theta - \hat\theta}{\sigma}
    \Leftrightarrow \\
\theta &= \frac{\sigma}{\sqrt{N}} \tau + \hat\theta
\end{align*}
$$

Note that $$\hat\theta$$, as the MLE is a function _only_ of the data.  So,
despite having a theta in the symbol, $$\hat\theta$$ is not a Bayesian
parameter. When we condition on the data, every funciton of the data is as good
as a constant.  Similarly, $$\mathcal{I}$$, as a limiting quantity, is simply a
number, as is $$\sqrt{N}$$.  So the maps between $$\theta$$ and $$\tau$$
are simply an affine reparameterization for any fixed $$x$$ and $$N$$.


### The Taylor series

Let us now examine the Taylor series expansion that will be the basis
of our approximation.  By the mean value theorem applied to a second-order
expansion around some value $$\tilde\theta$$,

$$
\begin{align*}
\sum_{n=1}^N \ell(x_n | \theta) - \sum_{n=1}^N \ell(x_n | \tilde\theta) =
    &\sum_{n=1}^N \ell_{(1)}(x_n | \tilde\theta) (\theta - \tilde\theta) + \\
    &\frac{1}{2}
        \sum_{n=1}^N \ell_{(2)}(x_n | \tilde\theta)
            (\theta - \tilde\theta)^2 + \\
    &\frac{1}{6}
        \sum_{n=1}^N \ell_{(3)}(x_n | \bar\theta)
            (\theta - \tilde\theta)^3,
\end{align*}
$$

for some $$\bar\theta$$ such that
$$|\bar\theta - \tilde\theta| \le |\bar\theta - \theta|$$.  Each of these
terms is still $$O(N)$$, so no LLN can be applied.  But by substititing
subsitituting $$\tilde\tau = \sqrt{N} (\theta - \tilde\theta)$$
(note that the derivatives do not change --- we are only substituting the
arguments), we get

$$
\begin{align*}
\sum_{n=1}^N \ell(x_n | \theta) - \sum_{n=1}^N \ell(x_n | \tilde\theta) =
    &\frac{1}{\sqrt{N}} \sum_{n=1}^N \ell_{(1)}(x_n | \tilde\theta)
        \tilde\tau + \\
    &\frac{1}{2} \frac{1}{N}
        \sum_{n=1}^N \ell_{(2)}(x_n | \tilde\theta) \tilde\tau^2 + \\
    &\frac{1}{6} \frac{1}{N^{3/2}}
        \sum_{n=1}^N \ell_{(3)}(x_n | \bar\theta)
            \bar\tau^3,
\end{align*}
$$

where $$\bar\tau = \sqrt{N} (\theta - \bar\theta)$$.  After this substitution,
we see that, for a particular $$\theta$$, though both individual terms on the
left hand side are $$O(N)$$, their difference is $$O(1)$$, since each term on
the right hand side converges.

Now, by an ordinary central limit theorem (CLT), we expect the first term,
$$\frac{1}{\sqrt{N}} \sum_{n=1}^N \ell_{(1)}(x_n | \tilde\theta)$$, to converge
to a non-degenerate random variable.  So centering at a generic $$\tilde\theta$$
leaves residual randomness.  However by taking $$\tilde\theta = \hat\theta$$,
this term disappears because, by the first-order condition, $$\hat\theta$$ is
the value that sets the gradient of $$\frac{1}{N} \ell_{(1)}(x_n | \theta)$$ to
zero.  This is the reason we center $$\tau$$ at $$\hat\theta$$ and not some
other value.  Plugging in, and recalling the additional re-scaling of $$\tau$$
by $$\sigma$$ gives

$$
\begin{align*}
\sum_{n=1}^N \ell(x_n | \theta) - \sum_{n=1}^N \ell(x_n | \hat\theta) =
    &\frac{1}{2} \sigma^2  \frac{1}{N}
        \sum_{n=1}^N \ell_{(2)}(x_n | \hat\theta) \tau^2 +
    \frac{1}{N^{3/2}} \sigma^3 \frac{1}{6}
        \sum_{n=1}^N \ell_{(3)}(x_n | \bar\theta)
            \bar\tau^3.
\end{align*}
$$

Finally, note that the coefficient of the $$\tau^2$$ term,
$$\sigma^2  \frac{1}{N} \sum_{n=1}^N \ell_{(2)}(x_n | \hat\theta)$$, is a
function of the data alone and, by definition of $$\sigma$$, converges to
$$-1$$, leaving

$$
\begin{align*}
\sum_{n=1}^N \ell(x_n | \theta) - \sum_{n=1}^N \ell(x_n | \hat\theta) =
    -\frac{1}{2} \tau^2 +
    \frac{1}{N^{3/2}} \sigma^3 \frac{1}{6}
        \sum_{n=1}^N \ell_{(3)}(x_n | \bar\theta)
            \bar\tau^3 + ??? \tau^2
\end{align*}
$$


### Uniform laws of large numbers

Uniform LLNs are a complicated subject in
their own right, and I recommend Chapter 18 of _Asymptotic Statistics_ by van
der Vaardt for an introduction.  The key term is that we require the functions
$$\theta \mapsto \ell(x_1, \theta)$$ (or whatever function for which we want
a uniform LLN) to be in a "Glivenko Cantelli" class of functions, meaning
that the map is not too expressive.



### Old text


Under certain circumstances, as $$N \rightarrow \infty$$, the posterior
approaches a particular normal distribution.  Let's take a minute to think about
why this might be.  The quantity inside the exponent, $$\sum_{n=1}^N \ell(x_n |
\theta)$$, in general may be quite complicated as a function of $$\theta$$.
However, for any fixed $$\theta$$, it grows linearly in $$N$$.  After
exponentiating, the highest parts of the function $$\theta \mapsto \sum_{n=1}^N
\ell(x_n | \theta)$$ will be _much_ higher than any surrounding areas.  After
normalizing by the integral in the denominator of $$p(\theta | x)$$, the
posterior will be negligibly small in any diminishing region except around the
maximum.

Consqeuently, the posterior _concentrates_---all its mass goes to a small
region. Now, if the function $$\theta \mapsto \sum_{n=1}^N \ell(x_n | \theta)$$
is smooth, then its behavior in a small region will be well-approximated by a
quadratic function of $$\theta$$, and log quadratic densities are Normal
densities.

Now, it does not suffice to show that quadaratic approximation is good merely
for any fixed $$\theta$$, since posterior quantities typically depend on the
global behavior of $$p(\theta | x)$$ for multiple $$\theta$$.  For example, a
posterior expectations is an integral over $$p(\theta | x)$$.  So we need to
show that integrating against the posterior is asymptotically equivalent to
integrating against a log quadratic density, after appropriate re-scaling.
This, then, will be the whole game---to show that we can Taylor expand the map
$$\theta \mapsto \sum_{n=1}^N \ell(x_n | \theta)$$ to second order, and show
that integrals against the second order approximation are asymptotically the
same as integrals against the full posterior.
