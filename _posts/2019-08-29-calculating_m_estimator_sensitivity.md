---
layout: post
title:  "A slick (but inefficient) way to calculate M-estimator sensitivity with automatic differentiation."
date:   2019-08-29 20:00:00 -0800
categories: sensitivity
---

I have some recent work [(A Higher-Order Swiss Army Infinitesimal Jackknife)](https://arxiv.org/abs/1907.12116) that is all about
calculating Taylor expansions of optima with respect to hyperparameters
using [automatic differentiation](http://www.jmlr.org/papers/volume18/17-468/17-468.pdf).

In the paper we talk about sensitivity to data weights, but the idea is
much more general.  Suppose you have a parameter which you're optimizing,
$$\theta$$, and some hyperparameter $$\epsilon$$.  For a fixed $$\epsilon$$,
you find $$\hat\theta(\epsilon)$$ to satisfy some vector of first-order
conditions,

$$
G(\hat\theta(\epsilon), \epsilon) = 0.
$$

For example, if you're optimizing $$F(\theta, \epsilon)$$ for a fixed
$$\epsilon$$, $$G$$ would be the vector of partial derivatives of $$F$$ with
respect to $$\theta$$.

Of course, $$\hat\theta(\epsilon)$$ depends on $$\epsilon$$, and you might
want to approximately calculate $$\hat\theta(\epsilon)$$ for different
$$\epsilon$$ without re-solving an optimiztion problem.  Specifically, you
might form a Taylor series expansion around some $$\epsilon_0$$:

$$
\hat\theta(\epsilon) \approx
    \hat\theta(\epsilon_0) +
    \left.\frac{d \hat\theta(\epsilon)}{d\epsilon}\right|_{\epsilon_0}
        (\epsilon - \epsilon_0) +
    \frac{1}{2}
    \left.\frac{d^2 \hat\theta(\epsilon)}{d\epsilon^2}\right|_{\epsilon_0}
        (\epsilon - \epsilon_0)(\epsilon - \epsilon_0) + ...
$$

The difficulty is how to calculate the derivatives, which are defined
implicitly through the solution of $$G(\hat\theta(\epsilon), \epsilon) = 0$$.
The computation section of [our paper](https://arxiv.org/abs/1907.12116)
describes one way to do so recursively, and I've implemented the solution
in the [ParametricSensitivityTaylorExpansion class](https://vittles-python.readthedocs.io/en/latest/api/sensitivity_functions.html#vittles.sensitivity_lib.ParametricSensitivityTaylorExpansion) of my Python package
[vittles](https://github.com/rgiordan/vittles).

My perspective works (and is amenable to theory), but is a bit complicated to
implement.  [Martin Jankowiak](https://www.linkedin.com/in/martin-jankowiak-16789589)
at Uber described to me his idea for an extremely elegant, though unfortunately
inefficient, implementation.  Let me demonstrate his idea in [autograd](https://github.com/HIPS/autograd) and discuss how it is inefficient.

First, we'll need to use the fact that the first derivative is given by

$$
\left.\frac{d \hat\theta(\epsilon)}{d\epsilon^T}\right|_{\epsilon_0} =
    -\left.\frac{\partial G(\theta, \epsilon)}{\partial \theta^T}
        \right|_{\epsilon_0} ^{-1}
    \left.\frac{\partial G(\theta, \epsilon)}{\partial \epsilon^T}
        \right|_{\epsilon_0} ^{-1}.
$$

(Recall that $$G(\theta, \epsilon)$$ is a vector of the same length as
$$\theta$$.)  Now, suppose we have implemented $$G(\theta, \epsilon)$$ in
Python, found a solution `theta0` at `epsilon0`:


```python
def g(theta, epsilon):
    ... your estimating equation here ...

g(theta0, epsilon0) # ...is a vector of zeros.
```

Using this, we can implement the optimal $$\hat\theta(\epsilon)$$ as the
following function, which only evaluates at $$\epsilon_0$$:

```python
def check_epsilon(epsilon):
    assert np.linalg.norm(epsilon - hyperpar0) < 1e-8

@primitive
def get_thetahat(epsilon):
    check_epsilon(epsilon)
    return theta0
```

As-is, this is a useless function.  It only returns what we already know,
which is that `theta0` is the optimum for `epsilon0`, and otherwise throws
an error.  However, we have marked it `@primitive`, which means we can
specify a custom derivative using the formula above.  We will only be
able to evaluate this derivative at `epsilon0`, of course, but that's all
we want.

```python
dg_dtheta = autograd.jacobian(g, argnum=0)
dg_depsilon = autograd.jacobian(g, argnum=1)

# Reverse mode AD is a "vector jacobian product", or "vjp".
def get_thetahat_vjp(ans, epsilon):
    def vjp(g):
        thetahat = get_thetahat(epsilon)
        return -1 * (dg_depsilon(thetahat, epsilon).T @
                     np.linalg.solve(dg_dtheta(thetahat, epsilon), g)).T
    return vjp

# Tell autograd to use get_thetahat_vjp to reverse-mode autodiff get_thetahat.
defvjp(get_thetahat, get_thetahat_vjp)
```

The fucnction `vjp` is simply the reverse-mode implementation of the formula
above for the first derivative.

Now, the magic is that this implemenatation of the derivative of `get_thetahat`
itself is composed of differentiable functions: linear algebra (`solve` and
matrix multiplication), derivatives of `g`, and ... `get_thetahat` itself, whose
derivative we have just defined! Consequently, this single definition suffices
for (reverse mode) automatic differentiation of `get_thetahat` of all orders.
A forward mode implementation is obviously similar.

So how does it work?  It works, but unfortunately, it's quite slow for higher
order derivatives, at least relative to `vittles`.  I believe that our paper
actually makes the reasons for our speedup clear.

* You'll need to repeatedly solve systems involving
$$\left.\frac{\partial G(\theta, \epsilon)}{\partial \theta^T}
\right|_{\epsilon_0} ^{-1}$$
, i.e., calculate
`np.linalg.solve(dg_dtheta(thetahat, epsilon), ...)`.
If you are implementing the derivatives in closed form (as I do in `vittles`)
you can Cholesky factorize this matrix once, but autograd doesn't know that ---
and can't, because it needs to differentate through the matrix evaluation for
this trick to work.
* Lower-order derivatives of $$\hat\theta(\epsilon)$$ appear in the expressions
for higher-order derivatives.  Again, if you are implementing the derivatives
in closed form you can avoid re-caluculating these lower-order derivatives,
but autodiff has no way to be aware of this redundant structure.
* Because every higher order derivative multiplies only terms of the form
$$\epsilon - \epsilon_0$$, there is a lot of redundancy due to symmetry.
Closed-form implementations can recognize this symmetry and avoid redundant
calculations but, again, automatic differentiation is not aware of this
structure.

Although `vittles` overcomes these difficulties, it does so at the cost of
considerable complexity --- as you might expect, since essentially the
benefits come from caching, which is always complicated.

For more details, you can see the notebook below, which is also available
[here](/assets/post_assets/thetahat_ad.ipynb) for download.

{% include_relative thetahat_ad.html %}
