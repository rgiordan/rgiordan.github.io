---
layout: post
title:  "Infant sleep training and model selection."
date:   2019-09-15 02:00:00 -0800
categories: bayes
---

We have a one-year old infant who is going through a sleep regression. As my
wife and I discussed how we might approach sleep training, it occurred to me
that choosing a strategy has a lot in common with model selection: even when all
the options are bad, you still need to choose one of them.

Should we follow the books and just let him cry alone for an hour?  But then
he'll probably cry so hard he pukes.  Should we just not worry and rock him
to sleep?  But then we'll be spending our evenings doing that until he's 16.
Should we go in and pick him up but not rock him?  Setting him down again
seems to be extra hard on him.  Should we not pick him up but lay there
so we at least haven't left the room?  That probably won't satisfy anyone,
and we still have to spend our evenings laying on the floor in the dark
for the forseeable future.

None of these options are good.  The problem is that there is no option to
do nothing, and we have to choose the least bad one.

This is precisely analogous with statistical model selection.  Suppose we have a
set of models, $$m = 1, \ldots, M$$, and we want to know which one generated
some data $$y$$ under a likelihood $$p(y | m)$$.  Given a prior probability
$$p(m)$$ of each model, Bayes' rule then says that the most likely model is the
one with the highest

$$
p(m | y) = \frac{p(y | m) p(m)}{\sum_{m'=1}^M p(m' | y)}.
$$

When all of our models are totally lousy --- when $$p(m | y) p(m)$$ is
infinitesimal for all $$m$$ --- by choosing the largest we will still have high
(no less than $$1 / M$$) posterior probability of being correct, because we
normalize by $$\sum_{m'=1}^M p(m' | y)$$ when calculating $$p(m | y)$$.

When thinking about sleep training (or, heaven forid, far more difficult
decisions with our baby), it can be tempting to think you have to choose a good
option when there are no good options.  This perspective can make you feel
trapped, like there is no way to do right by your baby.  But by remembering the
denominator might help you bear in mind that even the best bad thing is still
good from a decision theory perspective.

This whole scenario has an analogy to a classical Bayesian objection to using
naive frequentist hypothesis tests for model selection. Suppose I asked only:
``Should we let him cry alone for an hour?''  Let $$m = 1$$ correspond to this
sleep training strategy.  We find that the answer to this question is ``then
he'll probably cry so hard he pukes'', in other words, that $$p(y | m=1) p(m)$$
is very small, and reject this strategy.  The problem, of course, is that we
didn't consider any of the alternatives before rejecting.
