# Exact Conditional Poisson Regression: Theory and Implementation

------------------------------------------------------------------------

*This document explains the theory behind every exact-inference
algorithm in this package:*

- *[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md):
  a fast, closed-form method for the special case of one two-level
  effect with strata as the only other structure.*
- *[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md):
  a general method that allows any number of categorical covariates and
  strata.*
- *[`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md):
  a general method that computes the same mathematics again, using a
  different strategy that prunes more aggressively.*

*All three are independently implemented, yet compute exactly the same
conditional distribution; this document derives why that is true, and
works through exactly how each of the three arrives at it.*

*For a worked, runnable example that puts this theory to use, such as
simulating a low-incidence vaccine-efficacy trial, comparing the
package’s exact inference against
[`glm()`](https://rdrr.io/r/stats/glm.html)’s asymptotic inference, and
reusing an already-computed exact distribution to answer a different
question without recomputing it, see the companion vignette,
[`vignette("Example")`](https://gsk-biostatistics.github.io/poissonexact/articles/Example.md).*

------------------------------------------------------------------------

## 1. Motivation: why exact conditional inference

Poisson regression with a log link is the standard way to model counts
and event rates. The usual way of testing hypotheses and building
confidence intervals for it (Wald tests, likelihood-ratio tests,
asymptotic confidence intervals) relies on large-sample approximations.
These approximations describe what happens as the sample size grows very
large, and they can be poor when:

- sample sizes are small,
- the outcome is a rare event,
- the model has many categorical levels or many strata, each with few
  events,
- an estimated rate ratio is close to the edge of what the data can tell
  apart from zero or infinity (quasi-complete separation).

**Exact conditional inference** avoids large-sample approximations
entirely. Instead of relying on how an estimator behaves as the sample
grows, it builds the *exact* distribution of a test statistic for the
sample size actually observed. It does this by eliminating every
nuisance parameter through conditioning on that parameter’s sufficient
statistic. Every p-value, confidence interval, and point estimate this
produces is exact, not an approximation. This is the same idea behind
Fisher’s exact test for 2x2 tables and exact logistic regression (Cox
1970; Hirji, Mehta and Patel 1987; Mehta and Patel 1995).

This document works through the same idea for the Poisson case. Section
[3](#sec3) first states the general principle it rests on. Sections
[5](#sec5) through [7](#sec7) then describe three algorithms that
compute it: one fully general algorithm (Section [5](#sec5)), one
restricted to an important special case where a closed form is available
(Section [6](#sec6)), and one more that computes the same general-case
answer as Section [5](#sec5), but by a different, more aggressively
pruned route (Section [7](#sec7)).

------------------------------------------------------------------------

## 2. The model and notation shared by every algorithm

### 2.1 Indices and strata

Let i \in \\1,\dots,n\\ index the observations. Observations can
optionally be split into Q disjoint strata, I_s \subseteq \\1,\dots,n\\,
\qquad I_s \cap I\_{s'} = \varnothing \\ (s \ne s'), \qquad
\bigcup\_{s=1}^Q I_s = \\1,\dots,n\\. Write s(i) for the stratum that
contains observation i. If no stratification variable is supplied, Q = 1
and every observation belongs to one pseudo-stratum called “ALL.”
Because of this convention, every statement below about “strata” also
covers the unstratified case, and neither the code nor this document
needs a special case for it.

### 2.2 Response, offset, exposure

The observed response is a vector of non-negative integer counts
y^{\mathrm{obs}} = (y_1^{\mathrm{obs}}, \dots, y_n^{\mathrm{obs}}). An
optional offset o_i \in \mathbb{R} (a log-exposure) defines the exposure
E_i := \exp(o_i) \> 0. If no offset is supplied, every function in this
package sets o_i = 0 and E_i = 1 for every observation (one unit of
exposure each), so an unspecified offset never needs special handling
later.

### 2.3 Design: the interest column and the nuisance columns

Every algorithm in this package does exact inference for **one scalar
coefficient at a time**. (The general algorithms can repeat this for
several requested effects, one after another.)

**Notation convention used from here on.** A bold letter denotes a
vector or matrix; the same letter in ordinary (non-bold) type denotes a
single scalar entry of it, normally picked out by a subscript. For
example \mathbf{z} = (z_1,\dots,z_n) is the whole interest column and
z_i is its i-th entry; \boldsymbol{\gamma} = (\gamma_1,\dots,\gamma_p)
is the whole nuisance coefficient vector and \gamma_k is its k-th entry.
The handful of quantities that are conventionally written as plain
(non-bold) letters even though they collect several numbers (\tau =
(\tau_1,\dots,\tau_Q), the vector of stratum intercepts, and
\theta,\beta,\lambda_s, which are always scalars) follow the established
convention of the literature this document cites (Guess and Thomas 1990;
Hirji, Mehta and Patel 1987) rather than this rule, and are called out
individually wherever they first appear.

- \mathbf{z} \in \mathbb{Z}^n is the *interest* design column, with i-th
  entry z_i, typically a 0/1 dummy for one level of a categorical
  effect. Its coefficient \beta is the target of inference.
- X_N \in \mathbb{Z}^{n \times p} is the *nuisance* design matrix (“N”
  for “nuisance,” to distinguish it from the sample size n): every other
  model effect that is not currently being tested, collected as columns.
  This can include other categorical covariates, as well as other dummy
  columns of the same multi-level factor that \mathbf{z} belongs to. Row
  i of X_N, written \mathbf{x}\_{N,i} \in \mathbb{Z}^p, is the
  nuisance-column values for observation i; its k-th entry is
  x\_{N,i,k}. The nuisance coefficient vector paired with X_N is
  \boldsymbol{\gamma} \in \mathbb{R}^p.

**Scope: categorical effects only.** Both \mathbf{z} and every column of
X_N must be small-integer-valued: concretely, values in \\0,1,\dots,K\\
for some small K, almost always K=1 (a 0/1 dummy). In practice this
means the 0/1 dummy columns produced by treatment-coding a categorical
factor (Sections [5.1](#sec5-1) and [6.1](#sec6-1) give the exact
details). The dynamic program, convolution, and network recursion in
Sections [5](#sec5)-[7](#sec7) all work by explicitly listing every
value a column can take, one at a time; their cost grows with the number
of distinct values (roughly linearly per nuisance coordinate, and
multiplicatively across coordinates for the joint table of Section
[3.4](#sec3-4)). A handful of categorical levels (the ordinary case: a
handful of treatment arms or a handful of ordered bins) keeps this
listing cheap. **This implementation does not support continuous
covariates** (for example, a continuous biomarker, or age measured in
years) as either the interest effect or an adjustment covariate. A
continuous variable has to be converted to a categorical one first
(grouped into a small number of ordered levels) before it can be used
with any of these functions. Section [4.7](#sec4-7) notes the caveat
that comes with this: grouping a continuous covariate into bins throws
away information compared to modelling it on its original scale.

### 2.4 The stratified Poisson log-linear model

Y_i \sim \mathrm{Poisson}(\mu_i), \qquad \mu_i = E_i\exp(\eta_i), \qquad
\eta_i = \tau\_{s(i)} + \beta z_i + \boldsymbol{\gamma}^\top
\mathbf{x}\_{N,i}, equivalently \log \mu_i = o_i + \tau\_{s(i)} + \beta
z_i + \boldsymbol{\gamma}^\top \mathbf{x}\_{N,i}. \qquad\text{(2.1)}
\tau\_{s(i)} is not a separate intercept for each observation. It means
“observation i uses the intercept that belongs to its own stratum s(i).”
The full intercept vector \tau = (\tau_1,\dots,\tau_Q) has one entry per
stratum. Every stratum has its **own free parameter**: there is no
shared reference level, and no assumption that intercepts are related
across strata in any way. This point matters a great deal in Section
[4](#sec4). Every nuisance quantity here (the Q stratum intercepts and
the p nuisance coefficients) gets eliminated by conditioning, as
explained next.

**One modelling choice is easy to miss, because it is about what is
absent from (2.1) rather than what is present: \beta itself carries no
stratum index.** It is one shared number, used identically for every
observation no matter which stratum it falls in; it is not
\beta\_{s(i)}, the way \tau\_{s(i)} is. This is the general version of
an assumption that Guess and Thomas (1990) state explicitly for the
two-group case handled by this package’s closed-form algorithm (Section
[6.1](#sec6-1)): “the relative risk is assumed to be constant across
strata.” Every algorithm in this package makes this same assumption for
whichever effect is being tested, because every algorithm starts from
the same model (2.1). It is not a simplification specific to the
closed-form special case, and using the fully general multivariate-shift
or network algorithms instead (Sections [5](#sec5) and [7](#sec7)) does
not relax it.

It is worth being precise about why giving every stratum its own free
intercept \tau_s does *not*, by itself, relax this assumption. The two
are separate, independent pieces of (2.1). \tau_s sets each stratum’s
own baseline event rate at the reference level of every tested or
adjusted effect, and is free to differ arbitrarily from one stratum to
the next. \beta sets how much that rate changes when z_i changes, and
(2.1) says that this *change* is the same number in every stratum, no
matter how different their baselines are. Making the intercepts more
flexible (more strata, a finer partition, even one intercept per
observation) never touches \beta’s own single coefficient. It only
changes how much nuisance structure gets conditioned away (Section
[3](#sec3)). A model can have arbitrarily many, arbitrarily different
stratum baselines and still assume one constant effect on top of them:
that is exactly what (2.1) specifies. (This is the same “common effect”
assumption used in Mantel-Haenszel-type common-odds-ratio estimation and
in proportional-hazards models with stratified baselines. It is well
established elsewhere and is not a limitation specific to this package.)

If a genuinely stratum-varying effect *is* wanted (for example, the rate
ratio for \mathbf{z} differing by region), it has to be modelled
explicitly, as an interaction between \mathbf{z} and the variable it is
allowed to vary by. (`exact_vars` can name such an interaction directly;
see Section [5](#sec5)’s discussion of interaction-term naming, and
Section [4.3](#sec4-3)’s composite-variable construction, for how a term
like this turns into model-matrix columns.) Doing so replaces the single
\beta with one coefficient per level of the interacting variable. Each
of those coefficients is then just another effect being tested under
this package’s usual assumption: constant across whatever grouping
remains once the interaction’s own levels are accounted for. There is no
way to ask this package for a coefficient that is simultaneously “the
effect of \mathbf{z}” and “different in every stratum” without naming
that difference as its own term. Simply adding more strata never does
it.

### 2.5 Notation at a glance, and one deliberate disambiguation

Several symbols below carry more than one index. It is worth fixing,
once and for all, which quantities are single numbers per stratum and
which are numbers per stratum *and* per nuisance coordinate, because two
of them are easy to mix up at a glance:

| Symbol | Meaning | Indexed by |
|----|----|----|
| M_s^{\mathrm{obs}} | The **observed total count of stratum s**, \sum\_{i\in I_s} y_i^{\mathrm{obs}}. This is the one number each stratum’s intercept gets conditioned on. | stratum s only, one number per stratum. |
| \mathbf{t}\_{N,\mathrm{obs}} | The observed nuisance statistic vector, X_N^\top y^{\mathrm{obs}} (Section [2.3](#sec2-3)) | nuisance coordinate k=1,\dots,p only, shared across every stratum, not stratum-specific |
| S_i | The per-observation cap on y_i (Section [5.2](#sec5-2)) | observation i only |
| U\_{s,k} | The **most** stratum s could possibly contribute to nuisance coordinate k, given only its own observations and its own total M_s^{\mathrm{obs}} (defined in Section [5.5](#sec5-5)) | **both** stratum s and coordinate k |
| L\_{s,k} | The **least** stratum s could possibly be forced to contribute to nuisance coordinate k (defined in Section [5.5](#sec5-5); used only by the network algorithm, Section [7](#sec7)) | **both** stratum s and coordinate k |

Here is why this table is worth spelling out. U\_{s,k} and L\_{s,k} are
each defined as a `min`/`max` expression with M_s^{\mathrm{obs}} as one
of its two arguments (Section [5.5](#sec5-5) gives the equation for
U\_{s,k}). So M_s^{\mathrm{obs}} *appears inside* a formula indexed by
k, without itself depending on k. Read out of context, a formula like
“U\_{s,k} = \min(M_s^{\mathrm{obs}}, \dots)” can look as though
M_s^{\mathrm{obs}} has picked up a k subscript it does not actually
have. It has not. M_s^{\mathrm{obs}} is always the same single number
for every k within stratum s (the stratum’s own observed total, exactly
as defined in Section [2.4](#sec2-4)), and it only ever plays the role
of one candidate value inside a \min or \max whose *other* candidate
value depends on k (a sum of S_i over a subset of observations that
depends on k). Every formula in this document that mixes
M_s^{\mathrm{obs}} with a k-indexed quantity should be read this way:
one plain scalar meeting one per-coordinate quantity, not a hidden
second family of “M_s^{\mathrm{obs}} values.” Sections [5.5](#sec5-5)
and [7](#sec7) use the same pair of symbols, U\_{s,k} and L\_{s,k},
throughout for exactly this reason, rather than re-deriving and renaming
the same quantity elsewhere.

------------------------------------------------------------------------

## 3. The general principle: exponential families and eliminating nuisance parameters by conditioning

The Poisson-specific derivations in Sections [5](#sec5)-[7](#sec7) are
all examples of one classical fact about exponential families that holds
for any model, not just Poisson models. Stating that fact first, in
general terms, makes clear that “eliminate nuisance parameters by
conditioning on their sufficient statistics” is not a Poisson-specific
trick, and is not specific to any one of the three algorithms. All three
compute the *same* conditional distribution, just by different routes.
This section also fixes the vocabulary (“canonical parameter”,
“canonical sufficient statistic”) used for the rest of the document.

### 3.1 The model as a canonical exponential family

The Poisson log-likelihood for model (2.1) is \begin{aligned}
\ell(\tau,\beta,\boldsymbol{\gamma}) &= \sum\_{i=1}^n
\big(y_i\log\mu_i - \mu_i - \log(y_i!)\big) \\ &= \sum\_{i=1}^n
y_i\big(o_i+\tau\_{s(i)}+\beta z_i+\boldsymbol{\gamma}^\top
\mathbf{x}\_{N,i}\big) \\-\\ \sum\_{i=1}^n \mu_i \\-\\\sum\_{i=1}^n
\log(y_i!). \end{aligned} Collect the terms that are linear in y_i,
grouped by which parameter they multiply:
\ell(\tau,\beta,\boldsymbol{\gamma}) = \underbrace{\sum\_{s=1}^Q \tau_s
M_s(y)}\_{\text{stratum intercepts}} \\+\\ \underbrace{\beta\\
T(y)}\_{\text{interest}} \\+\\ \underbrace{\boldsymbol{\gamma}^\top
\mathbf{T}\_N(y)}\_{\text{nuisance covariates}} \\-\\ \sum_i \mu_i \\+\\
\log w(y), \qquad\text{(3.1)} Here M_s(y) := \sum\_{i\in I_s} y_i, T(y)
:= \sum_i y_i z_i, \mathbf{T}\_N(y) := X_N^\top y, and w(y) := \prod_i
E_i^{y_i}/y_i! is the Poisson base measure. (Section [5.1](#sec5-1)
spells out every algebraic step of this same factorization, for anyone
who wants to see it in full.) Equation (3.1) is exactly the canonical
form of a multi-parameter exponential family: a linear combination of
canonical parameters (\tau_1,\dots,\tau_Q,\beta,\gamma_1,\dots,\gamma_p)
against their canonical sufficient statistics
(M_1,\dots,M_Q,T,\mathbf{T}\_N), plus a normalizing term and a base
measure that do not depend on the parameters.

### 3.2 The conditioning theorem

Here is the general fact. Take a canonical exponential family, and split
its canonical parameter into a target piece and a nuisance piece. The
conditional distribution of the data, given the sufficient statistic
**of the nuisance piece alone**, no longer depends on the nuisance
parameters at all. It depends only on the target parameter, in the usual
tilted form. This is a standard result in the theory of exponential
families (Lehmann and Romano 2005, ch. 2 and ch. 4; Andersen 1970; Cox
1970; see also Agresti 2013, ch. 6-7, for its use with nuisance
intercepts in categorical/count data). Applied to (3.1), with
(\tau,\boldsymbol{\gamma}) as the nuisance piece and \beta as the
target, it says: \Pr\nolimits\_{\tau,\beta,\boldsymbol{\gamma}}\big(Y=y
\\ \big\|\\ M_s(Y)=M_s^{\mathrm{obs}}\\ \forall s,\\ \\
\mathbf{T}\_N(Y)=\mathbf{t}\_{N,\mathrm{obs}}\big) \\ \\ \text{does not
depend on } (\tau,\boldsymbol{\gamma}). \qquad\text{(3.2)} This is *why*
conditioning works in general, not only for this particular model and
not only for one of the three algorithms below. It follows directly from
(3.1) being in canonical exponential-family form. It does not need to be
re-checked by hand every time a new nuisance structure is added, or a
new algorithm is written to compute it. What *is* specific to the
Poisson distribution is that the resulting conditional distribution
turns out to have a convenient, explicitly computable form: a
base-measure-weighted count over a combinatorial reference set (Section
[3.4](#sec3-4)). That convenience is what the rest of this document uses
computationally. But the validity of conditioning itself follows from
(3.1) alone, with nothing Poisson-specific about it.

### 3.3 Finite conditional reference sets

Poisson counts are unbounded, so turning (3.2) into something we can
actually sum over requires a finite sample space, made explicit.
Introduce per-observation caps S_i \in \mathbb{Z}\_{\ge 0}, with 0 \le
y_i \le S_i (Section [5.2](#sec5-2) gives the exact rules used to set
them, together with a short proof of which are forced by conditioning
alone and which are extra modelling choices; see in particular the
paragraph beginning “Rule 2’s cap is implied by conditioning alone,”
which derives that claim directly from the stratum-total conditioning
event itself, with nothing further assumed). Define the **cap box**
\mathcal{B}(Q) := \\\\ y \in \mathbb{Z}\_{\ge0}^n : 0 \le y_i \le S_i\\
\forall i \\\\. The **intercept-conditioned feasible set** and the
**fully conditioned set** are \Omega := \Big\\\\ y \in \mathcal{B}(Q) :
\sum\_{i \in I_s} y_i = M_s^{\mathrm{obs}}\\ \\ \forall s \\\Big\\,
\qquad \Omega\_{\mathrm{cond}} := \\\\ y \in \Omega : X_N^\top y =
\mathbf{t}\_{N,\mathrm{obs}} \\\\. For any statistic value t and
nuisance value \mathbf{c}, the **fibers** (the sets of allocations that
achieve exactly those values) are \Omega(t,\mathbf{c}) := \\\\ y \in
\Omega : T(y)=t,\\ X_N^\top y = \mathbf{c} \\\\, \qquad
\Omega\_{\mathrm{cond}}(t) := \\\\ y \in \Omega\_{\mathrm{cond}} : T(y)
= t \\\\.

### 3.4 Coefficient tables and the conditional law of T

This section derives (3.3) below step by step from (3.1)-(3.2), rather
than only stating the result, because “collecting like terms into a
coefficient table” is exactly what every algorithm in Sections
[5](#sec5)-[7](#sec7) computes. They differ only in the bookkeeping
strategy used to do it, never in the target itself.

**Step 1: why \tau and \boldsymbol{\gamma} drop out of the conditional
probability, concretely.** By definition, conditioning restricts
attention to y \in \Omega\_{\mathrm{cond}} and renormalizes:
\begin{aligned} \Pr\nolimits\_{\tau,\beta,\boldsymbol{\gamma}}(Y=y \mid
Y \in \Omega\_{\mathrm{cond}}) =
\frac{\Pr\nolimits\_{\tau,\beta,\boldsymbol{\gamma}}(Y=y)}{\sum\_{y' \in
\Omega\_{\mathrm{cond}}}
\Pr\nolimits\_{\tau,\beta,\boldsymbol{\gamma}}(Y=y')}, \qquad y \in
\Omega\_{\mathrm{cond}}. \end{aligned} From the canonical form (3.1),
\Pr\_{\tau,\beta,\boldsymbol{\gamma}}(Y=y) \propto \exp\\\big(\sum_s
\tau_s M_s(y) + \beta T(y) + \boldsymbol{\gamma}^\top
\mathbf{T}\_N(y)\big)\\ w(y). The term -\sum_i \mu_i from (3.1) has
already been absorbed into this “\propto,” and it is worth saying
explicitly why, since the reason is more basic than the argument given
below for \tau and \boldsymbol{\gamma}. Unlike M_s(y), T(y), and
\mathbf{T}\_N(y), the quantity \mu_i = E_i\exp(\tau\_{s(i)}+\beta
z_i+\boldsymbol{\gamma}^\top \mathbf{x}\_{N,i}) does not depend on y at
all; it is a function only of the parameters and the fixed design and
offsets. So once (\tau,\beta,\boldsymbol{\gamma}) are held fixed, \sum_i
\mu_i is a single number, the same for *every* candidate y in the entire
sample space, not merely for y \in \Omega\_{\mathrm{cond}}. The factor
e^{-\sum_i\mu_i} therefore multiplies \Pr(Y=y) by an identical constant
no matter which y is plugged in, so it cancels top and bottom in any
ratio of probabilities and can be dropped the instant two probabilities
are compared, which is exactly what “\propto” records here. (In
exponential-family language, \sum_i\mu_i is the log-partition function
of the canonical form (3.1): the piece that depends on the parameters
alone and never on the data. It always drops out of statements about the
*shape* of the distribution across y, with or without any conditioning.)
Now look at the \tau and \boldsymbol{\gamma} terms for y ranging over
\Omega\_{\mathrm{cond}}. By the definition of \Omega\_{\mathrm{cond}}
(Section [3.3](#sec3-3)), *every* y in this set satisfies M_s(y) =
M_s^{\mathrm{obs}} for all s, and \mathbf{T}\_N(y) =
\mathbf{t}\_{N,\mathrm{obs}}. So \sum_s \tau_s M_s(y) = \sum_s \tau_s
M_s^{\mathrm{obs}} and \boldsymbol{\gamma}^\top \mathbf{T}\_N(y) =
\boldsymbol{\gamma}^\top \mathbf{t}\_{N,\mathrm{obs}} take the *same*
value for every single y in \Omega\_{\mathrm{cond}}. They never
distinguish one member of \Omega\_{\mathrm{cond}} from another. A factor
that is constant across every term of both the numerator and the
denominator of a ratio cancels exactly: \begin{aligned}
\Pr\nolimits\_\beta(Y=y \mid Y \in \Omega\_{\mathrm{cond}}) &=
\frac{e^{\sum_s \tau_s M_s^{\mathrm{obs}} + \beta T(y) +
\boldsymbol{\gamma}^\top\mathbf{t}\_{N,\mathrm{obs}}}\\
w(y)}{\sum\_{y'\in\Omega\_{\mathrm{cond}}} e^{\sum_s \tau_s
M_s^{\mathrm{obs}} + \beta T(y') +
\boldsymbol{\gamma}^\top\mathbf{t}\_{N,\mathrm{obs}}}\\ w(y')} \\ &=
\frac{e^{\beta T(y)}\\ w(y)}{\sum\_{y'\in\Omega\_{\mathrm{cond}}}
e^{\beta T(y')}\\ w(y')} \\ \propto\\ w(y)\\ e^{\beta T(y)}.
\end{aligned} \qquad\text{(3.3a)} This makes (3.2) concrete rather than
merely asserted: \tau and \boldsymbol{\gamma} vanish because
conditioning pins their canonical sufficient statistics to one exact
value, shared by every allocation still under consideration. Only \beta
survives in the conditional distribution, because its statistic T(y)
genuinely differs from one y \in \Omega\_{\mathrm{cond}} to another.

**Step 2: grouping allocations by the value of T.** We ultimately want
the distribution of the *scalar* T(Y), not the distribution of the full
vector Y. Many different allocations y can share the same value of T(y),
so partition \Omega\_{\mathrm{cond}} into **fibers**
\Omega\_{\mathrm{cond}}(t): the set of allocations that all give the
*same* value T(y)=t (Section [3.3](#sec3-3)). Inside one fiber, the
factor e^{\beta T(y)} = e^{\beta t} in (3.3a) is identical for every
member, so it can be pulled out of the sum over that fiber:
\Pr\nolimits\_\beta(T=t \mid Y\in\Omega\_{\mathrm{cond}}) = \sum\_{y \in
\Omega\_{\mathrm{cond}}(t)} \Pr\nolimits\_\beta(Y=y\mid
Y\in\Omega\_{\mathrm{cond}}) \\ \propto\\ e^{\beta t} \sum\_{y \in
\Omega\_{\mathrm{cond}}(t)} w(y). Define the **coefficient table** as
exactly that inner sum: add up the base-measure weight of every
allocation that produces this particular value of t: C(t) := \sum\_{y
\in \Omega\_{\mathrm{cond}}(t)} w(y). For example, suppose t=3 happens
to be achievable by exactly two allocations, y^{(1)} and y^{(2)}, with
base-measure weights w(y^{(1)}) and w(y^{(2)}). Then C(3) = w(y^{(1)}) +
w(y^{(2)}): nothing more than collecting like terms. C(t) is precisely
the coefficient of e^{\beta t} once every allocation’s contribution
w(y)e^{\beta T(y)} has been grouped by its exponent t; see Section
[6.3](#sec6-3) for further details, where this same idea is made
completely literal: the closed-form algorithm’s own coefficient table
is, quite concretely, the coefficient of a power of a single polynomial
variable.

**Step 3: the joint table C(t,\mathbf{c}), and why it gets computed
first.** \Omega\_{\mathrm{cond}} additionally requires \mathbf{T}\_N(y)
= \mathbf{t}\_{N,\mathrm{obs}} exactly. This constraint is awkward to
enumerate directly, since it depends on which nuisance total actually
happened to be observed. It is far more convenient to enumerate over the
*larger*, easier-to-describe set \Omega instead (only the stratum-total
constraints, no nuisance constraint yet; Section [3.3](#sec3-3)),
grouping allocations by **both** T(y)=t and \mathbf{T}\_N(y)=\mathbf{c}
at once: C(t,\mathbf{c}) := \sum\_{y \in \Omega(t,\mathbf{c})} w(y).
C(t) is then just one slice of this finer table, at the observed
nuisance value: C(t) = C(t,\mathbf{t}\_{N,\mathrm{obs}}), since
\Omega\_{\mathrm{cond}}(t) = \Omega(t,\mathbf{t}\_{N,\mathrm{obs}})
follows directly from the definition of \Omega\_{\mathrm{cond}}. This is
exactly why Sections [5](#sec5) and [7](#sec7) build the
two-dimensional-in-(t,\mathbf{c}) table first (by dynamic programming
and convolution, or by a network with continuous pruning), and only
*afterwards* filter down to \mathbf{c}=\mathbf{t}\_{N,\mathrm{obs}} and
drop the \mathbf{c} key (Section [5.6](#sec5-6), “filtering and
marginalizing”). Both algorithms enumerate over the convenient set
\Omega the whole time, and only narrow down to \Omega\_{\mathrm{cond}}
at the very last step (or, for the network algorithm, gradually before
it).

**Step 4: the conditional pmf of T.** Normalizing Step 2’s
proportionality (dividing by the same sum over every attainable t) gives
the exponential-tilt form used throughout the rest of this document:
\begin{aligned} p\_\beta(t) \\=\\ \Pr\nolimits\_\beta\big(T=t \mid
\mathbf{T}\_N = \mathbf{t}\_{N,\mathrm{obs}}, \\M_s =
M_s^{\mathrm{obs}}\\\big) \\&=\\ \frac{C(t)\\ e^{\beta
t}}{\displaystyle\sum\_{u \in \mathcal{T}} C(u)\\ e^{\beta u}}, \\
\qquad \mathcal{T} :&= \\\\t : \Omega\_{\mathrm{cond}}(t) \ne
\varnothing\\\\. \end{aligned} \qquad\text{(3.3)} Every downstream exact
quantity (p-value, confidence interval, conditional maximum likelihood
estimate, **MLE**) is a function of the one-dimensional table \\(t,
C(t))\\\_{t \in \mathcal{T}} alone (Section [8](#sec8)). (“Conditional”
MLE because it maximizes the conditional likelihood p\_\beta(t) of (3.3)
just derived, not the ordinary, unconditional likelihood of Y; Section
[8.2](#sec8-2) gives the precise equation it solves.) **The entire
computational problem of exact conditional Poisson regression is
therefore the problem of computing C(t)** (equivalently, computing
C(t,\mathbf{c}) and then slicing and marginalizing it, as in Step 3).
The three algorithms in this package differ only in *how* they compute
this table, never in what it means or in what is done with it once it is
obtained.

------------------------------------------------------------------------

## 4. Adjustment vs. stratification: when they agree, and when they don’t

Section [2.4](#sec2-4) treats “strata” and “nuisance covariates to
adjust for” as two visibly different mechanisms. Strata get their own
free intercept \tau_s per level, with no dummy coding and no reference
level. Adjustment covariates enter through \boldsymbol{\gamma}^\top
\mathbf{x}\_{N,i}, with the usual reference-coded dummy columns. A
natural question is whether these are truly different things, or the
same thing computed two different ways. It is worth answering precisely,
since the package’s API (the `strata=` argument versus adjusting for a
variable in the formula) forces a choice, for every algorithm in this
package equally.

The short answer: they give identical exact inference in one specific,
narrow case (Sections [4.1](#sec4-1)-[4.2](#sec4-2)), and that case
generalizes cleanly to several variables stratified jointly, *provided*
the adjustment side matches it term for term (Section [4.3](#sec4-3)).
Trouble starts when it doesn’t: when a variable is added to an existing
stratification as a plain adjustment covariate rather than as a further
stratum. Section [4.4](#sec4-4) walks through exactly why that produces
different numbers, using a worked example. Section [4.7](#sec4-7) then
gives practical guidance, in a clinical-trials setting, for choosing
between the two.

### 4.1 The narrow case where the two are the same model

Suppose G is a categorical variable with levels g=1,\dots,Q, and it is
the *only* nuisance structure being considered, with nothing else in the
model yet. Let its levels define the stratum partition, I_s := \\i : G_i
= s\\. Compare:

- **Model A (stratify by G alone).** \eta_i = \tau\_{s(i)} + \beta z_i,
  with Q free stratum intercepts \tau_1,\dots,\tau_Q and no shared
  global intercept.
- **Model B (adjust for G alone).** \eta_i = \tau_0 + \beta z_i +
  \delta^\top d_i, with a single global intercept \tau_0, and d_i \in
  \\0,1\\^{Q-1} the usual reference-coded dummy vector for G (level 1
  dropped as reference), with free coefficient vector \delta \in
  \mathbb{R}^{Q-1}.

These are the same statistical model, just written with different labels
for the same parameters. Setting \tau_0 := \tau_1 and \delta_g :=
\tau_g - \tau_1 for g=2,\dots,Q gives \eta_i^{(A)} = \eta_i^{(B)} for
every observation and every value of the free parameters. The map
(\tau_1,\dots,\tau_Q) \leftrightarrow (\tau_0,\delta_2,\dots,\delta_Q)
is an invertible linear change of variables from \mathbb{R}^Q to
\mathbb{R}^Q. So Model A and Model B produce exactly the same set of
possible linear predictors, and therefore exactly the same likelihood as
a function of the data.

### 4.2 The two conditioning events are the same event

Model A’s canonical sufficient statistics for its nuisance parameters
(\tau_1,\dots,\tau_Q) are the stratum totals (M_1,\dots,M_Q) (Section
[3.1](#sec3-1)). Model B’s canonical sufficient statistics for its
nuisance parameters (\tau_0,\delta_2,\dots,\delta_Q) are
\big(M\_{\mathrm{total}},\\ T\_{N,2},\dots,T\_{N,Q}\big), where
M\_{\mathrm{total}} := \sum\_{i=1}^n Y_i and T\_{N,g} := \sum_i Y_i\\
d\_{i,g} = \sum\_{i \in I_g} Y_i = M_g (since d\_{i,g}=1 exactly when
G_i=g). So Model B’s statistic vector is (M\_{\mathrm{total}},
M_2,\dots,M_Q) = \big(\textstyle\sum\_{s=1}^Q M_s,\\ M_2,\dots,M_Q\big):
an **invertible linear function** of Model A’s (M_1,\dots,M_Q), since
M_1 = M\_{\mathrm{total}} - \sum\_{g=2}^Q M_g recovers it.

Conditioning on a statistic, and conditioning on any invertible function
of that statistic, are the same operation: they condition on the same
event. If f is injective on the attainable values of (M_1,\dots,M_Q),
then \\(M_1,\dots,M_Q) = \mathrm{obs}\\ and \\f(M_1,\dots,M_Q) =
f(\mathrm{obs})\\ are literally the same set. So the conditional
distribution of anything is identical whichever way it is written
(Lehmann and Romano 2005, ch. 2; this is the same reasoning that lets
one condition on “the sample total” or “the sample mean”
interchangeably). So: \begin{aligned} \Pr\nolimits\_\beta\big(T&=t \mid
M_1=M_1^{\mathrm{obs}},\dots,M_Q=M_Q^{\mathrm{obs}}\big) \\=\\\\
\Pr\nolimits\_\beta\big(T&=t \mid
M\_{\mathrm{total}}=M\_{\mathrm{total}}^{\mathrm{obs}},\\
M_2=M_2^{\mathrm{obs}},\dots,M_Q=M_Q^{\mathrm{obs}}\big). \end{aligned}
\qquad\text{(4.1)} The left side is what any algorithm in this package
computes when G is the only thing passed as `strata=`. The right side is
what it would compute if G were instead the only thing adjusted for.
**In this narrow case, the two are provably the same conditional
distribution**: the same exact p-value, the same exact confidence
interval, the same conditional MLE.

### 4.3 The general principle: crossing variables into a stratification means adjusting for their full joint interaction (precisely defined)

Section [4.1](#sec4-1) covered a single variable G with nothing else in
the model. This section states the general version, for any number of
variables considered jointly, and fixes precisely what “the variables”
and “the full set of interaction dummy columns among them” mean, as
promised at the start of Section [4](#sec4). Every symbol below is
defined before it is used.

**Setup and notation.** Let V_1,\dots,V_k denote *any* set of
categorical variables that are being stratified **jointly**: that is,
`strata = c("V1", ..., "Vk")` in the package. This works whichever way
the variables arose. They might be a brand-new set of stratifiers
introduced together (for example, stratifying by age, gender, and region
all at once, in a model that had no strata before), or an existing set
of strata (V_1,\dots,V\_{k-1}, already in the model) with one further
variable V_k now being added to them. The statement below does not care
which; it is about the *set* \\V_1,\dots,V_k\\ taken together as the
thing being crossed into one joint stratification.

Define the **composite (crossed) variable** G^\star as the categorical
variable whose levels are every combination of (V_1,\dots,V_k) that
actually occurs in the data. If V_1,\dots,V_k have L_1,\dots,L_k levels
respectively, G^\star has up to L_1\times\cdots\times L_k levels (fewer
if some combinations are empty). “Stratify jointly by V_1,\dots,V_k”
means, exactly, *stratify by the single variable G^\star*: one free
intercept per combination of (V_1,\dots,V_k) that occurs, with no
relationship of any kind assumed between the intercepts of different
combinations.

**The precise statement.** Apply Section [4.1](#sec4-1)’s Model A/Model
B argument with G := G^\star. Model A is “stratify by G^\star” (one free
intercept per combination: this is exactly `strata = c("V1",...,"Vk")`).
Model B is “adjust for G^\star alone,” meaning: one shared global
intercept, plus a reference-coded dummy column for every combination of
(V_1,\dots,V_k) that occurs, except one reference combination. Sections
[4.1](#sec4-1)-[4.2](#sec4-2) apply without any change (they never used
anything about G except that it is a single categorical variable), so:

> **Stratifying jointly by V_1,\dots,V_k is an exact reparametrisation
> of adjusting for one global intercept plus the full set of dummy
> columns for the crossed variable G^\star = V_1\times\cdots\times
> V_k.**

This is what “the full set of interaction dummy columns among those
variables” means, precisely: dummy-coding the single composite variable
G^\star built by crossing *all* of V_1,\dots,V_k together, referenced
against one baseline combination, not each V_j’s own main-effect dummies
considered separately. A single set of G^\star dummy columns
automatically carries every V_j’s main effect *and* every interaction
among them (all pairwise interactions, all three-way interactions, and
so on up to the full k-way interaction) as one indivisible block. There
is no way to keep some of that structure and drop the rest and still
have the same reparametrisation. Section [4.1](#sec4-1) is the k=1 case
of this statement, where there is only one variable and so no
interaction terms to speak of. That is exactly why the distinction
between “main effect” and “full interaction” never came up there.

### 4.4 Worked contrast: why `fit23` and `fit232` disagree

Section [4.3](#sec4-3) makes it possible to say exactly why the two
calls below produce different results, and to measure the difference.
Take V_1= age (4 levels) and V_2= gender (2 levels) as the pre-existing
strata (so R = 4\times 2 = 8 age\timesgender combinations), and V_3=
region as the variable being added (Q=4 levels, so Q-1=3 main-effect
dummies with region 1 as reference).

``` r
fit23 <- ep_multi_shift(
  formula = cases ~ group + region, strata = c("age","gender"),
  exact_vars = "group", ...
)   # region as a main-effect adjustment, nested in 8 age x gender strata

fit232 <- ep_multi_shift(
  formula = cases ~ group, strata = c("age","gender","region"),
  exact_vars = "group", ...
)   # age, gender, and region stratified jointly: 32 age x gender x region strata
```

`fit232` stratifies jointly by V_1,V_2,V_3 together. By Section
[4.3](#sec4-3) this is the exact reparametrisation of one global
intercept plus the full set of dummy columns for G^\star =
\text{age}\times\text{gender}\times\text{region}: that is, R\times Q - 1
= 32-1=31 dummy columns, encoding the age main effect, the gender main
effect, the region main effect, the age\timesgender interaction, the
age\timesregion interaction, the gender\timesregion interaction, *and*
the age\timesgender\timesregion three-way interaction, all at once.
Equivalently, as `fit232` actually computes it: 32 freely-varying
stratum intercepts, one per combination, with \Omega\_{\mathrm{cond}}
built by conditioning on all 32 individual cell totals.

`fit23`, by contrast, has only R + (Q-1) = 8 + 3 = 11 free nuisance
parameters: the 8 existing age\timesgender intercepts (which fully
account for the age main effect, the gender main effect, and their
mutual interaction, per Section [4.1](#sec4-1) applied to the composite
age\timesgender variable), plus region’s 3 main-effect dummies, held
*fixed across all 8 strata*. That last part is the restriction: `fit23`
assumes region’s effect on the log-rate is exactly the same number in
every one of the 8 age\timesgender strata. It contains no
age\timesregion, gender\timesregion, or age\timesgender\timesregion
terms at all.

Comparing counts directly: 11 free parameters in `fit23` versus 32 in
`fit232`, a difference of 21 parameters. That number is not arbitrary.
It is exactly (R-1)(Q-1) = 7\times 3 = 21, the standard degrees of
freedom for the interaction term in an R\times Q factorial layout.
`fit23` is `fit232` with precisely that (R-1)(Q-1)-dimensional
interaction block assumed away. So `fit23` is a genuine sub-model of
`fit232`, a real, testable statistical restriction, not just a
relabelling. The equivalence of Sections [4.1](#sec4-1)-[4.2](#sec4-2)
requires matching parameter counts on both sides, and that simply does
not hold between `fit23` and `fit232`.

The same conclusion follows directly from the sufficient statistics.
`fit23` conditions on 8 age-by-gender totals plus 3 *global* region
totals (summed across all 8 strata): 11 numbers, matching its 11 free
parameters. `fit232` conditions on all 32 individual
age-by-gender-by-region cell totals. The 32 cell totals determine the 11
coarser numbers (each is a sum of some of the 32 cells), but the reverse
is false, exactly as row and column totals of a contingency table do not
determine its individual cell counts. `fit232` therefore conditions on
strictly more information, and conditioning on more detailed information
can, and generally will, shift the resulting conditional distribution of
`group`’s sufficient statistic. That is exactly why the two commands
give different answers. This is not a failure of Sections
[4.1](#sec4-1)-[4.2](#sec4-2): those sections proved an equivalence for
a case (a single variable, or a full joint crossing matched by full
joint adjustment) that simply does not cover `fit23` versus `fit232`,
since `fit23` never includes the 21 interaction dummy columns that would
have made it match.

### 4.5 Other conditions the equivalence needs

Even in the settings covered by Sections [4.1](#sec4-1) and
[4.3](#sec4-3), a few further conditions are worth stating explicitly:

- The composite variable’s levels (G in Section [4.1](#sec4-1), G^\star
  in Section [4.3](#sec4-3)) must *be* the stratum partition exactly. A
  continuous covariate that has been grouped into bins is different:
  adjusting for the bins conditions on the bin totals, which is a
  different sufficient statistic from adjusting for the covariate on its
  original continuous scale (a mode of adjustment this package does not
  implement at all; see Section [2.3](#sec2-3)).
- The layout must be saturated: one free parameter per level (Section
  [4.1](#sec4-1)) or per combination (Section [4.3](#sec4-3)), with
  nothing shared or borrowed across them. A model that shares
  information across levels (smoothed, penalized, or otherwise
  constrained) is a different model, with different sufficient
  statistics, and the equivalence does not apply to it.

Within these conditions, the equivalence is exact, not approximate. It
is a direct consequence of sufficiency being unchanged by invertible
reparametrisation (Section [4.2](#sec4-2)), with no large-sample
approximation anywhere in the argument.

### 4.6 Why the package still treats stratification and adjustment as different mechanisms

Given Sections [4.1](#sec4-1)-[4.2](#sec4-2), one might ask why the
package bothers to offer `strata=` as a separate argument at all, rather
than having users dummy-code everything into the adjustment side of the
formula. The reason is computational cost, not statistical validity. The
*answer* can be identical (in the narrow case where it applies), but the
*cost* of getting there is not, and the difference is large enough to
matter, for every algorithm in this package.

**Stratification is enforced locally, at no cost to the size of the
tracked table.** In the dynamic program shared by Sections [5](#sec5)
and [7](#sec7), the constraint \sum\_{i\in I_s} y_i = M_s^{\mathrm{obs}}
is enforced simply by keeping only the bucket m=M_s^{\mathrm{obs}} at
the end of stratum s’s own pass. This is bookkeeping entirely internal
to that one stratum’s loop, and it adds nothing to the (t,\mathbf{c})
key carried through the cross-stratum combination step (Section
[5.4](#sec5-4)). Whether there are 2 strata or 2,000, stratification
alone never changes the dimensionality of the joint coefficient table’s
key. It only changes how many small per-stratum contributions need to be
combined.

**Adjustment is enforced globally, and adds a dimension to every tracked
key.** If the same Q-level variable were instead folded into X_N as Q-1
dummy columns, the adjustment machinery would need to carry Q-1 extra
axes through *every* stratum’s computation and through the *entire*
cross-stratum combination step, applying Rules A and B (Section
[5.5](#sec5-5)) on all of them, and only filtering down to one observed
combination at the very end (Section [5.6](#sec5-6)). This is exactly
the “extra nuisance dimension” cost discussed in Section [10](#sec10).

**In practice, this means strata can be far finer than adjustment
covariates ever could be.** Because stratification is free in table
dimensionality, every algorithm in this package handles a design
stratified down to one observation per treatment arm per stratum (an
individually matched design, with almost as many strata as observations)
just as easily as it handles a handful of coarse strata. The per-stratum
computation is trivial, and combining strata just means combining more
(very small) contributions. Encoding that same fine partition as
adjustment dummy columns instead would be computationally hopeless,
since the coefficient table would need an axis for nearly every
observation. This is why exact conditional regression methods offer
stratification as a distinct mechanism from covariate adjustment (Hirji,
Mehta and Patel 1987; Mehta and Patel 1995), and it is why SAS’s
`PROC GENMOD EXACT` statement, whose cap and elimination rules this
package follows (Section [5.2](#sec5-2)), likewise treats a `STRATA`
statement as separate from the `CLASS`/model-effect covariates entering
the `EXACT` list.

### 4.7 Choosing between stratification and adjustment in practice: a clinical-trials perspective

Sections [4.1](#sec4-1)-[4.6](#sec4-6) settle the *mathematics*: when
the two mechanisms coincide, when they don’t, and what they cost
computationally. They do not say which one a given covariate *should*
get in a given analysis. That is a statistical modelling choice, not a
fact about the algorithm, and it deserves its own guidance.

**What each mechanism assumes.** Stratifying a covariate gives it a
fully free, unstructured effect: one parameter per level (or, jointly,
per combination, Section [4.3](#sec4-3)), with nothing assumed about how
levels relate to each other, and that parameter is eliminated entirely
rather than estimated. Choose this when the covariate’s own effect is
not of interest and there is no wish to assume anything about its shape.
Adjusting for a covariate instead imposes structure: typically a single
coefficient, or a small shared set of dummy coefficients, applied
identically everywhere else in the model. This borrows information
across levels and estimates the covariate’s own effect, but only under
the assumption that a shared, simple form actually holds (linear on the
log-rate scale, or one constant multiplicative shift per level, the same
in every existing stratum).

**Randomization stratification factors.** Variables that were literally
used to stratify the randomization (disease stage, prior therapy,
region, sometimes site) are conventionally carried into the primary
analysis as strata, not as adjustment covariates. This reflects the
trial’s design, not just computational convenience: regulatory guidance
on trial statistics (for example, ICH E9(R1)) ties the analysis to the
randomization scheme it was built for. Treating a randomization factor
as a stratum lets its baseline rate differ arbitrarily by level, with no
assumption about how, say, region A relates to region B. This is
appropriate when standard-of-care or population differences are
expected, but there is no wish to model that relationship explicitly.

**Continuous or many-level prognostic covariates.** Baseline severity
scores, biomarker levels, age, and similar variables are, in principle,
best handled by adjusting for them as covariates on their original
continuous scale rather than by stratifying. Binning a continuous
variable into strata discards information, and in a trial of modest size
it creates many strata with few events each, some contributing almost
nothing once conditioned on. As noted in Section [2.3](#sec2-3),
however, this package’s exact-conditional implementation only supports
categorical effects, so a continuous covariate must first be grouped
into a small number of ordered categories before it can be stratified on
or adjusted for. That grouping step reintroduces the information loss
just described, and is a real, separate limitation to weigh against the
assumptions discussed here, not a property of stratification versus
adjustment as such.

Once a continuous covariate has been binned, though, it is just an
ordinary categorical variable, and the choice of what to do with those
bins is governed by exactly the same rules already given in this
section, with no special case for the fact that the categories
originated from a continuous scale: Sections
[4.1](#sec4-1)-[4.3](#sec4-3)’s stratify/adjust equivalence, Section
[4.6](#sec4-6)’s cost trade-off, and the bias-variance and interaction
guidance above all apply to the binned variable directly. To be concrete
about what this package recommends in practice: prefer as few bins as
the clinical question allows: established, pre-specified clinical
cutpoints if they exist (for example, a standard severity-score
threshold), or otherwise tertiles or quartiles of the observed data,
rather than many narrow bins, since each additional bin is either an
extra stratum (more sparsely-populated strata, Section [4.6](#sec4-6))
or an extra adjustment dummy (more nuisance columns, at the cost
described in Section [4.6](#sec4-6)). Adjust for the binned covariate
(rather than stratifying on it) unless there is a specific reason to
think its association with the outcome differs across levels of the
effect under test, in which case stratifying jointly with it, per
Section [4.4](#sec4-4)’s worked contrast, is the honest choice.

Stratifying costs nothing in table dimensionality in this package’s
exact-conditional algorithms (Section [4.6](#sec4-6)), but it still
costs information statistically: each stratum’s contribution to
inference comes only from the observations inside it. Many strata with
sparse events mean a lot of unshared, individually estimated structure.
Adjustment assumes commonality across strata up front, which is more
efficient when that assumption holds and biased when it does not.

Section [4.4](#sec4-4)’s worked contrast is the sharpest version of this
trade-off. Stratifying a covariate jointly with existing strata
(`fit232`-style) lets its effect vary completely freely across every
existing stratum. Adjusting for it as a covariate nested inside existing
strata (`fit23`-style) assumes its effect is the same constant
everywhere. If there is reason to think a covariate’s association with
the outcome genuinely differs across levels of another factor already in
the model (a treatment-by-region interaction, say), stratifying jointly
(or explicitly fitting interaction terms) is the honest choice. Treating
it as a plain, main-effect-only adjustment covariate is only valid to
the extent that interaction can reasonably be ruled out, and Section
[4.4](#sec4-4) shows concretely how large a restriction that assumption
is: (R-1)(Q-1) degrees of freedom of interaction structure, assumed
away.

------------------------------------------------------------------------

## 5. The multivariate-shift algorithm (`ep_multi_shift()`)

This algorithm computes C(t) for any number of categorical covariates
(any number of nuisance columns p \ge 0) and any number of strata,
following the conventions of SAS `PROC GENMOD`’s `EXACT` statement. It
is the **multivariate-shift algorithm**: a per-stratum dynamic program
combined with explicit cross-stratum convolution. It adapts ideas from
Hirji, Mehta and Patel (1987), originally built for exact logistic
regression where every y_i \in \\0,1\\, to Poisson counts bounded by
SAS-style caps. (This is distinct from SAS’s own `METHOD=NETWORK`,
described in Section [7](#sec7), which prunes continuously, observation
by observation, rather than building an explicit joint table per
stratum.)

### 5.1 The Poisson base measure, explicitly

Section [3.1](#sec3-1) stated this compactly; here is the full
derivation. The Poisson pmf \Pr(Y_i=y_i\mid\mu_i) =
\mu_i^{y_i}e^{-\mu_i}/y_i! gives, using (2.1), \prod_i
\Pr(Y_i=y_i\mid\mu_i) = \prod_i \frac{\big(E_i\\e^{\tau\_{s(i)}+\beta
z_i+\boldsymbol{\gamma}^\top \mathbf{x}\_{N,i}}\big)^{y_i}}{y_i!}\\
e^{-E_i e^{\tau\_{s(i)}+\beta z_i+\boldsymbol{\gamma}^\top
\mathbf{x}\_{N,i}}}. Every factor e^{y_i(\tau\_{s(i)}+\beta z_i +
\boldsymbol{\gamma}^\top \mathbf{x}\_{N,i})} contributes to the
canonical linear form (3.1). Everything else (E_i^{y_i}/y_i! and the
normalizing exponential) either does not depend on the parameters, or is
constant given the conditioning statistics. Collecting the
parameter-free part defines the **base measure** w(y) := \prod_i
E_i^{y_i}/y_i!, with \log w(y) = \sum_i(y_i o_i - \log(y_i!)): exactly
the quantity accumulated by the dynamic program below.

### 5.2 SAS-style cap rules for S_i

Three mutually exclusive rules set the per-observation cap S_i, matching
SAS’s `OFFSET`/`EXACTMAX` conventions exactly:

1.  **`EXACTMAX = var` specified.** S_i = e_i, the observed value of the
    supplied variable, per observation.
2.  **`EXACTMAX` specified without a variable, or neither `OFFSET` nor
    `EXACTMAX` specified.** There is nothing else to bound y_i by, so
    the procedure uses the intercept’s own sufficient statistic within
    the observation’s stratum: S_i = M\_{s(i)}^{\mathrm{obs}} (single
    stratum: S_i = M^{\mathrm{obs}} for everyone).
3.  **`OFFSET` specified, `EXACTMAX` not specified.** S_i =
    \max\\\big(\lceil e^{o_i}\rceil,\\ 1\big), with the feasibility
    requirement S_i \ge y_i^{\mathrm{obs}}.

These three rules are not all the same kind of restriction, and it is
worth being precise about that. **Rule 2’s cap is implied by
conditioning alone; it is not an extra assumption.** If \sum\_{j \in
I_s} y_j = M_s^{\mathrm{obs}} and every y_j \ge 0, then for any single i
\in I_s, y_i = M_s^{\mathrm{obs}} - \sum\_{j\ne i} y_j \le
M_s^{\mathrm{obs}} automatically. The stratum-total conditioning event
already forces 0 \le y_i \le M_s^{\mathrm{obs}} with no further
assumption. So Rule 2’s cap simply restates a bound that was already
there. (Concretely, the per-observation loop bound in Section
[5.3](#sec5-3), \min(S_i,\\ M_s^{\mathrm{obs}}-m), reduces to
M_s^{\mathrm{obs}}-m under Rule 2, since S_i = M_s^{\mathrm{obs}} \ge
M_s^{\mathrm{obs}}-m whenever m\ge0.) **Rules 1 and 3, by contrast, are
genuine additional modelling assumptions**, adopted for tractability and
to match SAS’s own convention (SAS Institute, *SAS/STAT User’s Guide*,
the `GENMOD` procedure’s `EXACT` statement). They assert that individual
counts above a stated bound (a known capacity e_i, or a bound implied by
the observation’s own exposure) are a priori impossible. Conditioning
does not give that assumption for free; it would need to be justified on
subject-matter grounds in any given application. First rule only
restates what was already true, second and third rules genuinely narrow
the reference set further. That is exactly why the cap rules are
documented as rules the package resolves automatically from what was
supplied, not a single automatic default and not something a user tunes
directly.

These three rules, and only these three, determine \mathcal{B}(Q) and
hence \Omega and \Omega\_{\mathrm{cond}} (Section [3.3](#sec3-3)). In
the code, this logic lives in `parse_exactmax_spec()`,
`si_from_exactmax_var()`, `si_from_stratum_total()`, and
`si_from_offset()`. The choice among the three rules is made once,
automatically, from whether `exactmax`/`offset` were supplied. There is
no argument for substituting a different rule.

### 5.3 Stratum-level decomposition and the within-stratum dynamic program

Strata are disjoint, and both the base measure and the sufficient
statistics factor (or add) across strata, so the problem decomposes
stratum by stratum. Write y_s = (y_i)\_{i \in I_s} and, restricted to
stratum s, define w_s(y_s) := \prod\_{i \in I_s} \frac{E_i^{y_i}}{y_i!},
\qquad T_s(y_s) := \sum\_{i \in I_s} y_i z_i, \qquad
\mathbf{T}\_{N,s}(y_s) := \sum\_{i \in I_s} y_i\\ \mathbf{x}\_{N,i},
\Omega_s := \Big\\\\ y_s \in \mathbb{Z}\_{\ge0}^{\|I_s\|} : \sum\_{i \in
I_s} y_i = M_s^{\mathrm{obs}}, \\ 0 \le y_i \le S_i \\\Big\\, C_s(t_s,
\mathbf{c}\_s \mid M_s^{\mathrm{obs}}) := \sum\_{y_s \in \Omega_s}
\mathbf{1}\\T_s(y_s)=t_s,\\ \mathbf{T}\_{N,s}(y_s)=\mathbf{c}\_s\\\\
w_s(y_s). Each stratum’s coefficient table C_s is built by a dynamic
program that processes that stratum’s observations one at a time,
keeping track of:

- **buckets**, indexed by the partial stratum total m = 0, 1, \dots,
  M_s^{\mathrm{obs}} accumulated so far;
- within each bucket, a (sparse) map from the key (t, \mathbf{c}) (the
  partial values of the interest and nuisance statistics) to a **log
  coefficient** (the natural log of the accumulated base-measure weight,
  kept in log scale for numerical stability).

When processing observation i: for every current bucket m and every
feasible increment yy \in \\0,\dots,\min(S_i,\\ M_s^{\mathrm{obs}}-m)\\,
the recursion moves probability mass from bucket m to bucket m+yy,
updates the key by t \mathrel{+}= yy\\ z_i and \mathbf{c} \mathrel{+}=
yy\\ \mathbf{x}\_{N,i}, and adds the log-weight increment yy\\ o_i -
\log(yy!) (the log of E_i^{yy}/yy!, the i-th factor of w_s). States
landing on the same key within a bucket are merged by **log-sum-exp**
(Section [9](#sec9)). After the last observation of the stratum, only
the bucket m = M_s^{\mathrm{obs}} is kept: this is exactly what enforces
the stratum-total constraint that defines \Omega_s. In the code this
recursion is `.joint_dist_stratum_dp_vec()` (`ep-multivariate-shift.R`),
implemented directly as a vectorized base-R operation: a
row-replication-and-grouping step over plain named lists of vectors,
rather than a nested per-value loop (Section [5.8](#sec5-8) gives the
engineering detail).

### 5.4 Combining strata (convolution)

Statistics *add* and weights *multiply* across strata (they are
independent given the stratum totals), so the global joint table is the
multi-dimensional convolution of the per-stratum tables: C(t,\mathbf{c})
= \sum\_{\substack{(t_1,\dots,t_Q):\\ \sum_s t_s = t \\
(\mathbf{c}\_1,\dots,\mathbf{c}\_Q):\\ \sum_s \mathbf{c}\_s =
\mathbf{c}}} \\ \prod\_{s=1}^Q C_s(t_s,\mathbf{c}\_s \mid
M_s^{\mathrm{obs}}). \qquad\text{(5.1)} This is computed by repeated
**pairwise** convolution. Start from the first stratum’s table (Section
[5.7](#sec5-7) fixes the processing order), and combine it, one stratum
at a time, with each remaining stratum’s table. For every pair of keys
\big((t_1,\mathbf{c}\_1),(t_2,\mathbf{c}\_2)\big), one from each table,
form the new key (t_1+t_2,\\ \mathbf{c}\_1+\mathbf{c}\_2) and the new
log-weight (the sum of the two), then merge duplicate resulting keys by
log-sum-exp. In the code, `.convolve_joint_pnuis_vec()` performs this
pairwise convolution in row-chunks (size set by `chunk_size`), to bound
peak memory when both tables being combined are large.

### 5.5 Pruning for feasibility (Rule A and Rule B)

The joint table grows with every nuisance axis and can get large. Two
safe pruning rules discard states that provably cannot reach the
observed nuisance total, and never discard a state that could still be
needed. Both rules are valid whenever every nuisance increment is
non-negative, which is the usual case for 0/1 dummy columns and
non-negative counts:

- **Rule A (upper-bound pruning).** For each nuisance component k: if
  the partial value c_k already *exceeds* the observed total
  t\_{N,\mathrm{obs},k}, discard the state. Non-negative increments mean
  nuisance totals never decrease as more strata or observations are
  added, so overshooting can never be corrected later.
- **Rule B (reachability pruning).** For stratum s and nuisance
  component k, define the per-stratum upper bound on how much that
  stratum alone could still add (Section [2.5](#sec2-5)’s U\_{s,k}),
  U\_{s,k} := \min\\\Big(M_s^{\mathrm{obs}},\\ \sum\_{i \in I_s:\\
  x\_{ik}\ne 0} S_i\Big), and, once some initial strata have already
  been combined, the total remaining upper bound U^{\mathrm{rem}}\_k :=
  \sum\_{\text{strata not yet processed}} U\_{r,k}. A partial state with
  nuisance value c_k is discarded if c_k + U^{\mathrm{rem}}\_k \<
  t\_{N,\mathrm{obs},k}: even the best possible outcome from every
  remaining stratum could not reach the observed total. (Section
  [2.5](#sec2-5) explains why M_s^{\mathrm{obs}} appearing inside this
  formula does not make it k-indexed in its own right.)

Both rules are exact: they never remove a state that belongs in the
final answer. They are applied directly inside the per-stratum recursion
(Rule A, via the `tN_obs` cap) and inside the combination step (Rule B,
via `remaining_max` in `.convolve_joint_pnuis_vec()`). Every stratum’s
own table C_s is built regardless of processing order (step 4 below). So
once all Q tables exist, their exact row counts are known at no extra
cost. Strata are then always combined smallest-table-first (ordering by
[`nrow()`](https://rdrr.io/r/base/nrow.html) of each stratum’s own
table, automatically, not a user-configurable option), because this
keeps the running accumulator as small as possible for as long as
possible, and makes Rule B’s bound tightest early, when it matters most.
A table’s row count, not the stratum’s raw observed total
M_s^{\mathrm{obs}}, is used for this ordering because it is the quantity
that actually determines combination cost, and it is available for free
here. Contrast this with
[`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md)
(Section [7.6](#sec7-6)), which never builds a per-stratum table at all
and must therefore fall back to the cheaper (but less exact) proxy
M_s^{\mathrm{obs}}, decided before any per-stratum work begins.

### 5.6 Filtering and marginalizing to C(t)

After the full combination over all Q strata, one joint table
C(t,\mathbf{c}) remains. Two final steps produce the one-dimensional
table needed for inference. First, **filter** to the observed nuisance
value (keep only rows with \mathbf{c} = \mathbf{t}\_{N,\mathrm{obs}}).
Then **marginalize**: drop the now-constant nuisance key, summing via
log-sum-exp any duplicate t values left over: C(t) = C(t,
\mathbf{t}\_{N,\mathrm{obs}}). This is `filter_and_marginalize_pn_lw()`
in the code. Its output, the pairs \\(t, \log C(t))\\, is exactly the
one-dimensional table Section [8](#sec8) uses for every exact-inference
quantity.

### 5.7 Algorithmic summary: the exact steps performed

For each requested effect (a single model-matrix column \mathbf{z}) in
`exact_vars`,
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)
performs:

1.  **Resolve the caps S_i** for every observation, using whichever of
    the three rules in Section [5.2](#sec5-2) applies (`EXACTMAX=var`,
    `EXACTMAX` without a variable / no offset, or `OFFSET` alone).
2.  **Build the model matrix**, identify the column \mathbf{z} for the
    current effect, and identify the nuisance columns X_N: either every
    other model-matrix column (`condition_on = "model"`, the
    SAS-matching default) or only the other requested exact effects
    (`condition_on = "exact"`).
3.  **Compute the observed nuisance total** \mathbf{t}\_{N,\mathrm{obs}}
    = X_N^\top y^{\mathrm{obs}} and the observed interest statistic
    t\_{\mathrm{obs}} = z^\top y^{\mathrm{obs}}.
4.  **For each stratum**, run the within-stratum dynamic program
    (Section [5.3](#sec5-3)) to get that stratum’s own coefficient table
    C_s(t_s,\mathbf{c}\_s \mid M_s^{\mathrm{obs}}), applying Rule A
    pruning against \mathbf{t}\_{N,\mathrm{obs}} as states are
    generated.
5.  **Order the strata** by ascending row count of their own just-built
    table (Section [5.5](#sec5-5)), and compute each stratum’s maximum
    possible contribution to every nuisance component, for use by Rule
    B.
6.  **Combine the per-stratum tables pairwise**, in chunks, applying
    Rule B pruning (using the running “remaining maximum” bound) after
    each stratum is folded in, to get the joint table C(t,\mathbf{c}).
7.  **Filter** the joint table to \mathbf{c} =
    \mathbf{t}\_{N,\mathrm{obs}} and **marginalize** to the
    one-dimensional table \\(t,\log C(t))\\\_{t\in\mathcal{T}}.
8.  **Compute exact inference** from this one table: the conditional
    MLE, the exact confidence interval, and the exact p-value (Section
    [8](#sec8)), and record the parameter’s row of `exactparmest` and
    its distribution in `OUTDIST`.
9.  **Repeat** steps 2-8 for every remaining requested effect.

------------------------------------------------------------------------

## 6. The closed-form algorithm (`ep_closed_form()`)

### 6.1 The special case

[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
targets one specific, very common special case of the general model:
exactly **one** two-level (treatment vs. control) effect of interest,
with **strata as the only other structure in the model** (no additional
categorical covariates to adjust for: p = 0, no X_N at all). This is the
classical setting of Guess, Lydick, Small and Miller (1987) and Guess
and Thomas (1990): a stratified comparison of two incidence rates from
person-time data.

For stratum s = 1,\dots,Q, let y_s^{(1)}, y_s^{(0)} be the observed
event counts in the treatment and control groups (aggregated within the
stratum), and e_s^{(1)}, e_s^{(0)} the corresponding exposures
(person-time). The model is Y_s^{(1)} \sim \mathrm{Poisson}(e_s^{(1)}
\lambda_s \theta), \qquad Y_s^{(0)} \sim \mathrm{Poisson}(e_s^{(0)}
\lambda_s), \qquad s=1,\dots,Q, equivalently \log\mathbb{E}\[Y_s^{(1)}\]
= \log e_s^{(1)} + \log\lambda_s + \beta, \log\mathbb{E}\[Y_s^{(0)}\] =
\log e_s^{(0)} + \log\lambda_s, \beta := \log\theta, where \lambda_s \>
0 is a stratum-specific nuisance baseline rate and \theta \> 0 is the
rate ratio. \theta carries no stratum index s, exactly as Guess and
Thomas (1990) state it in their own model formulation (“the relative
risk \omega = r\_{1i}/r\_{0i} is assumed to be constant across strata,”
their notation for this section’s \theta and s). This is exactly the
general model of Section [2.4](#sec2-4), with the offset o_i = \log
e_s^{(1)} or \log e_s^{(0)} depending on group, the stratum intercept
\tau\_{s(i)} = \log\lambda_s, the treatment indicator as z_i, and the
nuisance design matrix X_N empty. So the “constant across strata”
assumption here is not a simplification specific to this closed-form
special case. It is the same assumption Section [2.4](#sec2-4) already
makes for \beta in the fully general model, together with the
observation there that a fully free, stratum-specific \lambda_s
(baseline rate) does nothing to relax it: \lambda_s is free to differ
arbitrarily from stratum to stratum, while \theta, the ratio on top of
it, is assumed to be the same number everywhere. Like every other
function in this package, the exposure/offset is optional: if none is
supplied, every observation is treated as contributing one unit of
exposure (o_i = 0).

### 6.2 Deriving the within-stratum conditional distribution (Poisson thinning)

The nuisance parameter here is \lambda_s, and its sufficient statistic
is the stratum total T_s := Y_s^{(1)} + Y_s^{(0)}. Y_s^{(1)} and
Y_s^{(0)} are independent Poisson variables, so their sum is Poisson
too, with the rates added: T_s = Y_s^{(1)} + Y_s^{(0)} \sim
\mathrm{Poisson}\big(e_s^{(1)}\lambda_s\theta + e_s^{(0)}\lambda_s\big).
The conditional probability that Y_s^{(1)} = r (so Y_s^{(0)} = t_s-r),
given T_s = t_s, is the joint probability divided by this marginal:
\begin{aligned} \Pr(Y_s^{(1)}=r \mid T_s=t_s) &=
\frac{\Pr(Y_s^{(1)}=r)\Pr(Y_s^{(0)}=t_s-r)}{\Pr(T_s=t_s)}\\ &=
\frac{\dfrac{(e_s^{(1)}\lambda_s\theta)^r
e^{-e_s^{(1)}\lambda_s\theta}}{r!}\cdot\dfrac{(e_s^{(0)}\lambda_s)^{t_s-r}e^{-e_s^{(0)}\lambda_s}}{(t_s-r)!}}{\dfrac{(e_s^{(1)}\lambda_s\theta+e_s^{(0)}\lambda_s)^{t_s}e^{-(e_s^{(1)}\lambda_s\theta+e_s^{(0)}\lambda_s)}}{t_s!}}.
\end{aligned} The two exponentials in the numerator multiply to
e^{-(e_s^{(1)}\lambda_s\theta+e_s^{(0)}\lambda_s)}, which cancels
exactly against the denominator’s exponential. This is where \lambda_s
starts to disappear (a direct, concrete instance of the general
cancellation argument in Section [3.4](#sec3-4)). What remains is
\begin{aligned} \Pr(Y_s^{(1)}=r\mid T_s=t_s) &=
\binom{t_s}{r}\frac{(e_s^{(1)}\lambda_s\theta)^r(e_s^{(0)}\lambda_s)^{t_s-r}}{(e_s^{(1)}\lambda_s\theta+e_s^{(0)}\lambda_s)^{t_s}}\\
&=
\binom{t_s}{r}\left(\frac{e_s^{(1)}\theta}{e_s^{(1)}\theta+e_s^{(0)}}\right)^{\\r}\left(\frac{e_s^{(0)}}{e_s^{(1)}\theta+e_s^{(0)}}\right)^{\\t_s-r},
\end{aligned} \qquad\text{(6.1)} Every remaining factor of \lambda_s has
also cancelled, leaving an expression that depends only on the data and
\theta. This is exactly \mathrm{Binomial}(t_s,\\ p_s(\theta)), with
p_s(\theta) := \frac{e_s^{(1)}\theta}{e_s^{(1)}\theta+e_s^{(0)}}.
Conditioning on the stratum total does two things at once. It removes
\lambda_s completely (guaranteed in general by Section [3.2](#sec3-2),
since T_s is \lambda_s’s sufficient statistic), and it leaves behind a
Binomial count, which is simple enough to have a closed-form generating
function. The general theorem (Section [3.2](#sec3-2)) guarantees the
conditional distribution is free of \lambda_s, but it does not, by
itself, guarantee that distribution will be simple; that extra
convenience is Poisson-specific. Independence across strata follows
because \big(Y_1^{(1)},\dots,Y_Q^{(1)}\big) given (T_1,\dots,T_Q)
factors into the product of these Q conditional distributions.

### 6.3 The coefficient table as a polynomial in \theta

Writing (6.1) as a function of \theta, \Pr(Y_s^{(1)}=r\mid T_s=t_s) =
\left(\frac{e_s^{(0)}}{e_s^{(1)}\theta+e_s^{(0)}}\right)^{\\t_s}\binom{t_s}{r}\left(\frac{e_s^{(1)}}{e_s^{(0)}}\right)^{\\r}\theta^r,
the prefactor does not depend on r. All of the dependence on which value
Y_s^{(1)} takes is carried by
\binom{t_s}{r}(e_s^{(1)}/e_s^{(0)})^r\theta^r: the r-th term of the
binomial expansion of
\big(1+\tfrac{e_s^{(1)}}{e_s^{(0)}}\theta\big)^{t_s}. Conditional on
every T_s, the whole-sample sufficient statistic K := \sum\_{s=1}^Q
Y_s^{(1)} is a sum of independent (though not identically distributed)
Binomial counts. Its distribution as a function of \theta therefore
comes from **multiplying the Q per-stratum polynomials** and reading off
the coefficient of \theta^k: c_k = \sum\_{\substack{r_1+\cdots+r_Q=k \\
0\le r_s\le t_s}}\\ \prod\_{s=1}^Q
\binom{t_s}{r_s}\left(\frac{e_s^{(1)}}{e_s^{(0)}}\right)^{r_s}, \qquad
\Pr\nolimits\_\theta(K=k) =
\left\[\prod\_{s=1}^Q\left(\frac{e_s^{(0)}}{e_s^{(1)}\theta+e_s^{(0)}}\right)^{t_s}\right\]
c_k\\\theta^k. \qquad\text{(6.2)} Normalizing (6.2) over k removes the
bracketed factor, since it does not depend on k. This leaves the same
exponential-tilt form as the general algorithm (equation 3.3):
\Pr\nolimits\_\theta(K=k) = \frac{c_k\\\theta^k}{\sum_u c_u\\\theta^u},
\qquad \beta = \log\theta. \qquad\text{(6.3)} Combining strata is
therefore literally multiplying Q polynomials of degree t_1,\dots,t_Q: a
one-dimensional convolution of vectors of length t_s+1, computed stratum
by stratum by an accumulate-and-merge recursion (`compute_ck_natlog()`
in the code, working with natural-log coefficients \log c_k throughout,
purely for numerical stability, with no change to the mathematics
above). Section [6.3.1](#sec6-3-1) gives that recursion in full,
together with the specific bound that makes it fast.

### 6.3.1 The Guess-Thomas (1990) recursive algorithm, and why it is fast

**A local notation note before starting.** To match Guess and Thomas’s
(1990) own presentation, this subsection indexes strata by i=1,\dots,Q
(so “stratum i,” “the first i strata”) rather than this document’s usual
stratum letter s. This i is scoped entirely to this subsection and never
refers to an individual observation, unlike every other use of i
elsewhere in this document (Section [2.1](#sec2-1)). Similarly, L_i and
H_i below are local summation bounds specific to this recursion,
unrelated to the per-stratum bounds L\_{s,k} and U\_{s,k} of Section
[2.5](#sec2-5)/[5.5](#sec5-5); both pairs happen to use the letters L
and U/H for a similar purpose (bracketing a sum), but they are otherwise
different quantities defined for different purposes, and neither is used
again outside its own subsection.

Equation (6.2) could in principle be evaluated by expanding the
polynomial product directly, one stratum at a time, tracking every
coefficient up to the full degree \sum_s t_s at every intermediate
stage. That direct approach is exactly the original algorithm of Guess,
Lydick, Small, and Miller (1987), and its running time grows like T^Q (T
the largest stratum total, Q the number of strata). This is fine for a
handful of small strata, but prohibitive once either grows. Guess and
Thomas (1990) instead adapt the network-based recursive-convolution
technique that Mehta and Patel had already developed for related
exact-inference problems (Mehta and Patel 1983; Mehta, Patel, and Gray
1985) to this specific binomial-thinning setting. This gives an
algorithm with the same theoretical definitions and the same numerical
answers as the 1987 method, but a running time that grows only like
(ST)^2. On their own worked example (a table with six non-null strata
and 37 cases), the direct method took several hours on the hardware of
the time, against about 50 seconds for the recursive method: better than
a hundredfold speed-up on that example, and the method this package’s
own `compute_ck_natlog()` implements.

Define M(0,0) := 1 and, for stratum i=1,\dots,Q, the *partial*
coefficient table built from the first i strata alone, \begin{aligned}
M(i,k) &:=
\big\[\theta^k\big\]\prod\_{s=1}^i\left(1+\frac{e_s^{(1)}}{e_s^{(0)}}\theta\right)^{t_s}
\\ &= \sum\_{\substack{j_1+\cdots+j_i=k\\ 0\le j_s\le t_s}}\\
\prod\_{s=1}^i
\binom{t_s}{j_s}\left(\frac{e_s^{(1)}}{e_s^{(0)}}\right)^{j_s}
\qquad\text{so that}\qquad c_k = M(Q,k). \end{aligned}, M(i,k) satisfies
the one-stratum-at-a-time recursion obtained by multiplying in stratum
i’s own polynomial, \begin{aligned} M(i,k) = \sum\_{r=L_i}^{H_i}
M(i-1,\\ k-r)\\
\binom{t_i}{r}\left(\frac{e_i^{(1)}}{e_i^{(0)}}\right)^{r}, \\ \qquad
L_i := \max\big(0,\\ k-V\_{i-1}\big), \quad H_i := \min\big(k,\\
t_i\big), \end{aligned} where V_i := \min\\\big(N,\\ \sum\_{s=1}^i
t_s\big) and V_0 := 0. Here N denotes how far the table actually needs
to be built for the purpose at hand: it is not a fixed constant, but
whichever of the following two values the caller actually needs: the
observed total k\_{\mathrm{obs}} = \sum_s y_s^{(1)} for the truncated
table used by the p-value and confidence interval (Section
[6.5](#sec6-5)), or the full \sum_s t_s for the complete table the
conditional MLE needs (Section [8.2](#sec8-2)). This recursion computes
exactly the same c_k as the direct expansion of (6.2). It is the
identical polynomial multiplication, one stratum at a time, with one
addition: at every intermediate stage i, only coefficients k \le V_i are
ever computed or stored. This is a safe truncation. It never discards a
coefficient the final table could still need, because every stratum’s
contribution to k is non-negative: a partial total already above N can
never be brought back down to N or below by any later stratum, so it can
never contribute to a final coefficient c_k with k\le N. Bounding the
table’s width by N at every stage, rather than letting it grow with the
full (and possibly much larger) \sum\_{s\le i} t_s, is the entire source
of the (ST)^2 speed-up over the 1987 method’s T^Q. The two algorithms
compute the identical sequence of polynomial multiplications; they
differ only in how much of each intermediate result is worth keeping.

`compute_ck_natlog()` implements this recursion directly and by name.
`vvi <- pmin(max_cases, cumsum(t))` is exactly V_i for every stratum at
once (with `max_cases` the caller-supplied N, `max_k` in the function’s
own argument list, Section [6.5](#sec6-5)), and its nested loop over `k`
and then `r` is exactly M(i,k)’s own sum over r=L_i,\dots,H_i. Guess and
Thomas (1990) note that the terms M(i,k), c_k, and c_k\theta^k are prone
to the same overflow and underflow difficulties as the analogous terms
in Mehta and Patel’s own network algorithms, and they handle them the
same way: every one of these quantities is carried on the natural-log
scale, with any sum of two log-scale terms computed via the log-sum-exp
identity (Section [9](#sec9)) rather than by exponentiating and
re-summing directly. `compute_ck_natlog()` follows this exactly, hence
its name, and hence why every quantity in this section is written on the
log scale throughout.

### 6.4 Correspondence with the general coefficient table

c_k above is exactly the general coefficient table C(t,\mathbf{c}) of
Section [3.4](#sec3-4), specialized to p=0 (no nuisance columns at all),
with a per-observation cap S_i equal to the stratum total: precisely
SAS’s own Rule 2 (Section [5.2](#sec5-2)) for eliminating stratum
intercepts.
[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
has no `exactmax` argument, and its derivation above (Sections
[6.1](#sec6-1)-[6.2](#sec6-2)) conditions only on the stratum total T_s:
there is no place in it for a per-observation cap at all. So this
function *always* uses Rule 2, regardless of whether `offset` is
supplied to the call. This is unlike the two general algorithms
(Sections [5](#sec5) and [7](#sec7)), which switch automatically to Rule
3’s offset-implied cap whenever `offset` is supplied without `exactmax`
(Section [5.2](#sec5-2)). That difference is not just cosmetic: when
Rule 3’s cap would actually bind, the closed-form algorithm and the two
general algorithms are constraining \Omega differently and can
legitimately disagree (Section [10](#sec10)). Where the caps coincide
(no `offset`, or an `offset` too generous to bind), every algorithm in
this package eliminates stratum intercepts the *same* way: a per-stratum
sub-computation constrained to that stratum’s own total, followed by
combining strata. (Compare Sections [5.3](#sec5-3)-[5.4](#sec5-4) and
Sections [7.2](#sec7-2)-[7.6](#sec7-6) directly with Sections
[6.2](#sec6-2)-[6.3](#sec6-3): they are the same operation, specialized
or restructured.) With no nuisance columns to track, all of them
collapse to the univariate binomial/polynomial picture above. With one
or more nuisance columns, the per-stratum computation and the
cross-stratum combination have to carry those extra axes along too, at a
cost the closed-form algorithm never has to pay (Section [10](#sec10)).

### 6.5 Exact inference specifics: the truncated-table trick

Everything downstream follows the same pattern as the general algorithm
(Section [8](#sec8)), applied to the one-dimensional table \\(k,c_k)\\
of (6.2). One extra trick, specific to this algorithm’s
reparametrisation on \beta, is worth deriving explicitly, since it
accounts for much of its speed (Section [10](#sec10)).

\Pr\_\beta(K=k) in (6.3) always sums to one over k’s *complete* support,
**at every value of \beta**. That means the total normalizing constant
has a closed form that never requires summing the coefficient table at
all: Z(\beta) := \sum_u c_u\\ e^{\beta u} =
\left\[\prod\_{s=1}^Q\left(\frac{e_s^{(0)}}{e_s^{(1)}\theta+e_s^{(0)}}\right)^{t_s}\right\]^{-1}
= \prod\_{s=1}^Q
\left(\frac{e_s^{(1)}\theta+e_s^{(0)}}{e_s^{(0)}}\right)^{t_s}, \qquad
\theta=e^\beta. Equivalently, on the log scale, with
\log(1-p\_\theta(s)) := \log e_s^{(0)} - \log(e_s^{(1)}\theta+e_s^{(0)})
the per-stratum log-probability of the “control” outcome, \log Z(\beta)
= -\sum\_{s=1}^Q t_s\\\log\big(1-p\_\theta(s)\big) =: -g(\beta). So any
*tail* probability of K (\Pr\_\beta(K\le k_0) or \Pr\_\beta(K\ge k_0),
at **any** \beta, not only \beta=0) can be computed from a **truncated**
coefficient table (only the c_k for k on the near side of k_0), plus
this analytic normalizing factor, using the complement for the far side:
\Pr\nolimits\_\beta(K\le k_0) = \sum\_{k=0}^{k_0} c_k\\ e^{\beta k +
g(\beta)}, \qquad \Pr\nolimits\_\beta(K\ge k_0) = 1 -
\Pr\nolimits\_\beta(K\le k_0-1). The exact p-value and the exact
confidence interval only ever need \Pr\_\beta(K\le k\_{\mathrm{obs}})
and \Pr\_\beta(K\ge k\_{\mathrm{obs}}) (Section [8](#sec8)). So the
coefficient table only ever has to be built up to k = k\_{\mathrm{obs}}
for those two quantities (`compute_ck_natlog(..., max_k = sum(x))`, the
default). Only the conditional MLE genuinely needs the mean of K over
its *entire* range, and it is the only one that requires the complete
table (`max_k = sum(x) + sum(y)`). **This same fact (that a tail
probability at any \beta can be recovered from the already-built table
without re-enumerating anything) is also exactly what lets a p-value be
recomputed against a different null hypothesis, or a confidence interval
at a different \alpha, from a fit that has already run: see the
companion vignette,
[`vignette("Example")`](https://gsk-biostatistics.github.io/poissonexact/articles/Example.md),
for a worked example.**

### 6.6 Why the search over \beta is done the way it is

Both the exact CI limits and the exact p-value require finding a \beta
at which a tail probability (evaluated with the
truncated-table-plus-analytic-normalizer formula in Section
[6.5](#sec6-5)) equals a target value (for example, \alpha/2). This is a
one-dimensional root-finding problem in \beta. The implementation
(`solve_tail_root()`, shared with every other algorithm, Section
[8](#sec8)) finds a bracket containing the root by doubling a symmetric
interval around \beta=0 until the target function changes sign at its
two ends, then refines the root with a standard bisection/Brent-style
search ([`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html)).
Searching in \beta (the natural, unbounded parameter) with a
sign-change-based bracket guarantees the search is mathematically valid
no matter how extreme the true root is. A one-directional multiplicative
bracket expansion searched directly on the rate-ratio scale \theta \in
(0,\infty) can, for at least one realistic dataset, silently fail to
bracket the true root and converge to a value many orders of magnitude
off, because expanding the bracket in only one direction does not
correctly detect that the root is actually smaller than the current
guess. The \beta-scale, sign-change-based search used throughout this
package does not have that failure mode.

### 6.7 Algorithmic summary: the exact steps performed

1.  **Aggregate to per-stratum arm counts**: for each stratum s, compute
    the four numbers (y_s^{(1)}, y_s^{(0)}, e_s^{(1)}, e_s^{(0)})
    (treatment/control event counts and exposures) by summing over
    whichever individual rows of `data` fall in that stratum and that
    level of `exact_vars`. `data` may be individual-level or already
    aggregated; the summation makes the two equivalent as inputs.
2.  **Build the truncated coefficient table** \\(k, \log
    c_k)\\\_{k=0}^{k\_{\mathrm{obs}}} via `compute_ck_natlog()`, where
    k\_{\mathrm{obs}} := \sum_s y_s^{(1)} is the observed value of the
    whole-sample sufficient statistic.
3.  **Evaluate the two one-sided tail probabilities at \beta=0**
    (\Pr_0(K \le k\_{\mathrm{obs}}) and \Pr_0(K \ge k\_{\mathrm{obs}}))
    using the truncated table plus the analytic normalizing factor
    g(\beta) of Section [6.5](#sec6-5), and combine them into the
    two-sided exact p-value (Section [8.4](#sec8-4), shared with every
    other algorithm).
4.  **Find the exact confidence limits** by root-finding in \beta
    (Section [6.6](#sec6-6)): for a two-sided interval, solve
    \Pr\_\beta(K \ge k\_{\mathrm{obs}}) = \alpha/2 for the lower limit
    and \Pr\_\beta(K \le k\_{\mathrm{obs}}) = \alpha/2 for the upper
    limit, each using the same truncated-table trick; a limit is set to
    \mp\infty directly, without root-finding, if k\_{\mathrm{obs}} sits
    at the boundary of its attainable range.
5.  **Build the complete coefficient table** \\(k,\log
    c_k)\\\_{k=0}^{k\_{\mathrm{obs}}+\sum_s y_s^{(0)}}
    (`compute_ck_natlog(..., max_k = sum(x)+sum(y))`) and use it to
    solve for the conditional MLE \hat\beta (Section [8.2](#sec8-2),
    shared with every other algorithm), the one quantity that genuinely
    needs the whole table rather than a truncation of it.
6.  **Assemble the result**: one row of `exactparmest` (`Parameter`,
    `Estimate` = \hat\beta, `StdError`, `Alpha`, `LowerCL`, `UpperCL`,
    `PValue`, `SuffStat` = k\_{\mathrm{obs}}) and one entry of `OUTDIST`
    (the complete distribution \\(k,\Pr_0(K=k))\\, for inspection or
    plotting, and, as Section [6.5](#sec6-5) notes, for direct reuse
    afterward).

------------------------------------------------------------------------

## 7. The network algorithm (`ep_network()`)

### 7.1 Where this algorithm sits, and why it exists

Section [3.4](#sec3-4) establishes that the entire computational problem
of exact conditional Poisson regression is the problem of computing one
coefficient table, C(t) (or, before the final filtering step, the joint
table C(t,\mathbf{c}) over the interest statistic t and the nuisance
statistic \mathbf{c}). Section [5](#sec5)’s algorithm computes that
table with the multivariate shift algorithm of Hirji, Mehta, and Patel
(1987): build each stratum’s own table by a per-observation dynamic
program, then combine strata by explicit pairwise convolution of whole
tables. SAS calls this `EXACTOPTIONS METHOD=DIRECT`, and documents it as
computing the exact distribution directly, at the cost of potentially
large intermediate memory.

[`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md)
computes the *same* table (same definition, same final numbers, same
p-values and confidence intervals) by a different route: the network
algorithm of Mehta, Patel, and Senchaudhuri (1992), which SAS invokes by
default as `EXACTOPTIONS METHOD=NETWORK`. SAS’s own documentation for
that option is precise about the relationship between the two:
`METHOD=NETWORK` “builds a network for each parameter that you are
conditioning out, combines the networks, then uses the multivariate
shift algorithm to create the exact distribution.” In other words, the
network algorithm does not replace the multivariate shift algorithm. It
wraps it in a different data structure: a network of nodes and arcs,
processed with a backward pass before the forward one, so that the same
underlying recursion can discard impossible partial results the instant
they arise, instead of only at the boundary between one stratum’s table
and the next.

This section works through exactly what that network is, how the
backward pass on it is computed, and how the resulting bounds are used
to prune the forward pass continuously. Section [10](#sec10) then lays
out, point by point, how this differs from what Section [5](#sec5)’s
algorithm does with the same mathematics.

### 7.2 The network representation of a stratum’s reference set

Fix a stratum s and a requested effect with interest column \mathbf{z}
and nuisance columns X_N (Section [2.3](#sec2-3)). Write the stratum’s
observations, in some fixed order, as 1,\dots,n_s (so n_s := \|I_s\|,
the number of observations in stratum s alone: read as “n, restricted to
stratum s,” exactly as Q, the total number of strata, and
M_s^{\mathrm{obs}}, the total restricted to stratum s, relate to each
other elsewhere in this document; n = \sum_s n_s recovers the global
sample size of Section [2.1](#sec2-1)), with per-observation caps
S_1,\dots,S\_{n_s} (Section [5.2](#sec5-2)) and required stratum total
M_s^{\mathrm{obs}}.

**Nodes.** For i = 0,1,\dots,n_s, a node (i,m) represents the statement
“the first i observations of this stratum have contributed a combined
total of m toward M_s^{\mathrm{obs}}.” The start node is (0,0). The only
node at stage n_s that matters is the terminal node (n_s,
M_s^{\mathrm{obs}}), since every other value of m at the last stage
fails to satisfy the stratum-total conditioning constraint and is
infeasible.

**Arcs.** An arc leaves node (i,m) for every candidate value y\_{i+1}
\in \\0,1,\dots,S\_{i+1}\\ that observation i+1 could take, and arrives
at node (i+1,\\ m+y\_{i+1}). Each arc carries three lengths, accumulated
additively along any path: a **rank length** y\_{i+1}\\ z\_{i+1}
(contribution to the interest statistic t), one **nuisance length**
y\_{i+1}\\ x\_{N,i+1,k} per nuisance coordinate k=1,\dots,p, and a
**log-weight length** y\_{i+1}\\ o\_{i+1} - \log(y\_{i+1}!) (the log of
that observation’s own base-measure factor, Section [5.1](#sec5-1)).

**Paths and allocations.** A path from (0,0) to (n_s,
M_s^{\mathrm{obs}}) is a sequence of arcs, one per observation, and it
corresponds exactly, one-to-one, to a feasible allocation y_s =
(y_1,\dots,y\_{n_s}) \in \Omega_s (Section [5.3](#sec5-3) defines
\Omega_s). Summing a path’s rank lengths gives T_s(y_s). Summing its
k-th nuisance lengths gives \mathbf{T}\_{N,s}(y_s)\_k. Summing its
log-weight lengths gives \log w_s(y_s). The stratum’s coefficient table
C_s(t_s,\mathbf{c}\_s\mid M_s^{\mathrm{obs}}) is therefore exactly the
(log-sum-exp’d) collection of log-weight lengths of every path whose
rank length is t_s and whose nuisance lengths are \mathbf{c}\_s. The
network is nothing more than a bookkeeping device for listing every
allocation in \Omega_s without ever writing them out one at a time.

This is precisely the same object that Section [5.3](#sec5-3)’s dynamic
program already computes: its buckets, indexed by partial stratum total,
are this network’s nodes, and its per-observation update is this
network’s arcs. Building the network changes nothing about *what* is
computed at this stage. What it adds is Section [7.3](#sec7-3).

### 7.3 Backward induction: longest and shortest paths

For a node (i,m), define the **sub-network** \Gamma(i,m) as the set of
all paths from (i,m) to the terminal node. Equivalently, this is every
way the remaining observations i+1,\dots,n_s could be assigned values
that bring the running total from m up to exactly M_s^{\mathrm{obs}}.
For each nuisance coordinate k, define the **longest-path** and
**shortest-path** quantities \mathrm{LP}\_{c_k}(i,m) :=
\max\_{\text{paths in }\Gamma(i,m)}\\ (\text{sum of remaining
}k\text{-th nuisance lengths}), \mathrm{SP}\_{c_k}(i,m) :=
\min\_{\text{paths in }\Gamma(i,m)}\\ (\text{sum of remaining
}k\text{-th nuisance lengths}). These are exactly Mehta, Patel, and
Senchaudhuri’s (1992) longest-path (LP) and shortest-path (SP)
quantities (their eqs. 3.4-3.5), specialized here to the nuisance
coordinates, which are the only ones this algorithm’s pruning test
(Section [7.5](#sec7-5)) ever evaluates. They satisfy the backward
recursion (the source paper’s eq. 3.7, here for LP; SP is identical with
\min in place of \max): \begin{aligned} \mathrm{LP}\_{c_k}(i,m) &=
\max\_{yy \\\in\\ \\0,\dots,\min(S\_{i+1},\\ M_s^{\mathrm{obs}}-m)\\}
\Big\[\\ yy\\ x\_{N,i+1,k} \\+\\ \mathrm{LP}\_{c_k}(i+1,\\ m+yy)
\\\Big\],\\ \qquad \mathrm{LP}\_{c_k}(n_s, M_s^{\mathrm{obs}}) :&= 0,
\end{aligned} The \min and \max here play two different, non-competing
roles, worth separating explicitly: \max\_{yy\\\in\\\\\dots\\}\[\cdots\]
is the actual optimization, choosing the best candidate value yy; the
\min(S\_{i+1},\\ M_s^{\mathrm{obs}}-m) inside the braces merely sets the
*upper end of the set* yy is drawn from (observation i+1’s own cap,
further capped so the running total m+yy never exceeds the stratum’s
required total), and is fixed before the maximization ever starts.
Nothing is ever minimized over yy; the recursion for \mathrm{SP}\_{c_k}
two lines above swaps the outer \max for a \min but keeps this same
domain-defining \min(S\_{i+1}, M_s^{\mathrm{obs}}-m) unchanged. Computed
in one backward sweep, i = n_s-1, n_s-2,\dots,0, before any forward
processing begins. The node (i,m) is itself feasible (has *any*
completion at all) exactly when 0 \le M_s^{\mathrm{obs}}-m \le
\sum\_{j=i+1}^{n_s} S_j. Infeasible nodes carry no LP/SP value and are
never generated.

#### 7.3.1 Intuition, in plain terms: what “remaining nuisance length” actually means

The formulas above are worth pausing on, because everything else in this
section is built on the single idea they encode. This subsection gives
that idea without any formulas at all, then works through a full
numerical example, before returning to why it matters.

**Where the words come from.** Mehta, Patel, and Senchaudhuri were
originally computing the distribution of a *linear rank statistic*.
Wilcoxon, logrank, and similar tests are built from ranks, so they
called each arc’s contribution to the tracked quantity its “rank
length.” When that machinery is carried over to Poisson regression,
nothing is literally a rank anymore: the tracked quantity is now a
sufficient statistic, a nuisance column’s c_k; but the vocabulary
(“nuisance length,” the multi-dimensional cousin of “rank length”) is
kept because it is the source paper’s own terminology, and matching it
makes the two easy to cross-reference.

**Committed versus remaining.** Picture processing a stratum’s
observations one at a time, left to right, exactly as the forward pass
in Section [7.5](#sec7-5) does. At the moment you are standing at node
(i,m), you have already made a decision, explicitly or as part of
building up a partial record, for every observation 1,\dots,i. Those
decisions are locked in: they already contributed a fixed amount, m, to
the running total, and nothing that happens later can change that.
Observations i+1,\dots,n_s, by contrast, have not been decided yet.
Whatever they eventually turn out to be, they will add *some more* to
every nuisance coordinate before the stratum is done. “Remaining
nuisance length” is simply that not-yet-decided amount: the contribution
still to come from the part of the stratum that hasn’t been looked at
yet.

Because enumerating the remaining observations’ actual values is exactly
the whole point of the forward pass, there is no single number for “the
remaining amount.” There is only a *range* of possibilities: the most it
could still turn out to be (the longest path, LP: the best case for
reaching a target from below) and the least it could still turn out to
be (the shortest path, SP: the best case for staying at or under a
target from above). Backward induction is nothing more than computing
that range once, for every node, before the forward pass starts, so the
forward pass never has to guess at it later.

**A complete worked example.** Take a single nuisance coordinate k in a
stratum of n_s = 4 observations, with coefficients x\_{N,\cdot,k} =
(1,0,1,0) and caps S = (2,3,1,2), so the stratum’s required total is
M_s^{\mathrm{obs}} = 4. Suppose the forward pass has already committed
to y_1=1,\\ y_2=1 (one particular way of using up the first two
observations), landing on node (i,m) = (2,\\2). The remaining budget is
r = M_s^{\mathrm{obs}} - m = 4-2 = 2, to be split between observations 3
and 4. Observation 3 is “active” for coordinate k (x=1,
\mathrm{Cap}(A)=S_3=1); observation 4 is “inactive” (x=0,
\mathrm{Cap}(B)=S_4=2). Only two completions of (y_3,y_4) satisfy
y_3+y_4=2 within their caps:

| Completion: y_3 | Completion: y_4 | Check: y_3+y_4 | Result: contribution to coordinate k (= y_3, since obs. 4 doesn’t count) |
|:--:|:--:|:--:|:---|
| 0 | 2 | 2 | **0** (this is \mathrm{SP}\_{c_k}(2,2)) |
| 1 | 1 | 2 | **1** (this is \mathrm{LP}\_{c_k}(2,2)) |

Every feasible completion is listed here (there are only two, small
enough to enumerate directly), and they bracket the remaining
contribution to coordinate k between 0 and 1, exactly matching the
closed-form formula (Section [7.4](#sec7-4)): \mathrm{SP}\_{c_k}(2,2) =
\max(0,\\ r-\mathrm{Cap}(B)) = \max(0,\\ 2-2) = 0, and
\mathrm{LP}\_{c_k}(2,2) = \min(r,\\ \mathrm{Cap}(A)) = \min(2,1) = 1.
The closed form is not approximating anything here; it is just naming
the two extreme rows of a table that, in this small case, could be
written out by hand.

**Why this is worth computing at all.** Think of a nuisance coordinate
as a running score that has to land on one *exact* target by the end of
the stratum. That is what conditioning on \mathbf{t}\_{N,\mathrm{obs}}
means (Section [3.4](#sec3-4)). At any point mid-stratum, two questions
can be asked of the current partial record: “even in the best remaining
case, can this still reach the target?” and “even in the worst remaining
case, will this still avoid overshooting the target?” If either answer
is no, this partial record has no future. No matter what the
not-yet-processed observations turn out to be, it cannot end where it
needs to. So it can be thrown away right now, mid-stratum, at this exact
observation, rather than being carried forward, expanded by every
subsequent observation, and only discovered to be dead much later (in
Section [5](#sec5)’s table-then-convolve algorithm, only at a stratum
boundary or after an explicit cross-join with another stratum’s whole
table). Every one of those wasted expansions is real, avoidable work.
The backward-induction bounds are what make it avoidable.

### 7.4 The closed form used by this package, and why it is exact

Evaluating the recursion of Section [7.3](#sec7-3) by dynamic
programming, node by node, would cost about as much as the forward pass
itself. This package avoids that entirely, because of a restriction
Section [2.3](#sec2-3) already imposes for an unrelated reason
(“categorical effects only”): every column of \mathbf{z} and X_N is a
0/1 dummy. Under that restriction, the recursion has a closed-form
solution.

**Claim.** Fix a node (i,m), let r := M_s^{\mathrm{obs}} - m (the exact
total still to be allocated among observations i+1,\dots,n_s), and split
those remaining observations by their own value of x\_{N,\cdot,k} into
an **active set** A (coefficient 1) with total remaining capacity
\mathrm{Cap}(A) := \sum\_{j \in A} S_j, and an **inactive set** B
(coefficient 0) with capacity \mathrm{Cap}(B) := \sum\_{j\in B} S_j.
Then \mathrm{SP}\_{c_k}(i,m) = \max\big(0,\\ r - \mathrm{Cap}(B)\big),
\qquad \mathrm{LP}\_{c_k}(i,m) = \min\big(r,\\ \mathrm{Cap}(A)\big), and
every integer value in between is also achievable by some feasible
completion.

**Proof.** *Upper bound.* \sum\_{j\in A} y_j \le \mathrm{Cap}(A)
trivially (each y_j \le S_j), and \sum\_{j\in A} y_j \le r since
\sum\_{j\in B} y_j \ge 0 and the two sums must add to r. So \sum\_{j\in
A} y_j \le \min(r,\mathrm{Cap}(A)), and this is achieved by filling A’s
observations up to their own caps until either A’s share reaches r or
every A-observation is saturated (whichever comes first), with the
remainder placed on B (feasible since, when A saturates before reaching
r, the leftover r - \mathrm{Cap}(A) \le \mathrm{Cap}(B) follows from r
\le \mathrm{Cap}(A)+\mathrm{Cap}(B), the node’s own feasibility
condition).

*Lower bound.* Symmetrically, \sum\_{j\in A} y_j \ge r - \mathrm{Cap}(B)
(since \sum\_{j\in B} y_j \le \mathrm{Cap}(B)) and \ge 0 trivially,
achieved by filling B first and placing any unavoidable remainder on A.

*Every intermediate integer.* Starting from the lower-bound allocation
and, while A’s current share is below \min(r,\mathrm{Cap}(A)), moving
one unit of the total from some B-observation with unused capacity to
some A-observation with unused capacity (both exist as long as A’s share
has not yet reached its ceiling) increases A’s share by exactly one,
keeps the total fixed at r, and never violates any individual cap.
Repeating this move reaches every integer between the two bounds.
\blacksquare

Because the achievable set is a contiguous integer interval, both
endpoints are always exactly, not approximately, the true LP/SP values
for a 0/1-coefficient column. This closed form is a specialization of
the paper’s general backward-induction recursion to columns restricted
to \\0,1\\, not a relaxation or approximation of it. The whole-stratum
versions of these same two bounds, U\_{s,k} and L\_{s,k} (used to bound
what a stratum not yet reached could still contribute; see Section
[7.6](#sec7-6)), are the identical formula applied with r :=
M_s^{\mathrm{obs}} and A/B taken over the stratum’s *entire* set of
observations rather than only those remaining after position i. Section
[5.5](#sec5-5) already introduces U\_{s,k} for exactly this
whole-stratum quantity, and Section [7.6](#sec7-6) reuses it and its
counterpart L\_{s,k} without renaming them.

### 7.5 Forward processing with continuous pruning

Before stating the forward algorithm, one scope distinction has to be
made explicit, because it is the detail most likely to be misread
otherwise: **a live record’s key (t,\mathbf{c}) is cumulative across
every stratum processed so far. It is not reset at a stratum boundary.**
Only the within-stratum counter m resets to 0 once a stratum’s own
required total M_s^{\mathrm{obs}} has been enforced (Section
[7.6](#sec7-6) spells out the mechanics). t and every c_k keep
accumulating from the very first observation of the very first stratum
processed, straight through to the last observation of the last one,
because t and \mathbf{c} are, in the end, the *global*, whole-sample
statistics T(y) and \mathbf{T}\_N(y) of Section [3.1](#sec3-1), not a
per-stratum quantity. \mathrm{LP}\_{c_k}(i,m) and
\mathrm{SP}\_{c_k}(i,m) (Section [7.3](#sec7-3)), by contrast, are
defined purely *within* stratum s: \Gamma(i,m) is the set of paths from
(i,m) to stratum s’s own terminal node, and never refers to any other
stratum.

Certifying a record dead therefore needs *two* pieces combined, not one:
the within-stratum-s bound on what observations i+1,\dots,n_s could
still contribute (Sections [7.3](#sec7-3)-[7.4](#sec7-4)), plus a bound
on what every stratum *not yet started* could still contribute on top of
that. Section [7.6](#sec7-6) derives this second piece and calls it
U^{\mathrm{fut}}\_k(s) (LP side) and L^{\mathrm{fut}}\_k(s) (SP side);
it can be read first if its derivation is wanted before seeing it used
below. Comparing the cumulative, whole-sample c_k against the
whole-sample target t\_{N,\mathrm{obs},k} using only the within-stratum
piece, and leaving out U^{\mathrm{fut}}\_k(s)/L^{\mathrm{fut}}\_k(s),
would silently compare two quantities of different scope and produce an
incorrect (indeed, far too aggressive) test. The correct test always
includes both pieces together, as shown in step 3 below.

With that in view, the network is processed forward exactly as Mehta,
Patel, and Senchaudhuri’s (1992) five-step algorithm describes (their
Section 3.1.1; full citation in Section [11](#sec11)), generalized from
their single tracked statistic to this package’s joint (t,\mathbf{c})
key, currently inside stratum s:

1.  Select a live record (t,\mathbf{c},\ell) at node (i,m) (\ell its
    log-weight; t and \mathbf{c} are the *cumulative*,
    across-every-stratum-so-far values, while m is this stratum’s own
    running total).
2.  Expand it along every arc out of (i,m) (one per candidate value yy
    \in \\0,\dots,S\_{i+1}\\), transmitting a candidate record to node
    (i+1, m+yy) with updated key (t + yy\\z\_{i+1},\\ \mathbf{c} +
    yy\\\mathbf{x}\_{N,i+1},\\ \ell + yy\\o\_{i+1} - \log(yy!)).
3.  At the destination node, **prune immediately** if the candidate
    record’s nuisance key can be certified dead using that node’s own
    within-stratum LP/SP (Sections [7.3](#sec7-3)-[7.4](#sec7-4))
    *combined with* the future-strata bound
    U^{\mathrm{fut}}\_k(s)/L^{\mathrm{fut}}\_k(s) (Section
    [7.6](#sec7-6)):
    - already over the observed total (the ordinary Rule A: c_k \>
      t\_{N,\mathrm{obs},k}), or
    - cannot reach it even in the best case achievable anywhere in the
      remaining computation: the rest of this stratum *and* every
      stratum still to come (c_k + \mathrm{LP}\_{c_k}(i+1,m+yy) +
      U^{\mathrm{fut}}\_k(s) \< t\_{N,\mathrm{obs},k}), or
    - cannot avoid overshooting it even in the least case achievable
      anywhere in the remaining computation (c_k +
      \mathrm{SP}\_{c_k}(i+1,m+yy) + L^{\mathrm{fut}}\_k(s) \>
      t\_{N,\mathrm{obs},k}).
4.  Otherwise, merge the candidate record into any existing record at
    the destination node with the same key by log-sum-exp (Mehta, Patel,
    and Senchaudhuri’s “clubbing,” Section 3.1.1), or insert it as new.
5.  Repeat until every record at every node of stage i is processed,
    then advance to stage i+1.

The point of doing this stage by stage, node by node, rather than
building one complete per-stratum table first, is that a record failing
test 3 is discarded the moment it is created. It is never carried
forward to be expanded by further observations, and (Section
[10](#sec10)) never carried into an explicit cross-join against another
stratum’s table.

### 7.6 Combining strata

Mehta, Patel, and Senchaudhuri combine per-stratum networks by pairwise
convolution with an early-elimination rule (their Section 3.1.2):
sorting records and skipping any pair whose combined value, plus the
maximum every *not-yet-combined* stratum could still contribute, cannot
reach the target.

This package’s implementation realizes the same idea without a separate,
discrete convolution subroutine (`.network_stratum_slice()` in
`ep-network.R`). Strata are processed in sequence on **one continuous
frontier** of live records (this is exactly what makes t and \mathbf{c}
cumulative across strata, as Section [7.5](#sec7-5) emphasizes), ordered
by ascending observed total M_s^{\mathrm{obs}}. This is the cheapest
proxy available before any per-stratum work begins, since this
algorithm, unlike Section [5](#sec5)’s, never builds a per-stratum table
whose exact row count could be measured instead (Section
[5.5](#sec5-5)’s closing remark contrasts the two orderings directly).
Two whole-stratum quantities, computed once per stratum before any
forward processing begins, give the per-stratum contribution to Mehta,
Patel, and Senchaudhuri’s “every stratum not yet combined” bound: the
same U\_{s,k} already introduced in Section [5.5](#sec5-5) for the LP
side, and its direct counterpart for the SP side: \begin{aligned}
&\text{(LP side)}\quad U\_{s,k} = \min\\\Big(M_s^{\mathrm{obs}},\\
\textstyle\sum\_{i \in I_s:\\ x\_{N,i,k}=1} S_i\Big),\\ \qquad
&\text{(SP side)}\quad L\_{s,k} = \max\\\Big(0,\\ M_s^{\mathrm{obs}} -
\textstyle\sum\_{i\in I_s:\\ x\_{N,i,k}=0} S_i\Big), \end{aligned} the
same closed form as Section [7.4](#sec7-4), applied to an entire stratum
at once instead of to the observations remaining within one (Section
[7.4](#sec7-4)’s closing remark makes this identification explicit).
Summed over every stratum still to come in the processing order, these
give exactly the future-strata bound used in Section [7.5](#sec7-5)’s
pruning test: U^{\mathrm{fut}}\_k(s) := \sum\_{r\\ \text{processed
after}\\ s} U\_{r,k}, \qquad L^{\mathrm{fut}}\_k(s) := \sum\_{r\\
\text{processed after}\\ s} L\_{r,k} (`rem_future`/`min_future` in the
code, zero once every stratum has been processed). At the end of each
stratum’s own observations, live records are filtered to those whose
running total exactly equals that stratum’s M_s^{\mathrm{obs}}
(eliminating that stratum’s own intercept, Section [3.2](#sec3-2)), the
within-stratum counter m resets to zero, and
U^{\mathrm{fut}}\_k/L^{\mathrm{fut}}\_k advance to the next stratum’s
own values. t and every c_k carry over completely unchanged, exactly as
Section [7.5](#sec7-5) describes, since it is only ever the running
total m that is stratum-local.

### 7.7 Correctness: why the extra pruning never discards a valid state

Every pruning rule used in Section [7.5](#sec7-5) is a *necessary*
condition for a record to be extendable to a feasible allocation over
\Omega\_{\mathrm{cond}} (Section [3.3](#sec3-3)), not merely a
heuristic. Because c_k is the cumulative, whole-sample nuisance value
(Section [7.5](#sec7-5)), each rule must bound what remains *everywhere
still to come*: the rest of the current stratum, plus every stratum not
yet started:

- Rule A (c_k \le t\_{N,\mathrm{obs},k}) is necessary because every arc
  length is non-negative, so c_k can never decrease as more observations
  (in this stratum or any later one) are processed.
- The LP-based rule (c_k + \mathrm{LP}\_{c_k}(i,m) +
  U^{\mathrm{fut}}\_k(s) \ge t\_{N,\mathrm{obs},k}) is necessary because
  \mathrm{LP}\_{c_k}(i,m) + U^{\mathrm{fut}}\_k(s) is, by definition,
  the *largest* combined amount the rest of stratum s and every stratum
  still to come could possibly still contribute. If even that best case
  falls short, no completion of this record (within this stratum or
  across any later one) can reach the target.
- The SP-based rule (c_k + \mathrm{SP}\_{c_k}(i,m) +
  L^{\mathrm{fut}}\_k(s) \le t\_{N,\mathrm{obs},k}) is necessary for the
  mirror-image reason: \mathrm{SP}\_{c_k}(i,m) + L^{\mathrm{fut}}\_k(s)
  is the *smallest* possible combined remaining contribution from the
  rest of this stratum and every later one. If even that already
  overshoots, no completion can avoid exceeding the target.

Because each test only ever removes a record that has *zero* surviving
completions, no record that belongs in the final coefficient table
C(t,\mathbf{c}) is ever discarded. The algorithm computes exactly the
same table as Section [5](#sec5)’s algorithm, just without ever
materializing the states that Section [7.4](#sec7-4)’s closed form can
already certify as dead. This is the same safety argument Section
[5.5](#sec5-5) gives for Rules A and B, extended to the tighter,
per-node bounds derived here.

------------------------------------------------------------------------

## 8. The shared exact-inference machinery

Sections [5](#sec5)-[7](#sec7) describe three different ways of arriving
at the same object: the one-dimensional coefficient table \\(t,
C(t))\\\_{t\in\mathcal{T}} (Section [3.4](#sec3-4)), or, for the
closed-form algorithm, its polynomial-coefficient equivalent \\(k,
c_k)\\ (Section [6.3](#sec6-3)). These are the same mathematical object
under a different name and a different truncation strategy. Everything
from here on is a function of that one table alone, and is therefore
identical across all three algorithms. This is implemented once, in
`ep-utils.R`, and reused by every main function in the package.

### 8.1 Log-scale arithmetic and log-sum-exp

Every coefficient is tracked as \log C(t) (or \log c_k), not C(t)
itself, because the underlying weights w(y) can be astronomically large
or small (products of many terms E_i^{y_i}/y_i!), well outside the range
representable by ordinary floating-point numbers. Merging two states
that land on the same key requires adding their (non-log) weights:
computing \log(e^{a}+e^{b}) from a=\log(\cdot) and b=\log(\cdot) without
ever forming e^a or e^b directly. This is the standard log-sum-exp
trick, \log(e^a+e^b) = m + \log(e^{a-m}+e^{b-m}) with m:=\max(a,b),
implemented once as `log_sum_exp()` (a single pair) and
`log_sum_exp_by()` (grouped merging of many keys at once, used
throughout Sections [5](#sec5) and [7](#sec7)).

### 8.2 The conditional maximum likelihood estimate (CMLE)

Given the *complete* table \\(t,\log C(t))\\ (every algorithm needs the
complete table for this one quantity, even the closed-form algorithm’s
truncated-table trick, Section [6.5](#sec6-5)), the conditional MLE
\hat\beta solves \mathbb{E}\_{\hat\beta}\\\[T\] = t\_{\mathrm{obs}},
\qquad \mathbb{E}\_\beta\\\[T\] := \sum_t t\\ p\_\beta(t) for \hat\beta,
using p\_\beta(t) from (3.3)/(6.3). This is a monotone equation in \beta
(the exponential family’s mean is increasing in its canonical
parameter), so it has a unique solution whenever t\_{\mathrm{obs}} is
strictly inside the attainable range of T. \hat\beta = \mp\infty by
definition if t\_{\mathrm{obs}} sits at the boundary.

The standard error reported alongside \hat\beta comes from the
conditional variance of T, \mathrm{Var}\_\beta\\\[T\] :=
\sum\_{t\in\mathcal{T}} \big(t - \mathbb{E}\_\beta\\\[T\]\big)^2\\
p\_\beta(t), evaluated at \beta=\hat\beta: since
\mathbb{E}\_\beta\\\[T\] is the derivative of the log-partition function
and \mathrm{Var}\_\beta\\\[T\] is its second derivative,
\mathrm{Var}\_{\hat\beta}\\\[T\] *is* the Fisher information for \beta
at \hat\beta, so \widehat{\mathrm{SE}}(\hat\beta) =
\big(\mathrm{Var}\_{\hat\beta}\\\[T\]\big)^{-1/2}. In the code,
`solve_cmle()` implements exactly this, shared verbatim by every
algorithm.

### 8.3 Root-finding scaffold

Both the CMLE (Section [8.2](#sec8-2)) and the confidence interval
(Section [8.4](#sec8-4)) reduce to the same shape of problem: find \beta
solving g(\beta) = \text{target} for some function g that is monotone in
\beta: g(\beta) = \mathbb{E}\_\beta\\\[T\] and target =t\_{\mathrm{obs}}
for the CMLE (Section [8.2](#sec8-2)); g(\beta) = \Pr\_\beta(T\ge
t\_{\mathrm{obs}}) or \Pr\_\beta(T\le t\_{\mathrm{obs}}) and target
=\alpha/2 (or \alpha, one-sided) for the confidence limits (Section
[8.4](#sec8-4) below).

`solve_tail_root()` is the one shared implementation of that root-find.
It first finds a bracket \[-b,b\] containing the root by doubling b
(starting from a small initial value) until
\big(g(-b)-\text{target}\big)\cdot\big(g(b)-\text{target}\big) \le 0,
i.e. until g(\beta)-\text{target} changes sign between the two ends,
then refines the root inside that bracket with a standard
bisection/Brent-style search
([`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html)). What
differs between algorithms is only *how* g(\beta) itself is evaluated at
a candidate \beta: from a complete table via `tilt_probs()` (Sections
[5](#sec5) and [7](#sec7)’s algorithms), or from a truncated table plus
an analytic normalizing factor (the closed-form algorithm, Section
[6.5](#sec6-5)). The search strategy itself never differs.

### 8.4 P-values and confidence intervals

Define the two one-sided tail probabilities of the conditional
distribution (3.3)/(6.3) at any \beta, P\_\beta^{\le} :=
\Pr\nolimits\_\beta(T \le t\_{\mathrm{obs}}) =
\sum\_{\substack{t\in\mathcal{T}\\ t\le t\_{\mathrm{obs}}}} p\_\beta(t),
\qquad P\_\beta^{\ge} := \Pr\nolimits\_\beta(T \ge t\_{\mathrm{obs}}) =
\sum\_{\substack{t\in\mathcal{T}\\ t\ge t\_{\mathrm{obs}}}} p\_\beta(t).
The two-sided exact p-value is, for every algorithm, twice the smaller
of these two tail probabilities at the null (\beta=0), capped at 1:
\text{p-value} = \min\\\Big(1,\\ 2\min\big(P_0^{\le},\\
P_0^{\ge}\big)\Big) \qquad (\texttt{two\\sided\\from\\tails()}). Which
side is smaller also fixes which one-sided confidence limit is
informative when `ci_type = "one-sided"`. The confidence interval itself
is built by inverting these same tail probabilities (a
Clopper-Pearson-style construction), and reduces to the same four cases
everywhere, implemented once as `exact_ci_from_tails()`:

- **Two-sided, level 1-\alpha.** Solve for the lower and upper limits
  separately, \beta_L:\\ \Pr\nolimits\_{\beta_L}(T \ge
  t\_{\mathrm{obs}}) = \alpha/2, \qquad \beta_U:\\
  \Pr\nolimits\_{\beta_U}(T \le t\_{\mathrm{obs}}) = \alpha/2, each via
  the root-finding scaffold of Section [8.3](#sec8-3).
- **One-sided, upper limit only.** Solve \Pr\_{\beta_U}(T \le
  t\_{\mathrm{obs}}) = \alpha (using the full \alpha, not \alpha/2),
  with the lower limit reported as -\infty.
- **One-sided, lower limit only.** Solve \Pr\_{\beta_L}(T \ge
  t\_{\mathrm{obs}}) = \alpha, with the upper limit reported as +\infty.
- **Boundary case, any of the above.** If t\_{\mathrm{obs}} sits at the
  edge of T’s attainable range, the corresponding limit is set to
  \mp\infty by definition rather than by root-finding (there is no
  finite \beta at which the relevant tail probability equals \alpha or
  \alpha/2).

`exact_ci_from_tails()` takes two caller-supplied tail-probability
callbacks rather than a table, so it can sit in front of either
representation described in Section [8.3](#sec8-3).

**Nothing in this section is specific to testing the null \beta=0, or to
any one value of \alpha.** \beta=0 and whatever `alpha` was passed to
the original call are simply the values plugged in when a fit first
runs. `tilt_probs()`, `two_sided_from_tails()`, and
`exact_ci_from_tails()`/`exact_ci_from_dist1d()` are all generic in
\beta_0 (the null being tested) and \alpha (the confidence level)
respectively; the *table* they operate on (`OUTDIST`’s
`SuffStat`/`LogWeight` columns) does not change when either of those is
changed. This is exactly the fact the companion vignette,
[`vignette("Example")`](https://gsk-biostatistics.github.io/poissonexact/articles/Example.md),
puts to use: once a fit has produced its `OUTDIST` table, a p-value
against a different null rate ratio, or a confidence interval at a
different \alpha, is cheap to recompute from that same table, without
repeating the combinatorial work of Sections [5](#sec5)-[7](#sec7) that
built it in the first place.

------------------------------------------------------------------------

## 9. Numerical implementation notes

- **Everything is computed in log space.** Every coefficient, every
  intermediate weight, and every probability calculation works with \log
  C(t) (or \log c_k) rather than C(t) itself, precisely because the
  underlying base-measure weights routinely underflow or overflow
  ordinary floating-point range (Section [8.1](#sec8-1)). Merges use
  `log_sum_exp()`/`log_sum_exp_by()` throughout.
- **Rule A/B pruning is applied as early as computationally
  convenient.** The multivariate-shift algorithm (Section [5](#sec5))
  applies Rule A inside the per-stratum recursion and Rule B once per
  stratum, at the combination step. The network algorithm (Section
  [7](#sec7)) applies both directions of the same idea (plus their
  SP-side counterpart) continuously, at every single observation, using
  the closed-form bound of Section [7.4](#sec7-4). Both are exact:
  neither ever discards a state that belongs in the final answer
  (Sections [5.5](#sec5-5) and [7.7](#sec7-7)). They differ only in how
  early the same provably-safe conclusion is acted on.
- **Strata are always combined smallest-first, automatically, never as a
  user-configurable option**, but “smallest” is measured differently in
  the two general algorithms, deliberately.
  [`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)
  builds every stratum’s own table regardless of processing order, so
  once built, each table’s exact row count is free to measure and is
  used directly (Section [5.5](#sec5-5)).
  [`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md)
  never builds a per-stratum table at all, so it must decide its
  processing order beforehand from the cheaper proxy M_s^{\mathrm{obs}},
  the raw observed stratum total (Section [7.6](#sec7-6)). Both
  orderings serve the same goal: keep the running accumulator or
  frontier small for as long as possible, and Rule B’s bound tight
  early. Each algorithm uses the best signal its own design makes
  available.
- **Combining two strata’s tables is done in row-chunks, controlled by
  `chunk_size`.** Combining strata (Section [5.4](#sec5-4)) requires
  pairing every row of one stratum’s table with every row of the
  other’s. `chunk_size` bounds how many such pairs are materialized
  simultaneously, trading memory for speed.
  [`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md)
  (Section [7](#sec7)) never forms this explicit cross-join at all: its
  continuous frontier is bounded by the same Rule A/B/SP-side pruning
  applied per observation instead. So it has no `chunk_size` argument.
- **The per-observation cap rule is resolved automatically, not chosen
  by the user.** Section [5.2](#sec5-2)’s three SAS-style rules are
  selected from whichever of `exactmax`/`offset` was actually supplied
  to the call. There is no argument for substituting a different rule
  for computing S_i.
- **Cross-algorithm agreement is verified directly, not assumed.** All
  three algorithms are independently implemented, yet compute the same
  table by construction (Section [3.4](#sec3-4)) whenever they are given
  the same per-observation caps S_i (Section [10](#sec10); this holds
  automatically unless `offset` is supplied to the general algorithms
  without `exactmax` and the resulting Rule 3 cap actually binds, in
  which case
  [`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
  is solving a different, cap-free problem by design). The package’s
  test suite checks agreement directly on cases in scope:
  [`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md),
  [`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md),
  and
  [`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md)
  are run on the same model, and their `exactparmest`/`OUTDIST` outputs
  compared to numerical tolerance (see
  `test-cross-algorithm-agreement.R`), rather than relying on the
  mathematics of Section [3.4](#sec3-4) alone; that same file also
  documents the expected exception with a worked example.

------------------------------------------------------------------------

## 10. Comparing the three algorithms, and choosing among them

All three algorithms compute the same object, C(t,\mathbf{c}) (Section
[3.4](#sec3-4)), given the same per-observation caps S_i, and therefore
give identical exact inference (the same p-value, the same confidence
interval, the same conditional MLE, up to ordinary floating-point
agreement) **whenever they are actually solving the same problem**. They
differ only in *how* they get there, and in what that costs.

That qualification is not a formality, and it is worth stating precisely
because it is easy to miss.
[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
has no `exactmax` argument and never applies Rule 1 or Rule 3 (Section
[5.2](#sec5-2)): its own derivation (Section [6.2](#sec6-2)) conditions
only on the stratum total, so it always uses Rule 2’s cap (S_i = stratum
total), regardless of whether `offset` is supplied to the call (Section
[6.4](#sec6-4)).
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)
and
[`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md),
by contrast, apply Rule 3 automatically whenever `offset` is supplied
without `exactmax`: S_i = \max(\lceil e^{o_i}\rceil,\\ 1). When that cap
does not bind (every S_i at least as large as the stratum total, so it
excludes nothing beyond what Rule 2 already implies), all three
algorithms are constraining \Omega identically and therefore agree
exactly. When it does bind (small exposures relative to the observed
counts, typically in unstratified or lightly-stratified two-group data),
the closed-form algorithm and the two general algorithms are, by design,
no longer answering the same question: the closed-form algorithm’s
reference set is every non-negative allocation consistent with the
stratum total alone, while the general algorithms’ reference set
additionally excludes any allocation where some y_i exceeds its own
offset-implied bound. The two reference sets differ, so the exact
p-value, confidence interval, and conditional MLE can legitimately
differ too. This is expected behaviour, a genuine difference in which
modelling assumption each algorithm encodes (Section [5.2](#sec5-2)),
not a bug in either one.
`tests/testthat/test-cross-algorithm-agreement.R` includes a small
worked example (a two-observation, unstratified fit where the
offset-implied caps bind asymmetrically) that reproduces exactly this
divergence and locks it in as expected, not regression, behaviour.

### 10.1 Why the closed-form algorithm is faster, precisely, in the case it covers

The closed-form algorithm (Section [6](#sec6)) is restricted to p=0 (no
nuisance columns beyond the stratum intercepts) and a single two-level
effect. Within that restriction, it is faster than the general
algorithms for two independent reasons. First, it never tracks a
nuisance key at all: there is no \mathbf{c} dimension to carry through
every stratum’s computation and every combination step, only the scalar
t, while the general algorithms must carry one axis per nuisance
coordinate through every stratum and every pairwise combination (Section
[4.6](#sec4-6)). Second, its exponential-tilt family has a closed-form
total (Section [6.5](#sec6-5)), so its p-value and confidence interval
only ever need a *truncated* table plus an analytic normalizing factor,
never the complete table that the general algorithms’ tail probabilities
are computed from directly. Both effects are specific to p=0 and the
binomial-thinning structure of Section [6.2](#sec6-2). Neither survives
once a genuine nuisance column is added, which is exactly why the
closed-form algorithm cannot be generalized to that case, and why the
algorithms of Sections [5](#sec5) and [7](#sec7) exist.

### 10.2 A worked example where the algorithms genuinely disagree

The case below is deliberately the smallest one that can trigger the
disagreement introduced above, namely a single, unstratified two-level
comparison with an offset and no `exactmax`:

``` r
binom5 <- data.frame(group = c(1, 0), cases = c(3, 2), fu = c(4, 3))
ep_multi_shift(
  cases ~ group, data = binom5, offset = "fu", offset_log = FALSE, exact_vars = "group"
)
ep_closed_form(
  cases ~ group, data = binom5, offset = "fu", offset_log = FALSE, exact_vars = "group"
)
```

There is one stratum and two observations: raw exposures E_1 = 4 (the
`group = 1` row’s `fu`) and E_0 = 3 (the `group = 0` row’s), with
observed counts x_1^{\mathrm{obs}} = 3 and x_0^{\mathrm{obs}} = 2, so
M^{\mathrm{obs}} = x_1^{\mathrm{obs}}+x_0^{\mathrm{obs}} = 5.

**Same base-measure weights, different support.** Both algorithms build
their reference distribution for the group-1 count r from the
*identical* per-value weight (Section [5.1](#sec5-1)), w(r) := E_1^r/r!
\cdot E_0^{5-r}/(5-r)! = 4^r/r! \cdot 3^{5-r}/(5-r)!:

| r    | 0     | 1    | 2   | 3   | 4   | 5      |
|------|-------|------|-----|-----|-----|--------|
| w(r) | 2.025 | 13.5 | 36  | 48  | 32  | 8.5333 |

[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
places no cap beyond the stratum total (Section [6.4](#sec6-4)), so its
reference set is every value r \in \\0,\dots,5\\: at \beta=0 these six
weights sum to 140.0583, and the observed value r=3 carries null
probability 48/140.0583 = 0.3427 (`fit_closed$OUTDIST` reports
0.34271435, up to rounding, and its `LogWeight` column is this same
table on the algebraically reduced Guess-Thomas scale,
\ln\\\big(C(5,r)(4/3)^r\big), a constant, r-independent shift in log
space that leaves every probability unchanged).

[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)
and
[`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md),
given `offset` without `exactmax`, apply Rule 3 to each observation
separately (Section [5.2](#sec5-2)): S_1 = \lceil 4 \rceil = 4, and S_0
= \lceil 3\rceil, which is 3 mathematically, but the value the code
actually computes is 4. `si_from_offset()` computes S_i as
`ceiling(exp(log_offset))`, and in IEEE double precision `exp(log(3))`
is not exactly `3`; it evaluates to `3.0000000000000004`, whose ceiling
is `4`, not `3`. So the true caps used here are S_1 = S_0 = 4.
Feasibility requires r \le S_1 = 4 and 5-r \le S_0 = 4 (equivalently r
\ge 1) at once, so the reference set is r \in \\1,2,3,4\\: the *same
four weights* 13.5, 36, 48, 32, now summing to 129.5, with the observed
value r=3 carrying null probability 48/129.5 = 0.3707
(`fit_shift$OUTDIST` reports `SuffStat` 1{:}4 with `prob` 0.1042471,
0.2779923, 0.3706564, 0.2471042, an exact match, and its `LogWeight`
column, 2.602690, 3.583519, 3.871201, 3.465736, is \ln w(r) for
r=1,\dots,4 directly, on the *un*-reduced base-measure scale of Section
[5.1](#sec5-1), which is why it differs from
[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)’s
own `LogWeight` column by the same constant mentioned above).

This is a genuine, if narrow, numerical subtlety worth being direct
about: whether Rule 3’s cap for an exactly-integer-valued offset lands
on that integer or one more than it can depend on sub-10^{-15} rounding
in the log/exp round trip, not on anything in the data. It is not a
logic error (`si_from_offset()` computes precisely the formula in
Section [5.2](#sec5-2), `ceiling(exp(log_offset))`, faithfully), and it
has no practical bearing on real exposure data (person-time is
essentially never exact to 16 significant digits in the first place). It
simply means a hand-derivation like the one two paragraphs above cannot
stop at “\lceil 3 \rceil = 3” and expect it to match what the code
actually did for an integer-valued raw offset; the actual S_i used is
the safer thing to check (e.g. against `OUTDIST`’s own support) rather
than assumed from the offset value alone.

Whichever exact caps apply, the mechanism is the same: the weights
themselves never change: only which of them are eligible to be summed.
That alone is enough to change everything downstream: the conditional
MLE solves a mean-matching equation (Section [8.2](#sec8-2)) against a
different reference distribution in each case, and the exact tail
probabilities that define the p-value and confidence limits (Section
[8.3](#sec8-3)) are integrated over a different support in each case.
Running the two calls above gives exactly this:

|  | [`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md) (= [`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md)) | [`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md) |
|----|----|----|
| Support of r | \\1,2,3,4\\ | \\0,1,2,3,4,5\\ |
| Estimate (\hat\beta) | 0.2828976 | 0.1177830 |
| StdError | 1.117743 | 0.9128709 |
| 95% CI | (-2.386108,\\ 4.081686) | (-2.048941,\\ 2.600418) |
| PValue | 1 | 1 |
| SuffStat | 3 | 3 |

The p-value happens to agree here (the observed value is unremarkable
under either reference distribution), but the point estimate and
confidence interval do not, by a wide margin, and this is not a bug in
either function. The two calls are answering genuinely different
questions about which values r was even permitted to take, because they
encode different modelling assumptions (Rule 2 vs. Rule 3, Section
[5.2](#sec5-2)), not because either implementation made a mistake.

If the offsets are widened just enough that the cap stops excluding
anything: `fu = c(5, 5)`, so S_1 = S_0 = 5, both algorithms build the
identical full support \\0,\dots,5\\ again, and every quantity above
matches to the last reported digit:

|  | [`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md) (= [`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md)) | [`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md) |
|----|----|----|
| Support of r | \\0,\dots,5\\ | \\0,\dots,5\\ |
| Estimate (\hat\beta) | 0.4054651 | 0.4054651 |
| StdError | 0.9128709 | 0.9128709 |
| 95% CI | (-1.761258,\\ 2.8881) | (-1.761258,\\ 2.8881) |
| PValue | 1 | 1 |
| SuffStat | 3 | 3 |

`tests/testthat/test-cross-algorithm-agreement.R` locks in both halves
of this example directly against `OUTDIST`’s support (`0:5` for
[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
throughout; `1:4` when the cap binds and `0:5` when it does not, for the
two general algorithms), so this is exercised by the test suite, not
only described here.

### 10.3 Multivariate-shift versus network, precisely

Where the closed-form algorithm does not apply (any nuisance columns
beyond stratum intercepts, or more than one two-level effect at a time),
Sections [5](#sec5) and [7](#sec7) compute the identical table by
different strategies:

|  | [`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md) (multivariate shift / SAS `DIRECT`) | [`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md) (network / SAS `NETWORK`) |
|----|----|----|
| Source | Hirji, Mehta, and Patel (1987) | Mehta, Patel, and Senchaudhuri (1992) |
| Per-stratum computation | Per-observation recursion builds one **complete** table C_s(t_s,\mathbf{c}\_s) | Same per-observation recursion, but every record is checked against node-level LP/SP the instant it is created |
| Where pruning happens | Rule A inside the recursion; Rule B only at stratum/combination boundaries (once per stratum) | Rule A **and** the LP/SP-based generalization of Rule B, continuously, at every observation |
| How strata are combined | Explicit pairwise cross-join of two whole tables, chunked (`chunk_size`), **then** masked by Rule A/B | One continuous frontier across all strata; between-strata bound folded into the same continuous test, no cross-join ever materialized |
| Bound direction used | Upper bound only (cannot exceed the observed total; cannot exceed what remains reachable) | Both directions: upper (LP, reachability) and lower (SP, cannot-avoid-overshoot) |
| Extra precomputation required | None beyond what the recursion itself produces | One backward, closed-form sweep per stratum (Section [7.4](#sec7-4)), before the forward pass |
| Output | `exactparmest`, `OUTDIST` | Identical in definition and (up to floating-point agreement) in value, checked directly by the test suite (`test-cross-algorithm-agreement.R`) |
| When it should be faster | Simpler to reason about; fine when tables stay small | Large or many strata, several nuisance coordinates, imbalanced data, exactly SAS’s own stated advantage for `METHOD=NETWORK` |

[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)’s
cost is dominated by the size of the explicit cross product formed at
each pairwise combination. `nrow(table_1) x nrow(table_2)` candidate
rows are always materialized before Rule A/B ever removes any of them.
[`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md)
never forms that product: a record is checked against the LP/SP bound
the instant it is created, at essentially no cost beyond the arithmetic
in Section [7.4](#sec7-4)’s closed form (two comparisons per nuisance
coordinate, using quantities already precomputed once per stratum). The
live-state count at any point in the computation is therefore governed
by genuine reachability rather than by the size of an intermediate cross
product. This is precisely the property that lets SAS’s `METHOD=NETWORK`
outperform `METHOD=DIRECT` on the same problem, and it is the property
[`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md)
is built to reproduce. The gain is largest exactly where SAS’s own
documentation says it is largest: heavily imbalanced strata and models
with several nuisance coordinates, where a large fraction of
naively-generated states are in fact unreachable and would otherwise
survive all the way to a cross-join before being wasted.

### 10.4 Practical guidance

- If the model has exactly one two-level exact effect and strata as the
  only other structure, use
  [`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md).
  It computes the identical answer at substantially lower cost (Section
  [10.1](#sec10-1)), and there is no reason to prefer either general
  algorithm in this case, *provided* it is actually being asked the same
  question as they are (see the next point).
- If you supply `offset` without `exactmax` and compare against
  [`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md),
  be aware of Section [10.2](#sec10-2):
  [`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
  never applies the offset-implied per-observation cap that the general
  algorithms apply automatically in that situation. The two are then
  deliberately solving different problems, and comparing them is only
  meaningful once you have decided which cap (Rule 2, the plain stratum
  total, or Rule 3, the offset-implied bound) matches the assumption you
  actually want to make about individual counts.
- Otherwise (any additional categorical covariates to adjust for, or
  more than one exact effect), use one of the general algorithms.
  [`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)
  is the simplest to reason about and a good default: its per-stratum
  recursion and cross-stratum combination are already implemented in
  vectorized form (Section [5.8](#sec5-8)), so there is no separate,
  slower variant to graduate away from if it turns out to be too slow.
  [`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md)
  computes the identical result again, via a structurally different,
  more aggressively pruned strategy (Section [7](#sec7)), and should be
  tried when the model has several nuisance coordinates or heavily
  imbalanced strata: precisely the regime where SAS’s own
  `METHOD=NETWORK` is documented to outperform `METHOD=DIRECT`.
- The package’s own test suite (`test-cross-algorithm-agreement.R`)
  checks that all three algorithms agree with each other, across a range
  of models, to numerical tolerance, and, alongside that, locks in
  Section [10.2](#sec10-2)’s example as an *expected exception* to that
  agreement. Running a given model through more than one algorithm and
  comparing `exactparmest` is cheap insurance before relying on a result
  for anything important:
  [`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)
  and
  [`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md)
  are independently implemented yet apply identical cap rules, so they
  should never disagree except through a genuine bug;
  [`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
  should additionally agree with both of them whenever its own no-cap
  assumption coincides with whatever cap the general algorithms actually
  used (Section [10.2](#sec10-2)), but is not expected to agree with
  them when it does not. The companion vignette,
  [`vignette("Example")`](https://gsk-biostatistics.github.io/poissonexact/articles/Example.md),
  does exactly this kind of cross-check on a worked example.

------------------------------------------------------------------------

## 11. References

- Agresti A. *Categorical Data Analysis*, 3rd ed. Wiley, 2013. ISBN
  978-0-470-46363-5.
- Andersen EB. Sufficiency and exponential families for discrete sample
  spaces. *Journal of the American Statistical Association*,
  65(331):1248-1255, 1970.
  <https://doi.org/10.1080/01621459.1970.10481160>
- Cox DR. *The Analysis of Binary Data*. Methuen, 1970. ISBN
  978-0-416-10400-4.
- Guess HA, Lydick EG, Small RD, Miller LP. Exact binomial confidence
  intervals for the relative risk in follow-up studies with sparsely
  stratified incidence density data. *American Journal of Epidemiology*,
  125(2):340-347, 1987.
  <https://doi.org/10.1093/oxfordjournals.aje.a114535>
- Guess HA, Thomas JE. A rapidly converging algorithm for exact binomial
  confidence intervals about the relative risk in follow-up studies with
  stratified incidence-density data. *Epidemiology*, 1(1):75-77, 1990.
  <https://doi.org/10.1097/00001648-199001000-00016>. PMID: 2081244.
- Hirji KF, Mehta CR, Patel NR. Computing Distributions for Exact
  Logistic Regression. *Journal of the American Statistical
  Association*, 82(400):1110-1117, 1987.
  <https://doi.org/10.1080/01621459.1987.10478547>
- ICH E9(R1). *Addendum on Estimands and Sensitivity Analysis in
  Clinical Trials to the Guideline on Statistical Principles for
  Clinical Trials*. International Council for Harmonisation, 2019.
  <https://database.ich.org/sites/default/files/E9-R1_Step4_Guideline_2019_1203.pdf>
- Lehmann EL, Romano JP. *Testing Statistical Hypotheses*, 3rd
  ed. Springer, 2005. ISBN 978-0-387-98864-1.
- Mehta CR, Patel NR. A network algorithm for performing Fisher’s exact
  test in r x c contingency tables. *Journal of the American Statistical
  Association*, 78(382):427-434, 1983.
  <https://doi.org/10.1080/01621459.1983.10477989>
- Mehta CR, Patel NR, Gray R. Computing an exact confidence interval for
  the common odds ratio in several 2x2 contingency tables. *Journal of
  the American Statistical Association*, 80(392):969-973, 1985.
  <https://doi.org/10.2307/2288562>
- Mehta CR, Patel NR. Exact logistic regression: theory and examples.
  *Statistics in Medicine*, 14(19):2143-2160, 1995.
  <https://doi.org/10.1002/sim.4780141908>
- Mehta CR, Patel NR, Senchaudhuri P. Exact Stratified Linear Rank Tests
  for Ordered Categorical and Binary Data. *Journal of Computational and
  Graphical Statistics*, 1(1):21-40, 1992.
  <https://doi.org/10.1080/10618600.1992.10474574>
- SAS Institute. *SAS/STAT User’s Guide*: the `GENMOD` procedure’s
  `EXACT`/`EXACTOPTIONS` statements,
  `METHOD=DIRECT`/`NETWORK`/`NETWORKMC`.
