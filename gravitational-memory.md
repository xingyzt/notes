---
title: Gravitational memory, a brief introduction
date: 2025-05-17
---

_This was originally a term paper I wrote for Prof. Ori Ganor’s general relativity class at UC Berkeley._

The changing mass-energy quadrupole moments in a gravitational scattering or merger event can be decomposed into an oscillatory __AC__ component, and a monotonic __DC__ one.
While the AC component sources the transient gravitational waves that are now heard on a weekly basis,
the DC component is theorized to imprint an additional, permanent __gravitational memory__ on spacetime, which has so far escaped detection.

Zel'dovich, Polnarev, Braginsky, and Grischuk
first discovered gravitational memory in the framework of linearized gravity;
we present a pedagogical derivation of their arguments in [Section 1](#1-gravitational-memory-in-linear-theory) [[Ze74]](#ref-Ze74) [[Br85]](#ref-Br85) [[Br87]](#ref-Br87).
However, 
detailed analyses of nonlinear general relativity 
by Christodoulou, Blanchet, and Damour
in the 1990s uncovered an additional memory contribution from the gravitational wave radiation itself,
which we discuss in [Section 2](#2-nonlinear-memory) [[Ch91]](#ref-Ch91) [[Bl92]](#ref-Bl92).
Astonishingly,
in the case of binary black holes,
the nonlinear term is much larger than the linear memory effect
and of comparable order to the gravitational waves themselves
[[Th92]](#ref-Th92).
Across the frequency bands swept by ground-based and space-based interferometers, 
as well as pulsar timing arrays,
future observational prospects are promising
[[Fa09]](#ref-Fa09) [[Ha10]](#ref-Ha10) [[Ta21]](#ref-Ta21).

On the theoretical side,
gravitational memory also turns out to possess very deep connections to the physics of quantum gravity and black hole information theory.
In the 2010s,
five decades after
the early studies by Bondi, van der Burg, Metzner, and Sachs (BMS)
on the structure of asymptotically-flat spacetime
([Section 3](#3-bms-coordinates)) [[Bo62]](#ref-Bo62) [[Sa62]](#ref-Sa62),
Strominger and Zhiboedov demonstrated the DC effect as a __supertranslation__ on the BMS vacuum
[[St14]](#ref-St14) [[Bo16]](#ref-Bo16)
which we review in [Section 4](#4-supertranslations-and-memory).
Finally,
in [Section 5](#5-a-soft-theorem),
we discuss the ensuing works
which link gravitational memory to Weinberg's soft graviton theorem
[[We65]](#ref-We65) [[St16]](#ref-St16) [[St17]](#ref-St17)
completing the first of many __IR triangles__
in the study of asymptotic gauge symmetries.

_Notation:_
In this review, we set $c = G = 1$;
use mostly-positive metric signatures;
and implicitly sum over repeated, lowercase indices --- Greek for $0,1,2,3$ in spacetime
and Latin for $1,2,3$ in space.
The symbols are most similar to Andrew Strominger's lecture notes [[St17]](#ref-St17).


## 1. Gravitational memory in linear theory

We will first derive gravitational waves and memory as perturbations to the metric,
then examine them in both the linear and nonlinear regimes.

Einstein's equations read
\\[
    R_{\mu\nu} - \f12 R g_{\mu\nu} = 8\pi T_{\mu\nu},
    \tag{1}
\\]
where $R_{\mu\nu}, R$ are the Ricci curvature tensor and trace,
$g_{\mu\nu}$ is the metric,
and $T_{\mu\nu}$ is the stress density.
In linear post-Newtonian theory with the transverse-traceless gauge,
we expand to first order in perturbation from the Minkowskian metric,
$g_{\mu\nu} = \eta_{\mu\nu} + h_{\mu\nu}$,
then enforce the gauge constraints $\Box  x^\mu = 0$ (Lorentz),
$\p^\mu h_{\mu\nu} = 0$ (transverse), and
$\eta^{\mu\nu} h_{\mu\nu} = 0$ (traceless).
This simplifies the equations to
[[Ca97]](#ref-Ca97)
\\[
    \Box h_{\mu\nu}
    = -16\pi T_{\mu\nu}^{TT},
    \tag{2}
\\]

where $^{TT}$ means taking the transverse, traceless component.

Since the nonlinear memory effect is caused by the self-coupling of gravitational radiation,
we would like our results to be Lorentz-covariant.
Following the method employed in relativistic electrodynamics [[Ja98]](#ref-Ja98),
we solve [Eq. 2](#tag-2) by convolution with the retarded Green's function.
The spatial metric perturbation at $x^\mu$
due to a stress-energy-momentum distribution $T_{\mu\nu}(\tilde x^\mu)$ is then

\\[
\tag{3}
    h_{jk}
    (x^\mu)
    =
    \brk{
    8\int \d^4 \tilde x\,
    \Theta(x^0 - \tilde x^0)\,
    \delta((\tilde x^\nu - x^\nu)^2)\,
    T_{jk}(\tilde x^\mu)
    }^{TT},
\\]
where $\Theta$ is the Heaviside step
and $\delta$ is the Dirac delta.

Say we have some piece of matter with rest mass $M$ and rest volume $V$,
moving with four-momentum $P^\mu$ in the observer frame.
If its pressure is negligible, $\mathcal P \ll M/V$,
we can reasonably model its stress density  by

\\[
\tag{4}
    T^{\mu\nu}
    = \f 1V \f{P^\mu P^\nu}M.
\\]

Now,
the characteristic length-scale of gravitational waves,
set by the orbital period,
is typically much larger than the radii of the masses themselves.
As such, we can approximate the number distribution as a thin world-line,
$1/V \simeq \gamma\delta^3 = \int \d\tau\, \delta^4$.
For multiple masses indexed by $A$,
we then have

$$
\tag{6}
    T^{jk}(\tilde x^\mu) = 
    \sum_A 
    \f1{M_A}
    \int \d\tau\,
    \delta^4(\tilde x^\mu - x^\mu_A(\tau)) \,
    P^j_A(\tau) \, P^k_A(\tau)
$$

$$
    \implies
    h^{jk}(x^\mu)
    =
    \brk{
    8
    \sum_A \f1{M_A}
    \int\d\tau\,
    \Theta(x^0 - x_A^0(\tau))\,
    \delta((x^\nu - x_A^\nu(\tau))^2)\,
    P_A^j(\tau)\, P_A^k(\tau)
    }^{TT}.
    \tag{7}
$$

The retarding Heaviside restricts the Dirac delta to spike at the
intersection of the $x_A^\mu$ world-line with $x^\mu$'s past light-cone.
We denote this in mass $A$'s proper time by
$\tau_A(x^\mu; x^\mu_A)$,
implicitly defined such that $x^0_A(\tau_A) < x^0$ and $(x^\mu - x_A^\mu(\tau_A))^2 = 0$.

Since $\delta(f(\tau)) = \delta(\tau-\tau_A)/|f'(\tau-\tau_A)|$ 
for any differentiable function $f$ with a singular zero at $\tau_A$,
the metric perturbation then simplifies to

\\[
\tag{7}
    h^{jk}(x^\mu)
    =
    \brk{
    4
    \sum_A 
    \f{P_A^j(\tau_A)\, P_A^k(\tau_A)}
    {\eta_{\rho\sigma}\,(x^\rho - x_A^\rho(\tau_A)) P^\sigma_A(\tau_A)}\,
    }^{TT}.
\\]

Furthermore,
if the particles are all localized far from the observer,
we can take the average of the line-of-sight 4-vector
$r^\mu \equiv \braket{x^\mu - x_A^\mu(\tau_A)}_A$
and approximate

\\[
\tag{8}
    h^{jk}
    \simeq
    \brk{
    4
    \sum_A 
    \f{P_A^j\, P_A^k}
    {r_\sigma P^\sigma_A}\,
    }^{TT}
\\]

with retarded time $\braket{\tau_A} = x^0 - r^0$.
This is the expression analyzed by Braginsky and Thorne [[Br87]](#ref-Br87).
As it is Lorentz covariant,
we can extend its domain to massless particles.

For oscillatory events such as gravitational inspirals,
the AC component of $h_{jk}$ oscillates
on the order of the orbital period.
However, at least in linear theory,
momenta are negligible in these events' prelude and aftermath,
so there is relatively little DC change to the metric
[[Fa09]](#ref-Fa09).
On the other hand,
the DC effect during recoil events such as 
three-body gravitational scattering and supernova neutrino ejections
can permanently deform the metric by
[[Br87]](#ref-Br87)

\\[
    \Delta h^{jk}
    =
    \Delta\brk{
    4
    \sum_A
    \f{P_A^j P_A^k}{r_\sigma P_A^\sigma}
    }^{TT}.
\tag{9}
\\]

This is what Braginsky & Grischuk termed __gravitational memory__
[[Br85]](#ref-Br85).
Intuitively speaking,
spacetime is not 'elastic' as the electromagnetic fields are;
instead, it behaves more like a gauge potential
---
after the passage of gravitational radiation,
it is completely free to enter a new vacuum configuration.

It is also useful to write the above
in terms of the three-velocity $\vec v_A$,
the line-of-sight distance $r \equiv r^0$,
and their relative angle $\theta_A$,

$$
\tag{10}
    \Delta h^{jk}
    =
    \Delta \brk{
    \f4r
    \sum_A
    \f{P_A^0\, v^j_A\, v^k_A}{1 - v_A\cos\theta_A}
    }^{TT}\!\!\!,\quad
    P_A^0 = \begin{cases}
        \gamma_A M_A, & \text{matter with Lorentz factor }\gamma_A;\\
        E_A, & \text{radiation with energy }E_A.
    \end{cases}
$$

If we align the line-of-sight with the $z$ axis,
the transverse gauge forces all $h_{jk}$ to vanish except for $h_\times\equiv h_{xy} = h_{yx}$ and $h_+\equiv h_{xx} = -h_{yy}$.
Averaging over the velocity directions,
we find the expected root-mean-square change in strain

\\[
\tag{11}
    \Delta h_+^\text{rms}
    = \sqrt3\, \Delta h_\times^\text{rms}
    = \f{\sqrt{8}}r
    \Delta 
    \sum_A P_A^0 v_A^2
    \cdot
    \overbrace{
    \sqrt{
    \f{3}{v_A^4} - \f2{v_A^2}
    - \f{3}{v_A^5\gamma_A^2}\tanh^{-1}v_A
    }}^{\to\sqrt{2/5}\text{ as }v\to0\,;\ \to1\text{ as }v\to1.}.
\\]

Interestingly,
beyond the isotropic Lorentz factor,
the Doppler beaming of ultra-relativistic matter only contributes a $\sqrt{5/2}$ amplification on average.
However,
chance alignments of the velocity vector towards the line-of-sight ($|\theta_A|\ll1$) should significantly dominate observations for events on the edge of detectability.

For non-relativistic matter,
[Eq. 11](#tag-11) reduces to $(8/\sqrt5)(\Delta E_\text{kin}/r)$,
where $E_\text{kin}$ is the matter's total center-of-mass kinetic energy.
Thorne & Braginsky use this expression in their initial order-of-magnitude estimates [[Br87]](#ref-Br87).
They assume a significant fraction of $E_\text{kin}$ is transformed during the recoil interaction,
so that the memory strain is of order $\Delta h_+ \sim (0.1\text{ to }1)E_\text{kin}/r$.
A supernova,
releasing much of its progenitor star's gravitational binding energy $\sim 0.1\, M_\odot$ in seconds via asymmetric neutrino ejection,
would then produce a memory $\sim 10^{-19}$
at a galactic distance of $\sim 10$ kpc;
Whereas a massive black hole $\sim10^5\, M_\odot$,
ejected at relativistic speeds over seconds 
from a merging galactic nucleus,
would produce a memory $\sim 10^{-21}$
 at Hubble distance $\sim 10$ Gpc.

Although ground-based laser interferometers like LIGO 
have maximum strain sensitivities between $10^{-20}\sim10^{-22}$,
they are only designed for high-frequency signals between $10\sim1000$ Hz
--- among many confounding factors,
their mirror pendulums strongly dampen lower-frequency signals.
Conversely,
pulsar timing arrays (PTAs) like NanoGRAV are only sensitive to strain changes on the order of pulsar light travel time
between $10^7\sim10^9$ s,
not to mention they are ‘only’ sensitive to strains above $10^{-15}\sim10^{-17}$.
The best hope for detecting linear memory lies
in space-based laser timing probes such as LISA,
which resolve strains above $10^{-19}\sim10^{-21}$ over timescales of seconds to hours
[[Fa09]](#ref-Fa09) [[Mo15]](#ref-Mo15).

## 2. Nonlinear memory

A codominant contribution was hiding in plain sight:
the change in quadrupole moment due to the gravitational radiation itself.
Not to be confused with the gravitational waves that directly pass by the observer
---
rather,
this strain is an imprint of the momentum anisotropy within the gravitational radiation in general,
$\Delta h_+^\text{rad} \sim \Delta E_\text{rad}/r$.
The immense amount of AC gravitational radiation emitted in the final moments of compact object inspirals,
of order $\Delta E_\text{rad} \sim (0.01\text{ to }0.1)\, M$,
therefore incite DC gravitational memory _of the same order_.
This is rather remarkable ---
as the name __nonlinear memory__ implies,
Chrisdoulou, Blanchet, and Damour only discovered the effect
as they systematically examined wave propagation beyond linear post-Newtonian order,
and encountered the contributions from gravitational wave energy of order $h\p h$
[[Ch91]](#ref-C91) [[Bl92]](#ref-Bl92).

With strains on the same order as AC effects
$\sim10^{-19} (M/M_\odot)(10\text{ kpc}/r)$,
the remaining question for estimating nonlinear memory's observational prospects is timescale.
Focusing on black hole binaries (BBHs),
we can divide their interactions into three phases:
the inspiral from hydrodynamic friction and gravitational radiation,
plunging into a merger
as the orbits become gravitationally unstable,
and concluding with the __ringdown__ of quasi-normal modes
from the newborn black hole [[Bu07]](#ref-Bu07) [[Fa09]](#ref-Fa09).
The merger and ringdown last around $\sim 70M \sim (M/3000M_\odot)$ seconds,
falling into the frequency ranges of both space- and ground-based interferometers,
[[Fa09]](#ref-Fa09).
On the other hand,
memory produced during the inspirals of massive and supermassive BBHs,
which last from minutes to years, may be detectable by space-based laser experiments as well as pulsar timing arrays.
This is particularly exciting for PTAs,
which are typically insensitive to the orbital-period-timescale AC effects [[Fa09]](#ref-Fa09).

In their updated analysis,
Braginsky & Thorne take [Eq. 10](#tag-10)'s discrete sum [[Th92]](#ref-Th92) literally,
and connect this nonlinear effect to the graviton emission terms that they had initially overlooked.
Taking the sum's continuum limit,
we can also interpret the overlooked term classically as the piece of quadrupole moment carried away by the gravitational waves,

\\[
    \Delta h_\text{rad}^{jk}
    = \brk{\f 4r \int \f{\d E_\text{rad}}{\d\Omega }
    \f{\xi^j \xi^k}{1 - \cos\theta}
    \d\Omega
    }^{TT},
\tag{12}
\\]

where $\xi^i$ is a unit three-vector directed from the source to the solid angle $\d\Omega$.
This is the typical starting point for more careful analytical modeling
[[Fa09]](#ref-Fa09).

## 3. BMS coordinates

Bondi, van der Burg, Metzner, and Sachs (BMS) offered an alternative study of gravitational waves from the framework of asymptotically flat spacetimes.
Their picture allows for non-perturbative gravitational phenomena to occur as long as they are localized in the spacetime ‘interior’, 
so that the metric approaches Minkowski spacetime at ‘asymptotic infinity’,
where we position our observers.

Before we impose the asymptotic conditions,
we should first construct a coordinate system $x^\mu$ fit for gravitational wave analysis near future infinity.
Recall that due to the four Bianchi identities,
the Einstein field equations have four gauge degrees of freedom [[Ca97]](#ref-Ca97).
A coordinate system with four gauge constraints can therefore represent any Lorentzian metric.
Following Mitra's pedagogical approach [[Mi23]](#ref-Mi23), 
we then begin by foliating spacetime with the family of future light-cones at spatial origin,
parametrized by the retarded time coordinate $x^0 \equiv u$.
The tangent vector of every null geodesic is thus normal to $n_\mu = -\p_\mu u$.

Demanding each null geodesic to be constant in every coordinate but a radial $x^1 \equiv r$
gives us the first three gauge constraints:
Because future null geodesics then have tangent vectors
$n^\mu = e^{-2\beta} \p_r x^\mu$ for some arbitrary parameterization $\beta(r)$,
we have $g_{r\nu} \to (-e^{2\beta}, 0, 0, 0)$ from $n_\mu = g_{\mu\nu} n^\nu$.
This leaves $x^A$ as angular coordinates on $S^2$, 
where capital Latin indices range $2,3$.
Bondi uses the final gauge freedom to demand flux to fall off as $1/r^2$, 
constraining $\p_r\det(g_{AB}/r^2) = 0$.
With $U \equiv g^{rr}$ and $U^A \equiv g^{rA}$ as arbitrary functions of the coordinates,
we can therefore cast all Lorenztian metrics into the __retarded Bondi gauge__,

\\[
    \d s^2
    = -U\d u^2
    - 2e^{2\beta}\d u\d r
    + g_{AB}\prn{
    \d x^A + \tfrac12 U^A\d u
    }\prn{
    \d x^B + \tfrac12 U^B\d u
    }.
    \tag{13}
\\]

If we foliate with past light-cones instead,
parametrized by the advanced time coordinate $x^{\hat 0} \equiv v$,
we can also write these metrics in the __advanced Bondi gauge__,

\\[
    \d s^2
    = -U\d v^2
    - 2e^{2\beta}\d v\d r
    + g_{AB}\prn{
    \d x^A + \tfrac12 U^A\d v
    }\prn{
    \d x^B + \tfrac12 U^B\d v
    }.
    \tag{14}
\\]


In flat spacetime,
$g_{AB}$ is the standard Euclidean metric on a two-sphere of radius $r$.
Following Strominger and other modern authors,
we choose to explicitly parametrize this by the Riemann sphere 
$x^A \to (z,\bar z)$ due to its coordinate symmetry and holomorphic behavior.
This is equivalent to the more familiar $(\theta,\phi)$ coordinates 
by stereographic projection of $z = \tan(\theta/2)\, \exp(i\phi)$ from the extended complex plane.
The metric on the unit Riemann sphere is given by

\\[
    2\gamma_{z\bar z}\d z\d\bar z
    = \f4{(1 + z\bar z)^2}\d z\d\bar z,
    \tag{15}
\\]

and we compute its differential properties in the [appendix](#appendix-differential-geometry-on-the-riemann-sphere).

Now,
to make precise the notions of ‘interior’ versus an ‘asymptotic infinity’,
we can perform a conformal compactification on this spacetime,

\\[
    u \mapsto \tilde u \equiv \arctan u,\quad 
    v\mapsto \tilde v \equiv \arctan v,
    \tag{16}
\\]

to bring all finite spacetime 
(up to angular coordinates)
into an open diamond,
then realize the asymptotic points as the diamond's boundary.
The future and past null infinities
$\mathcal I^+$, $\mathcal I^-$
then respectively form the hypersurfaces $S^2\times\R$ at $\tilde v = +\pi/2$
and $\tilde u = -\pi/2$.

In their study of gravitational waves,
BMS are specifically interested in asymptotically flat spacetimes created by localized sources.
We could naively enforce this physical constraint 
by defining a metric perturbation relative to the Minkowskian spacetime,
then demanding the metric perturbation to fall-off as $1/r$.
However, this perturbative formulation is obfuscated by the many gauge freedoms,
and at the time BMS even had doubts as to its physical validity [[Bo62]](#ref-Bo62) [[Wa84]](#ref-Wa84).
Instead, modern works on the BMS metric
prefer to state this constraint using the conformally-invariant Weyl tensor

\\[
    C_{\rho\sigma\mu\nu}
    \equiv
    R_{\rho\sigma\mu\nu}
    -
    (g_{\rho[\mu} R_{\nu]\sigma} - g_{\sigma[\mu} R_{\nu]\rho})
    + \f13 R g_{\rho[[\mu} g_{\nu]]\sigma}
    \tag{17}
\\]

so that constraint is rigorous at the conformally-compactified null infinities [[Wa84]](#ref-Wa84) [[Ca97]](#ref-Ca97).
Strominger, for example, proscribes that [[St17]](#ref-St17)

\\[
    C_{rzrz},
    C_{rurz},
    C_{rur\bar z}
    \sim \O(1/r^3).
    \tag{18}
\\]

As the Weyl tensor involves two derivatives of the metric,
this condition does qualitatively agree with demanding a $1/r$ fall-off
in the metric perturbation.

Metrics whose Weyl tensors satisfy this flatness condition have retarded Bondi form

$$\begin{align*}
    \d s^2
    &= -\prn{
    1 - \f{2m}r
    }\d u^2
    - 2\d u\d r
    \notag
    \\
    &\quad
    + r^2\prn{
    \gamma_{z\bar z}\d z\d\bar z
    + \f1r c_{zz}\d z^2
    }
    \notag
    \\
    &\quad 
    + r\prn{
    \f1{r} D^z c_{zz}
    + \f1{r^2}\brk{
    \f43(N_z + u\p_z m)
    - \f14\p_z(c_{zz} c^{zz})
    }
    }\d u\d z
    \notag
    \\
    &\quad
    + \text{c.c.}
    + \O(1/r^3),
    \tag{19}
\end{align*}$$

where $m$, $c_{zz}$, $c^{zz} = c_{\bar z\bar z}/(\gamma_{z\bar z})^2$, and $N_z$ are functions only of $u$, $z$, and $\bar z$ [[St17]](#ref-St17).
This form is very alien to those of us trained in linear theory.
To help see the connection, 
we can switch to a timelike coordinate

\\[
    t \equiv u + r
    + m\log\f{m}r
    - \prn{1 - \f{m} r}^{-1} D^z c_{zz} z
    + \text{c.c.}
     + \O\prn{\f1{r^2}},
     \tag{20}
\\]

then take the Cartesian flat-sky limit around the north pole of the Riemann sphere,

\\[
    z \sim \f{x + iy}{r\sqrt2} + \O\prn{\f1{r^2}}.
    \tag{21}
\\]

Now, we recover a more familiar expression at large $r$,

$$ \begin{align*}
    \d s^2
    &= -\prn{
    1 - \f{2m}r
    }\d t^2
    + \prn{
    1 + \f{2m}r
    }\d r^2
    \notag
    \\
    &\quad
    + \prn{
    1 + \f{c_{zz} + c_{\bar z\bar z}}{2r}
    }\d x^2
    + \prn{1
    - \f{c_{zz} + c_{\bar z\bar z}}{2r}
    }\d y^2
    \notag
    \\
    &\quad
    - \f{c_{zz} - c_{\bar z\bar z}}{ir}
    \d x\d y
    + \O\prn{\f1{r^2}}.
    \tag{22}
\end{align*}$$

Directly reading this against the metric of linear theory,
we can identify $m$ as the mass term,
and changes in $(c_{zz} + c_{\bar z\bar z})/{2r}$, $(c_{zz} - c_{\bar z\bar z})/{2ir}$ 
as the linear gravitational wave modes.
In fact, $m$, termed the __Bondi mass aspect__,
does behave like mass in Kerr-Newmann metrics;
While the retarded time derivatives of $c_{zz}$ and $c_{\bar z\bar z}$,
termed the __Bondi news tensor__ $N_{AB} \equiv \p_u c_{AB}$,
correspond to the two circular gravitational wave modes [[St17]](#ref-St17).

Notice that, just as in the post-Newtonian linear theory
([Section 1](#1-gravitational-memory-in-linear-theory)),
expanding to first order in $1/r$ misses the gravitational radiation's self-interaction $\p_z(c_{zz} c^{zz})$,
as well as higher moments $N_z$, $\p_z m$, $D^z c_{zz}$.
Nonlinear gravitational memory therefore only reveals itself at $\O(1/r^2)$.
It would be an interesting exercise to see the analogues to all these terms in 
the second-order post-Newtonian expansion.
Unfortunately, as I am short on time,
it suffices to say that $N_A$ corresponds to angular momentum [[St17]](#ref-St17).

Because $m$, $c_{AB}$, and $N_A$ are all theoretically measurable using inertial test masses at null infinity,
and moreover can combine to describe any metric in the Bondi gauges which satisfies the BMS asymptotic flatness criterion up to $\O(1/r^3)$,
they are together known as the metric's __asymptotic data__ up to $\O(1/r^3)$
[[St14]](#ref-St14) [[Bo16]](#ref-Bo16).
Spacetimes with different asymptotic data are therefore physically inequivalent.
However, the arrival of gravitational radiation (in the form of ‘Bondi news’) can change the asymptotic data,
and physically shift us from one spacetime to another
[[Bo62]](#ref-Bo62).
Furthermore,
as long as these news preserve the metric's asymptotic flatness,
the spacetime will continue to satisfy the BMS constraints.
Thus,
even though the Bondi news have physically changed the physical situation,
they still represent,
in a sense,
symmetries of BMS spacetimes.

## 4. Supertranslations and memory

Typically, when we are examining the symmetries of a spacetime,
we focus on isometries --- diffeomorphisms which leave the metric globally invariant.
However,
since our role in the BMS metric ([Eq. 19](#tag-19)) are only as observers at infinity,
so we can relax the isometry requirement and examine diffeomorphisms which,
though not necessarily symmetries in the spacetime interior,
are still __asymptotic symmetries__ in that
they preserve the metric's flatness at infinity.
As they are not necessarily isometries,
there is no demand that they leave the physical system invariant.
The group of asymptotic symmetries of the BMS metric,
quotiented up to Bondi gauge invariance and other physically trivial transformations,
form the __BMS group__ [[Sa62]](#ref-Sa62) [[Bo16]](#ref-Bo16) [[St14]](#ref-St14).
It exhibits far greater complexity
than the Poincar\'e group's action on globally flat spacetime.
This strange behavior is key to revealing gravitational memory's
connection to quantum information theory.

Every diffeomorphism is infinitesimally generated by some vector field $\xi^\nu(x^\mu)$.
We can derive $\xi$ by relating its flow of the metric,
given by the Lie derivative [[Ca97]](#ref-Ca97) [[St17]](#ref-St17)

\\[
    \mathcal L_\xi g_{\mu\nu} 
    = 2\nabla_{(\mu} g_{\nu)\rho} \xi^\rho
    + \xi^\rho \nabla_\rho g_{\mu\nu},
    = 2\p_{(\mu} g_{\nu)\rho} \xi^\rho
    + \xi^\rho \p_\rho g_{\mu\nu},
    \tag{23}
\\]

(where $\nabla_\mu$ is the covariant derivative)
to the constraints we are imposing on the metric.
Following Strominger [[St17]](#ref-St17),
we will first limit restrict ourselves to irrotational,
translation-like transformations.
This lets us ignore the $\O(1/r^2)$ terms in the metric,
as well as set $\xi^u, \xi^r, r\xi^z$, and $r\xi^{\bar z}$ all to be constant in $r$.

Under these restrictions, 
$g_{ur}$, $g_{rz}$, and $g_{z\bar z}$ should be invariant under $\xi$ flow.
Thus,

\\[
    \mathcal L_\xi g_{ur} = -\p_u \xi^u + \O(1/r)
    \tag{24}
\\]

means to leading order,
$\xi^u$ is purely some function $f$ of the angular coordinates $(z,\bar z)$.
Then,

$$\begin{align*}
\mathcal L_\xi g_{rz}
&=
-\p_z \xi^u  + r\gamma_{z\bar z}(2\xi^{\bar z} + r\p_r \xi^{\bar z})
+ \O(1/r)
\\
&=
r\gamma_{z\bar z}\prn{\f1r D^{\bar z} f + 2\xi^{\bar z} + r\p_r \xi^{\bar z}}
+ \O\prn{\f1r}
\tag{25}
\end{align*}$$

is solved by

\\[
    \xi^{\bar z} = -\f1r D^{\bar z} f + \O\prn{\f1{r^2}},
    \tag{26}
\\]

and likewise for $\xi^z$.
Finally,

$$
\begin{align*}
\mathcal L_\xi g_{z\bar z}
&= r\gamma_{z\bar z} (2\xi^r + r D_z \xi^z + r D_{\bar z}\xi^{\bar z}) 
+ \O(1/r)
\\
&= r\gamma_{z\bar z} (2\xi^r + D\cdot D f)
+ \O(1/r),
\end{align*}
\tag{27}
$$

with $D\cdot D$ as the spherical Laplacian (see the [appendix](#appendix-differential-geometry-on-the-riemann-sphere)),
is solved by $\xi^r = \frac12(D\cdot D) f + \O(1/r)$.
We have thus found the generators of ‘supertranslations’,

\\[
	\xi(f) = f\p_u - \f1r(D^zf\p_z + D^{\bar z}f\p_{\bar z}) + \f12(D\cdot D) f \p_r
    + \O\prn{\f1r}.
    \tag{28}
\\]

Unlike the three-dimensional family of translational isometries in flat spacetime,
these supertranslations are parametrized by an infinite-dimensional family of functions $f$ on the sphere;
and perform transformations between physically inequivalent spacetimes:
Asymptotic data is not preserved, but to order $\O(1/r)$ flow under $\xi(f)$ as

$$\begin{align*}
\mathcal L_f N_{zz} &= f \p_u N_{zz},\\
	\mathcal L_f m &= f \p_u m + \f18\brk{ N^{zz} (D\cdot D) f + 4 D_z N^{zz} D_z f + \text{c.c.} },\\
    \tag{29}
	\mathcal L_f c_{zz} &= f \p_u c_{zz} - (D\cdot D) f,
\end{align*}$$

To model the passage of a transient gravitational wave,
Strominger & Zhiboedov suppose that outside some time interval $u_i < u < u_f$ 
at null future infinity $\mathcal I^+$,
we are in some __BMS vacuum__ free of gravitational radiation ---
equivalently, with vanishing Bondi news $N_{AB}$.
Earlier, BMS had shown that the lack of news implies the conservation of Bondi mass $m$,
and that BMS vacuums are related by supertranslations
[[Bo62]](#ref-Bo62).
Building upon their work, Strominger & Zhiboedov demonstrate how
the arrival of Bondi news then heralds the transition from one vacuum to another,
and imprints any permanent change it brings in the perturbations $c_{AB}$ as gravitational memory
[[St16]](#ref-St16),
In the framework of BMS,
we can therefore view gravitational memory as an asymptotic vacuum transition.


## 5. A soft theorem

Now we turn to quantum field theory for another surprise.
In the 60s,
Weinberg proved that scattering processes are measurably affected by their interaction
with undetectably low-energy gauge bosons,
even if these bosons are, in the ‘soft limit’, taken to have zero momenta
[[We65]](#ref-We65).
More formally,
whenever we consider some process with matrix element $i\mathcal M_0$,
we need to also consider modifications where we attach soft bosons
to arbitrary real legs in its Feynman diagram
[[Sc13]](#ref-Sc13).

Suppose a leg initially had a particle with four-momenta $p^\mu$ entering the process;
then we modify it to have incoming momenta $p^\mu + q^\mu$,
transferring some momentum $q$ to a gauge boson interaction vertex,
and entering the process with momentum $p^\mu$.
In the case of a graviton with polarization $\epsilon_{\mu\nu}$,
the Feynman rule for its interaction vertex would need to be some rank-2 tensor $i\Gamma^{\mu\nu}$
built covariantly from $p^\mu$ and $q^\mu$.
Generically then,

\\[
	-i \Gamma^{\mu\nu} = p^\mu p^\nu F_1 + p^\mu q^\nu F_2 + q^\mu q^\nu F_3,
    \tag{30}
\\]

where the __form factors__ $F_n$ 
are functions of the invariant quantities $q_\rho p^\rho$ and $q_\rho q^\rho$.
Combined with the on-shell graviton propagator,
the modified process would read

\\[
	i\mathcal M_0 \mapsto
	i\mathcal M_0' = 
	-i \frac{\Gamma^{\mu\nu} \epsilon_{\mu\nu} }{(p^\rho - q^\rho)^2 + m^2} \mathcal M_0.
    \tag{31}
\\]

By the Ward identity,
the graviton four-momentum is transverse to its polarization,
$q^\mu \epsilon_{\mu\nu} = 0$,
so the $F_2$ and $F_3$ form factors vanish.
Substituting $p_\rho p^\rho = -m^2$ and taking the soft limit $q_\rho q^\rho\to 0$,
we then obtain

\\[
	i\mathcal M_0' = -i F_1(0,0)\, \f{p^\mu p^\nu \epsilon_{\mu\nu} }{q_\rho p^\rho} \mathcal M_0,
    \tag{32}
\\]

This form is suspiciously similar to what we obtained for gravitational memory ([Eq. 2](#tag-2)),
with $q_\rho$ serving as our line-of-sight null vector.
Strominger goes on to analyze the quantum mechanics of this connection in great depths
[[St17]](#ref-St17).
Unfortunately, as we are again short on time,
let us content ourselves with the qualitative conclusion that, in the infrared regime of general relativity,
stars behave much the same as elementary particles.


## 6. Outlook

From its humble start as a curious astrophysical phenomenon,
the study of gravitational memory has blossomed into a rich tangle of topics
that seemingly has something interesting prepared for everyone.
For mathematical physicists,
the existence of soft theorems for other gauge theories --- notably QED and QCD ---
have inspired searches for similar asymptotic symmetries and memory effects.
Within classical relativity,
Pasterski, Strominger, Zhiboedov have also found a higher-order __gravitational spin memory__
associated with angular momentum conservation ---
the puzzling nature of angular momentum in general relativity makes this the more intriguing.
At the same time,
new ideas in the BMS formalism is still yielding improvements to the numerical relativity simulations
that serve as the backbone for the gravitational wave observatories.
One key prediction from these numerical relativists
--- 
a sharp rise in gravitational memory at the final moments of merger and ringdown
---
may be confirmed in the coming decades by space-based laser interferometers.


## Appendix: Differential geometry on the Riemann sphere

As with any Euclidean metric on the unit 2-sphere,
the Riemann sphere metric $\gamma$ ([Eq. 15](#tag-15)) equals its own Ricci curvature tensor.
It has non-vanishing Christoffel symbols 
$\Gamma^z_{zz} = \overline{\Gamma^{\bar z}_{\bar z\bar z}} = -2\bar z/(1 + z\bar z)$,
so its covariant derivative acts on vectors by

$$\begin{align*}
    D_z V^z &= \prn{\p_z - \f{2\bar z}{1 + z\bar z}} V^z,
    & D_{\bar z} V^{\bar z} &= \p_{\bar z} V^{\bar z},  & \text{(c.c.)}.
\end{align*}$$

Its Laplacian also has the property that

\\[
    \f12D\cdot D 
    = D_z D^z
    = \gamma_{z\bar z} D^{\bar z} D^z
    = D_{\bar z} D^{\bar z}
    = \f1{\gamma_{z\bar z}}\p_z \p_{\bar z}.
\\]

An interesting exercise would be finding the spherical harmonics in these coordinates.

## References

1. <span id="ref-Ze74"> Y. B. Zel’dovich and A. G. Polnarev. </span> “Radiation of gravitational waves by a cluster of superdense stars”. In: _Sov. Astron._ 18 (1974), p. 17.

1. <span id="ref-Br85"> V. B. Braginsky and L. P. Grishchuk. </span> “Kinematic Resonance and Memory Effect in Free Mass Gravitational Antennas”. In: _Sov. Phys. JETP_ 62 (1985), pp. 427–430. 

1. <span id="ref-Br87"> V. B. Braginsky and Kip S. Thorne. </span> “Gravitational-wave bursts with memory and experimental prospects”. In: _Nature_ 327.6118 (May 1987), pp. 123–125. doi: [10.1038/327123a0](https://doi.org/10.1038/327123a0).

1. <span id="ref-Ch91"> Demetrios Christodoulou. </span> “Nonlinear nature of gravitation and gravitational-wave experiments”. In: _Phys. Rev. Lett._ 67 (12 Sept. 1991), pp. 1486–1489. doi: [10.1103/PhysRevLett.67.1486](https://doi.org/10.1103/PhysRevLett.67.1486).

1. <span id="ref-Bl92"> Luc Blanchet and Thibault Damour. </span> “Hereditary effects in gravitational radiation”. In: _Phys. Rev. D_ 46 (1992), pp. 4304–4319. doi: [10.1103/PhysRevD.46.4304](https://doi.org/10.1103/PhysRevD.46.4304).

1. <span id="ref-Th92"> Kip S. Thorne. </span> “Gravitational-wave bursts with memory: The Christodoulou effect”. In: _Phys. Rev. D_ 45 (2 Jan. 1992), pp. 520–524. doi: [10.1103/PhysRevD.45.520](https://doi.org/10.1103/PhysRevD.45.520).

1. <span id="ref-Fa09"> Marc Favata. </span> “Gravitational-wave memory revisited: memory from the merger and recoil of binary black holes”. In: _J. Phys. Conf. Ser._ 154 (2009). Ed. by Alberto Lobo and Carlos F. Sopuerta, p. 012043. doi: [10.1088/1742-6596/154/1/012043](https://doi.org/10.1088/0004-637X/696/2/L159). arXiv: [0811.3451 [astro-ph]](https://arxiv.org/abs/0902.3660).

1. <span id="ref-Ha10"> Rutger van Haasteren and Yuri Levin. </span> “Gravitational-wave memory and pulsar timing arrays”. In: _Mon. Not. Roy. Astron. Soc._ 401 (2010), p. 2372. doi: [10.1111/j.1365-2966.2009.15885.x](https://doi.org/10.1111/j.1365-2966.2009.15885.x). arXiv: [0909.0954 [astro-ph.IM]](https://arxiv.org/abs/0909.0954).

1. <span id="ref-Ta21"> Stephen R. Taylor. </span> _The Nanohertz Gravitational Wave Astronomer._ Taylor & Francis, May 2021. arXiv: [2105.13270 [astro-ph.HE]](https://arxiv.org/abs/2105.13270).

1. <span id="ref-Bo62"> H. Bondi, M. G. J. van der Burg, and A. W. K. Metzner. </span> “Gravitational Waves in General Relativity. VII. Waves from Axi-Symmetric Isolated Systems”. In: _Proceedings of the Royal Society of London Series A_ 269.1336 (Aug. 1962), pp. 21–52. doi: [10.1098/rspa.1962.0161](https://doi.org/10.1098/rspa.1962.0161).

1. <span id="ref-Sa62"> R. Sachs. </span> “Asymptotic symmetries in gravitational theory”. In: _Phys. Rev._ 128 (1962), pp. 2851–2864. doi: [10.1103/PhysRev.128.2851](https://doi.org/10.1103/PhysRev.128.2851).

1. <span id="ref-St14"> Andrew Strominger. </span> “On BMS Invariance of Gravitational Scattering”. In: _JHEP_ 07 (2014), p. 152. doi: [10.1007/JHEP07(2014)152](https://doi.org/10.1007/JHEP07(2014)152). arXiv: [1312.2229 [hep-th]](https://arxiv.org/abs/1312.2229).

1. <span id="ref-Bo16"> Michael Boyle. </span> “Transformations of asymptotic gravitational-wave data”. In: _Phys. Rev. D_ 93.8 (2016), p. 084031. doi: [10.1103/PhysRevD.93.084031](https://doi.org/10.1103/PhysRevD.93.084031). arXiv: [1509.00862 [gr-qc]](https://arxiv.org/abs/1509.00862).

1. <span id="ref-We65"> Steven Weinberg. </span> “Infrared Photons and Gravitons”. In: _Phys. Rev._ 140 (2B Oct. 1965), B516–B524. doi: [10.1103/PhysRev.140.B516](https://doi.org/10.1103/PhysRev.140.B516).

1. <span id="ref-St16"> Andrew Strominger and Alexander Zhiboedov. </span> “Gravitational Memory, BMS Super- translations and Soft Theorems”. In: _JHEP_ 01 (2016), p. 086. doi: [10.1007/JHEP01(2016)086](https://doi.org/10.1007/JHEP01(2016)086). arXiv: [1411.5745 [hep-th]](https://arxiv.org/abs/1411.5745).

1. <span id="ref-St17"> Andrew Strominger. </span> _Lectures on the Infrared Structure of Gravity and Gauge Theory._ Mar. 2017. isbn: 978-0-691-17973-5. arXiv: [1703.05448 [hep-th]](https://arxiv.org/abs/1703.05448).

1. <span id="ref-Ca97"> Sean M. Carroll. </span> _Lecture notes on general relativity._ Dec. 1997. Chap. 3–6. arXiv: [gr-qc/9712019](https://arxiv.org/abs/gr-qc/9712019).

1. <span id="ref-Ja98"> John David Jackson. </span> _Classical Electrodynamics._ Wiley, 1998. Chap. 12 and 14. isbn: 978-0-471-30932-1.

1. <span id="ref-Mo15"> C. J. Moore, R. H. Cole, and C. P. L. Berry.</span> “Gravitational-wave sensitivity curves”. In: _Class. Quant. Grav._ 32.1 (2015), p. 015014. doi: [10.1088/0264-9381/32/1/015014](https://doi.org/10.1088/0264-9381/32/1/015014). arXiv: [1408.0740 [gr-qc]](https://arxiv.org/abs/1408.0740). url: [https://www.sr.bham.ac.uk/˜cplb/GWplotter/](https://www.sr.bham.ac.uk/~cplb/GWplotter/).

1. <span id="ref-Bu07"> Alessandra Buonanno, Gregory B. Cook, and Frans Pretorius. “Inspiral, merger and ring-down of equal-mass black-hole binaries”. In: _Phys. Rev. D_ 75 (2007), p. 124018. doi: [10.1103/PhysRevD.75.124018](https://doi.org/10.1103/PhysRevD.75.124018). arXiv: [gr-qc/0610122](https://arxiv.org/abs/gr-qc/0610122).

1. <span id="ref-Mi23"> Prahar Mitra. </span> “Derivation of Bondi metric”. Physics Stack Exchange. (version: 2023-06- 22). url: [https://physics.stackexchange.com/q/636961](https://physics.stackexchange.com/q/636961).

1. <span id="ref-Wa84"> Robert M. Wald. </span> General Relativity. Chicago, USA: Chicago Univ. Pr., 1984. Chap. 11. doi: [10.7208/chicago/9780226870373.001.0001](https://doi.org/10.7208/chicago/9780226870373.001.0001).

1. <span id="ref-Sc13"> Matthew D. Schwartz. </span> _Quantum Field Theory and the Standard Model_. Cambridge University Press, 2013. Chap. 9.5. doi: [10.1017/9781139540940](https://doi.org/10.1017/9781139540940).
