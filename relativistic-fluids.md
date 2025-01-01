---
title:  "Relativistic fluids from relativistic dust"
---

This article demonstrates a method to construct the energy-momentum tensor
of a perfect fluid from “dust”, and aims to give an intuition of this tensor's meaning.
It builds upon the lectures of Professor Ori Ganor.
*[See my other notes here](/notes)*. 

## Background

*(This section is included as a brief review for those familiar with the basics of special relativity.
For a much more thorough introduction, I highly recommend 
[Eigenchris’s YouTube series on relativity](https://www.youtube.com/watch?v=bEtBncTEc6k&list=PLJHszsWbB6hqlw73QjgZcFh4DrkQLSCQa)
and
[tensor algebra](https://www.youtube.com/watch?v=8ptMTLzV4-I&list=PLJHszsWbB6hrkmmq57lX8BV-o-YIOFsiG),
as well as the first part of
[Sean Carroll’s lecture notes on general relativity](https://arxiv.org/abs/gr-qc/9712019).
Conversely, feel free to 
[skip this section](#derivations)
if you want to get right to the derivations.)*

A key postulate of special relativity is

### 1. The speed of light in vacuum is constant in every reference frame.

In other words, $c$ is **invariant**. 
This motivates the definition of invariant **spacetime intervals**:
\\[
|\d s|^2 \equiv -|c\,\d t|^2 + |\d x|^2 + |\d y|^2 + |\d z|^2.
\\]
Just like how the lengths of everyday objects do not change under spatial rotations,
the spacetime interval between events do not change under Lorentz boosts.
Even as strange effects like length-contraction and time-dilation become visible when objects are 
moving relative to each other near the speed of light, the spacetime interval stays constant.
It’s like having a new and improved, *relativistic* ruler.

To produce the invariant spacetime interval with this ruler, it must be able to measure time alongside the spatial lengths.
Geometrically, this motivates adding a zeroth $t$ element to our vectors’ existing three $x$, $y$, and $z$ components.
For example, a position vector is now

$$ \vec r = r^\alpha \vec e_\alpha, $$

where $\vec e_\alpha \to (\uvec t\ \uvec x\ \uvec y\ \uvec z)$ are the temporal and spatial **basis vectors**.
(To simplify the writing, Greek-letter indices like $\alpha, \beta, \mu, \nu$ implicitly range over 0, 1, 2, 3;
Repeated upper-lower index pairs like ${^\beta}{_\beta}$ are implicitly summed over. 
Do not confuse them with exponential powers, which will always be outside vertical brackets $|v|^2$.
This notation is known as the **Einstein summation convention**.)

This way, a spacetime displacement
\\[
\d\vec s = s^\alpha \vec e_\alpha \equiv \d t\,\vec e_0 + \d x\,\vec e_1 + \d y\,\vec e_2 + \d z\,\vec e_3
\\]
with components 
\\[
\d s^\alpha \to 
\begin{pmatrix} \d s^0 \\ \d s^1 \\ \d s^2 \\ \d s^3 \end{pmatrix} \equiv
\begin{pmatrix} \d t \\ \d x \\ \d y \\ \d z \end{pmatrix}
\\]
has the spacetime interval (again, summing over repeated index pairs)
\\[
(\d\vec s)^2 =  g_{\alpha\beta} \d s^\alpha \d s^\beta,
\\]
where
\\[ 
\vec g =  g_{\alpha\beta} \ (\vec e^\alpha \otimes \vec e^\beta)
\\]
is the **Minkowskian metric** with components
\\[
g_{\alpha\beta} \to \begin{pmatrix}
-1 \\\ & 1 \\\ & & 1 \\\ & & & 1
\end{pmatrix}.
\\]
(When dealing with objects having two indices, the left index corresponds to the row, 
and the right index the column. For instance, because the metric is symmetric, $ g_{\alpha\beta} =  g_{\beta\alpha}$.)

The metric is useful for a lot of things (see the links above for more details), but most fundamentally, 
it underscores *an interweaving of space and time* that is thematic throughout
special and general relativity. This means

### 2. To gain full intuition, we must treat space and time on an equal footing.

I strongly believe that to fully understand special relativity, 
it must feel natural to mix quantities derived from space and those derived from time.
For example:
* Velocity is both motion through space and motion through time.
  In your **inertial reference frame** $\mathcal K$ (also known as your **rest frame**), 
  you are stationary in space, but still traveling in time.
  Your **four-velocity** $\vec u$ thus has rest-frame components
\\[
u^\alpha \to \begin{pmatrix}1\\\0\\\0\\\0\end{pmatrix}.
\\]

* Likewise, momentum is the motion of mass-energy through space and through time.
  A baseball flying through space can pack a punch with its spatial momentum components,
  but sitting on a table, it is still carrying all of its mass-energy in its flight through time.
  In the baseball’s rest frame, its **four-momentum** is thus
\\[
p^\alpha = mu^\alpha \to \begin{pmatrix}m\\\0\\\0\\\0\end{pmatrix}.
\\]
*(We can create new **four-vectors** by multiplying them with scalar invariants.
Here, $m$ always refers to “the mass as measured in the rest frame,” which is by-definition invariant.)*

* When we multiply this baseball’s four-momentum by its rest-frame number density $n$,
  we obtain its **four-momentum flux**. In the rest frame, this is just
\\[
n p^\alpha \to \begin{pmatrix}n m\\\0\\\0\\\0\end{pmatrix}.
\\]

  In other words, the time component of four-momentum flux is mass-energy density!

<br>

These intertwined interpretations are key to understanding the stress-energy-momentum tensor. More on that in a bit.

For now, notice that we have only defined these quantities in the rest frames of their respective objects.
With all the length-contraction and time-dilation stuff, how do we calculate how these quantities might look 
when we are moving relative to these objects near the speed of light? Luckily,

### 3. We can transform between any inertial reference frames via Lorentz boosts

Suppose another inertial reference frame $\mathcal K'$ is traveling with (ordinary) velocity $v\vec e_1$ 
as measured in $\mathcal K$.
Given some four-vector $\vec a$, it will have a set of components $a^\alpha$ in $\mathcal K$, 
and another set $a^{\alpha^\prime}$ in $\mathcal K'$.
We can try to relate them by a **Lorentz boost** $\Lambda$:
\\[
a^{\alpha^\prime} = \Lambda{^{\alpha^\prime}}_\alpha a^\alpha
\\]
(Notice how in the primed reference frame $\mathcal K'$,
the indices are annotated with $'$ as well.)

Here's what we know about $\Lambda$:
- It must be a *linear transformation* on the components of $\vec a$,
  since we expect adding or scaling vectors are actions independent of any choice of reference frames.
- It must *leave the metric components unchanged*,
  since the speed of light is the same in every reference frame.
- It must be *invertible by the reverse boost*
  (in other words, $(\Lambda^{-1}){^{\alpha}}_{\alpha^\prime} = \Lambda{\_{\alpha^\prime}}^\alpha$),
  since boosting from $\mathcal K$ to $\mathcal K'$ then back to $\mathcal K'$ again should not change anything.

If you work through the algebra,
these conditions uniquely constrain the components of the Lorentz boost by $v$ in the $\vec e\_1$ direction to be
\\[
\begin{align\*}
    \Lambda{^{\alpha^\prime}}\_\alpha
    \equiv \frac{\d x^{\alpha^\prime}}{\d x^\alpha} 
    \to \begin{pmatrix} \gamma & -\gamma v\\\ -\gamma v & \gamma\\\ & & 1\\\ & & & 1 \end{pmatrix},\quad
    \gamma \equiv \f1{\sqrt{1 - |v|^2}}.
\end{align\*}
\\]

Try acting this on some four-vectors. 
Notice how it mixes components of time ($\vec e\_0$) and space ($\vec e\_1$)!

Similarly, the Lorentz boost when $\mathcal K'$ is moving in the $\vec e_2$ direction relative to $\mathcal K$ is
\\[\begin{align\*}
    \Lambda{^{\alpha^\prime}}\_\alpha
    \to \begin{pmatrix} \gamma & & -\gamma v\\\ & 1\\\ -\gamma v & & \gamma\\\ & & & 1 \end{pmatrix}
\end{align*}.\\]

Can you figure out what it should look like when $\mathcal K'$ is moving in the $\vec e_3$ direction?


### Summary of the notation

* We work in natural units, where the speed of light $c=1$. 
* The metric $\vec g$ is spacelike: 
a spacetime displacement $\d\vec s$ has a positive interval $(\d\vec s)^2$
if its spatial displacement exceeds its temporal displacement in magnitude.
* Greek indices ($\alpha,\beta,\mu,\nu,\cdots$) ranging 0, 1, 2, 3 cover the four spacetime dimensions.
They are annotated according to the reference frame.
* Contravariant components have upper indices; covariant components have lower indices.
* Implicit Einstein summation convention is assumed for matching upper-lower indices.

With these things in mind, let’s get dusty!

## Intuitions

A fluid is a macroscopic phenomenon composed of particles 
--- packets of mass-energy --- 
traveling at significant speeds relative to neighboring particles.
It has emergent properties, such as density and pressure, which are not well-defined for its individual constituents.
Examples include water, neutron star interiors, and photons scattering in a sparse medium (but not, for instance, photons concentrated in a laser beam).

In particular, we are interested in the behavior of perfect fluids, which are isotropic in their rest frame, and have 
negligible interactions (except some scattering to keep them isotropic).

Fluid volumes are locally modeled by a stress-energy-momentum tensor $\vec T$, 
whose components qualitatively mean
\\[
T^{\alpha\beta} \to \begin{pmatrix} \text{mass-energy density} 
& \text{mass-energy flux flowing next to }\beta \\\ \text{momentum density in }\alpha 
& \text{momentum flux in }\alpha\text{ flowing next to }\beta \end{pmatrix}.
\\]

This is a four-by-four, **rank-two tensor** (four-vectors are rank-one tensors),
which we can build out of products and sums of scalars and four-vectors.
At first glance, it looks rather convoluted,
but that's only because we aren't treating space and time on an equal footing yet.
When we start thinking of density as the time component of flux,
and mass-energy as the time component of four-momentum, 
all of this reduces to
\\[
T^{\alpha\beta} \to \prn{\text{momentum flux in }\alpha\text{ flowing next to }\beta}.
\\]

In mechanics class,
you may have stumbled upon the fact that pressure, shear, and energy density all have the same dimensions.
In fact, all of them are in this tensor!
The diagonal components, representing the $\alpha$-component of four-momentum flowing next to the $\alpha$ normal vector,
are nothing but pressure, created by particles transmitting momentum directly against the fluid volume’s surface.
While the off-diagonal components, representing the components of four-momentum flowing in directions $\alpha$ orthogonal to the $\beta$ normal vector, 
are nothing but shear, transmitting force parallel to the surface.
This is why we need to describe stress-energy-momentum as a rank-two tensor:
it measures the interaction of *two* directed quantities, 
while a rank-one four-vector only measures one direction.

Still, there seems to be an asymmetry between the $\alpha$ and $\beta$ directions here.
But this inbalance is also illusory.
If we consider fluid elements to always hold the same amount of rest mass,
then any flux in $\alpha$ going to $\beta$ must equal any flux in $\beta$ going to $\alpha$,
lest they accumulate.
Thus $T^{\alpha\beta} = T^{\beta\alpha}$: the stress-energy-momentum tensor is fully symmetric!

## Derivations

We want to find the stress-energy-momentum tensor for a small volume $\Delta V$ of a perfect fluid.

Because extensive quantities of the same inertial reference frame are additive, 
our method of attack is as follows: 
1. Find the stress-energy-momentum tensor of each constituent particle in their respective rest frames.
2. Lorentz boost every constituent tensor to the rest frame of the fluid volume.
3. Make their components extensive, then add them up.

Let’s do this!

### 1. The building block: stress-energy-momentum tensor for dust

As noted in the previous section, fluids are composed of many particles,
each with their own stress-energy-momentum tensor.
Since we are using letter indices to denote Lorentz geometry already,
we will use the $\\#$ symbol to index these particles.

In the inertial reference frame $\mathcal K^\\#$ of some particle $\\#$,
it is completely stationary, 
so all flux components in its stress-energy-momentum tensor
$T_{\\#}^{\alpha^\\#\beta^\\#}$ vanishes.
(We annotate indices in $\mathcal K^\\#$ with $^\\#$.)
Because it exerts no intrinsic pressure on its surroundings, we call this particle **dust**.
Only the mass-energy density component is nonzero:
$T_\\#^{0^\\#0^\\#} = \rho_\\# = m_\\# / V_\\#$,
the particle’s rest mass divided by its rest volume.
Only when we transform to some other reference frame $\mathcal K$ 
do we see the momentum flux components in $T_\\#^{\alpha\beta}$ emerge.

How can we find out how the components of $\vec T_\\#$ transform?
Well, we know how four-vectors transform, 
so we should try to write this stress-energy-momentum tensor in terms of those things.
A good candidate is the particle’s four-velocity $\vec u_\\#$,
whose components are $(1\ 0\ 0\ 0)$ in its rest frame.
Since only the time-time component of $\vec T_\\#$ is nonzero in its rest frame,
we can therefore write
\\[
T\_\\#^{\alpha^\\#\beta^\\#} = \rho\_\\#\, u\_\\#^{\alpha^\\#} u\_\\#^{\beta^\\#}
\to \begin{pmatrix} \rho\_\\#\\\ & 0\\\ & & 0\\\ & & & 0 \end{pmatrix},
\\]
or equivalently,
\\[
\vec T_\\# = \rho_\\# \vec u_\\#\otimes\vec u_\\#.
\\]

### 2. Putting dust in the fluid frame

If particle $\\#$ moves with three-velocity $v_\\#\vec e_1$ in the rest frame of its fluid element $\mathcal K$,
then the components of its four-velocity $\vec u_\\#$ in $\mathcal K$ are boosted to

$$\begin{align*}
u_\#^\alpha
&= \Lambda{^{\alpha}}_{\alpha^\#} u_\#^{\alpha^\#}
\to \begin{pmatrix}
\gamma & \gamma v\\
\gamma v & \gamma\\
& & 1\\
& & & 1
\end{pmatrix}
\begin{pmatrix}
1\\ 0\\ 0\\ 0
\end{pmatrix}
= \gamma_\# \begin{pmatrix}
1\\ v_\#\\ 0\\ 0
\end{pmatrix}.
\end{align*}$$

It is too tedious to derive here,
but for an arbitrary relative three-velocity $v\_\\#^1\vec e\_1 + v\_\\#^2 \vec e\_2 + v\_\\#^3\vec e_3$
with total three-magnitude $v_\\#$,
the four-velocity components in $\mathcal K$ are

$$\begin{align*}
u_\#^\alpha
\to \gamma_\# \begin{pmatrix}
1\\ v^1_\#\\ v^2_\#\\ v^3_\#
\end{pmatrix},\quad
\gamma_\# \equiv \f1{
\sqrt{1 - |v_\#|^2}
}.
\end{align*}$$

And since $\rho_\\#$ is a scalar invariant,
the entire stress-energy-momentum tensor has components

$$\begin{align*}
T_\#^{\alpha\beta}
&\to \rho_\#\ 
\gamma_\# 
\begin{pmatrix}
1 \\ v^1_\# \\ v^2_\# \\ v^3_\#
\end{pmatrix}
\otimes
\gamma_\# 
\begin{pmatrix}
1 \\ v^1_\# \\ v^2_\# \\ v^3_\#
\end{pmatrix}
\\
&=
\rho_\# (\gamma_\#)^2
\begin{pmatrix}
1 & v^1_\# & v^2_\# & v^3_\# \\
v^1_\# & |v^1_\#|^2 & v^2_\# v^1_\# & v^3_\# v^1_\# \\
v^2_\# & v^1_\# v^2_\# & |v^2_\#|^2 & v^3_\# v^2_\# \\
v^3_\# & v^1_\# v^3_\# & v^2_\# v^3_\# & |v^3_\#|^2 \\
\end{pmatrix}.
\end{align*}$$

Pretty busy at first glance, 
but you can see the symmetries if you look carefully!

### 3. Summing together: stress-energy-momentum tensor for fluids

The components of the stress-energy-momentum tensor are all **intensive** quantities,
so it makes physical sense to add them together 
*if* we premultiply each with their particle’s Lorentz-contracted volume in $\mathcal K$,
\\[
\Delta V\_\\# = \f{m\_\\#}{\rho\_\\# \gamma\_\\#},
\\]
thereby making an **extensive** tensor
\\[
T_\\#^{\alpha\beta}
\Delta V\_\\#
\to 
m_\\# \gamma\_\\# 
\begin{pmatrix} 1 & v^1\_\\# & v^2\_\\# & v^3\_\\# \\\ v^1\_\\# & |v^1\_\\#|^2 & v^2\_\\# v^1\_\\# & v^3\_\\# v^1\_\\# \\\ v^2\_\\# & v^1\_\\# v^2\_\\# & |v^2\_\\#|^2 & v^3\_\\# v^2\_\\# \\\ v^3\_\\# & v^1\_\\# v^3\_\\# & v^2\_\\# v^3\_\\# & |v^3\_\\#|^2 \\ \end{pmatrix}.
\\]

After summation, we then divide by the fluid element volume $\Delta V$
to recover the intensive stress-energy-momentum tensor for the fluid:
\\[\begin{align\*}
T^{\alpha\beta} &= \f1{\Delta V} \sum\_\\# T\_\\#^{\alpha\beta} \Delta V\_\\# \\\ &\to 
\f1{\Delta V} \sum\_\\# 
m\_\\# \gamma\_\\# 
\begin{pmatrix} 1 & v^1\_\\# & v^2\_\\# & v^3\_\\# \\\ v^1\_\\# & |v^1\_\\#|^2 & v^2\_\\# v^1\_\\# & v^3\_\\# v^1\_\\# \\\ v^2\_\\# & v^1\_\\# v^2\_\\# & |v^2\_\\#|^2 & v^3\_\\# v^2\_\\# \\\ v^3\_\\# & v^1\_\\# v^3\_\\# & v^2\_\\# v^3\_\\# & |v^3\_\\#|^2 \end{pmatrix}
\\\ &=
\boxed{
\f{N}{\Delta V} \ang{
m\_\\# \gamma\_\\# 
\begin{pmatrix} 1 & v^1\_\\# & v^2\_\\# & v^3\_\\# \\\ v^1\_\\# & |v^1\_\\#|^2 & v^2\_\\# v^1\_\\# & v^3\_\\# v^1\_\\# \\\ v^2\_\\# & v^1\_\\# v^2\_\\# & |v^2\_\\#|^2 & v^3\_\\# v^2\_\\# \\\ v^3\_\\# & v^1\_\\# v^3\_\\# & v^2\_\\# v^3\_\\# & |v^3\_\\#|^2 \end{pmatrix}
}
}.
\end{align\*}\\]
where the angle brackets $\ang{\cdots}$ denote an average over the ensemble of $N$ particles.
This is the formula of $\vec T$ for an arbitrary relativistic fluid composed of dust!

### 4. The case for a perfect fluid

If we make some assumptions about this dust, then we can simplify the expression some more.
Suppose every particle has the same rest mass $m$, then

$$\begin{align*}
T^{\alpha\beta} \to \f{Nm}{\Delta V}
\ang{
\gamma_\# 
\begin{pmatrix}
1 & v^1_\# & v^2_\# & v^3_\# \\
v^1_\# & (v^1_\#)^2 & v^2_\# v^1_\# & v^3_\# v^1_\# \\
v^2_\# & v^1_\# v^2_\# & (v^2_\#)^2 & v^3_\# v^2_\# \\
v^3_\# & v^1_\# v^3_\# & v^2_\# v^3_\# & (v^3_\#)^2 \\
\end{pmatrix}
},
\end{align*}$$

Furthermore, if the particle velocities are isotropic --- that is, uniformly distributed in direction ---
then the average value of every velocity component vanishes:
$\ang{ v^1 } =\ang{ v^2 } =\ang{ v^3 } = 0$,
as do products of independent components.
Their magnitudes are equipartitioned between the three spatial directions:
$\ang{ |v^1\_\\#|^2 } = \ang{ |v^2\_\\#|^2 } = \ang{ |v^3\_\\#|^2 } = \tfrac13\ang{|v\_\\#|^2} $.

$$\begin{align*}
T^{\alpha\beta} \to \f{Nm}{\Delta V}
\ang{
\gamma_\#
\begin{pmatrix}
1 & \\
& \tfrac13|v_\#|^2 & \\
& & \tfrac13|v_\#|^2 \\
& & & \tfrac13|v_\#|^2 \\
\end{pmatrix}
}.
\end{align*}$$

We can identify the 00 component as the energy-momentum density
\\[
\boxed{
\rho = \f{Nm}{\Delta V}\ang{\gamma_\\#}
},
\\]
and the 11, 22, 33 components as the isotropic pressure
\\[
\boxed{
P = \f{Nm}{3\Delta V}\ang{\gamma_\\# |v_\\#|^2}
}.
\\]
The off-diagonal entries are zero: **A perfect fluid has no shear!**
It can be completely described by $\rho$ and $P$.
In tensor form,
\\[
\boxed{
\vec T = (\rho + P)\vec u\otimes\vec u + P\vec g
},
\\]
where $\vec u$ is the four-velocity of the fluid element and $\vec g$ is the metric.
In the fluid element’s rest frame, this combination of $\vec u$ and $\vec g$ produces the components
\\[
T^{\alpha\beta}
= (P + \rho) u^\alpha u^\beta + P g^{\alpha\beta}
\to \begin{pmatrix} \rho\\\ & P\\\ & & P\\\ & & & P \end{pmatrix},
\\]
exactly what we had before.

To make further simplifications, we would need to know the distribution of velocity magnitudes $v_\\#$.
That is a problem for statistical mechanics.
Then, we can determine the fluid’s **equation of state** $w \equiv P/\rho$.


