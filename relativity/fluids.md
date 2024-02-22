---
title:  "Relativistic fluids from relativistic dust"
---

*[Xing’s notes](/notes)*. This article demonstrates a method to construct the energy--momentum tensor
of a perfect fluid from “dust”, and aims to give an intuition of this tensor's meaning.
It builds upon the lectures of Professor Ori Ganor.

## Background

*(This section is meant as a brief review for those familiar with the basics of special relativity.
For a much more thorough introduction, I highly recommend 
[Eigenchris’s YouTube series on relativity](https://www.youtube.com/watch?v=bEtBncTEc6k&list=PLJHszsWbB6hqlw73QjgZcFh4DrkQLSCQa)
and
[tensor algebra](https://www.youtube.com/watch?v=8ptMTLzV4-I&list=PLJHszsWbB6hrkmmq57lX8BV-o-YIOFsiG),
as well as the first part of
[Sean Carroll’s lecture notes on general relativity](https://arxiv.org/abs/gr-qc/9712019).
Conversely, feel free to skip this section if you want to get right to the derivations.)*

A key postulate of special relativity is

### 1. The speed of light in vacuum is constant in every inertial reference frame.

In other words, $c$ is **invariant**. 
This motivates the definition of invariant **spacetime intervals**:
\\[
(\d s)^2 \equiv -(c\,\d t)^2 + (\d x)^2 + (\d y)^2 + (\d z)^2.
\\]
Just like how the lengths of everyday objects do not change under spatial rotations,
the spacetime interval between events do not change under Lorentz boosts.
Even as strange effects like length--contraction and time--dilation become visible when objects are 
moving relative to each other near the speed of light, the spacetime interval stays constant.
It’s like having a new and improved, *relativistic* ruler.

To produce the invariant spacetime interval with this ruler, it must be able to measure time alongside the spatial lengths.
Geometrically, this motivates adding a zeroth $t$ element to our vectors’ existing three $x$, $y$, and $z$.
For example, a position vector is now

$$ \vec r = r^\alpha \vec e_\alpha. $$

(To simplify the writing, Greek-letter indices like $\alpha, \beta, \mu, \nu$ implicitly range over 0, 1, 2, 3;
Repeated upper--lower index pairs like ${^\beta}{_\beta}$ are implicitly summed over. 
Do not confuse them with exponents, which will always be outside parentheses!
This notation is known as the **Einstein summation convention**.)

This way, a spacetime displacement

$$ \d\vec s = s^\alpha \vec e_\alpha \equiv \d t\,\vec e_0 + \d x\,\vec e_1 + \d y\,\vec e_2 + \d z\,\vec e_3 $$

with components 

$$ s^\alpha \to 
\begin{pmatrix} s^0 \\ s^1 \\ s^2 \\ s^3 \end{pmatrix} \equiv
\begin{pmatrix} \d t \\ \d x \\ \d y \\ \d z \end{pmatrix} $$

has the spacetime interval (again, summing over repeated index pairs)

$$ (\d\vec s)^2 = \eta_{\alpha\beta} s^\alpha s^\beta, $$

where

$$ \vec\eta = \eta_{\alpha\beta} \ (\vec e^\alpha \otimes \vec e^\beta) $$

is the **Minkowskian metric** with components

$$ \eta_{\alpha\beta} \to \begin{pmatrix}
-1 \\ & 1 \\ & & 1 \\ & & & 1
\end{pmatrix}. $$

(When dealing with objects having two indices, the left index corresponds to the row, 
and the right index the column. For instance, because the metric is symmetric, $\eta_{\alpha\beta} = \eta_{\beta\alpha}$.)

The metric is useful for a lot of things (see the links above for more details), but most fundamentally, 
it underscores *an interweaving of space and time* that is thematic throughout
special and general relativity. This means

### 2. To gain full intuition, we must treat space and time on an equal footing.

I strongly believe that to fully understand special relativity, 
it must feel natural to mix quantities derived from space and those derived from time.
For example:
* Velocity is both motion through space and motion through time.
  In your **rest frame** $K$, you are stationary in space, but still traveling in time.
  Your **four--velocity** $\vec u$ thus has rest-frame components

$$ u^\alpha \to \begin{pmatrix}1\\0\\0\\0\end{pmatrix}. $$

* Likewise, momentum is the motion of mass--energy through space and through time.
  A baseball flying through space can pack a punch with its spatial momentum components,
  but sitting on a table, it is still carrying all of its mass--energy in its flight through time.
  In the baseball’s rest frame, its **four--momentum** is thus

$$ p^\alpha = mu^\alpha \to \begin{pmatrix}m\\0\\0\\0\end{pmatrix}. $$

*(We can create new **four--vectors** by multiplying them with scalar invariants.
Here, $m$ always refers to “the mass as measured in the rest frame,” which is by-definition invariant.)*

These intertwined interpretations are key to understanding the stress--energy--momentum tensor. More on that in a bit.

For now, notice that we have only defined these quantities in the rest frames of their respective objects.
With all the length--contraction and time--dilation stuff, how do we calculate how these quantities might look 
when we are moving relative to these objects near the speed of light? Luckily,

### 3. We can transform between any inertial reference frames via Lorentz boosts

Suppose another inertial reference frame $K'$ is traveling with (ordinary) velocity $v\vec e_1$ as measured in $K$.
If we know that the four--velocity $\vec u$ has components $u^\alpha$ in $K$, then it will have the components

$$ u^{\alpha^\prime} = \Lambda{^{\alpha^\prime}}_\alpha u^\alpha $$

in $K'$, where the **Lorentz boost** in $v\vec e_1$ is

$$\begin{align*}
    \Lambda{^{\alpha^\prime}}_\alpha
    \equiv \frac{\d x^{\alpha^\prime}}{\d x^\alpha} 
    \to \begin{pmatrix}
        \gamma_v & -\gamma_v v\\\
        -\gamma_v v & \gamma\\\
        & & 1\\\
        & & & 1
    \end{pmatrix},\quad
    \gamma_v \equiv \frac1{\sqrt{1 - v^2}}.
\end{align*}$$

(This section is unfinished because I have a midterm coming up lol)

### Summary of the notation

* We work in natural units, where the speed of light $c=1$. 
* The metric is spacelike: 
* Greek indices ($\alpha,\beta,\mu,\nu,\cdots$) range 0, 1, 2, 3.
* Contravariant components have upper indices; covariant components have lower indices.
* Implicit Einstein summation convention is assumed for matching upper--lower indices.

With these things in mind, let’s get dusty!

## Derivation

A fluid is a macroscopic phenomenon composed of mass--energy packets traveling at significant speeds relative to neighboring mass--energy packets.
It has emergent properties, such as density and pressure, which are not well-defined for its individual constituents.
Examples include water, neutron star interiors, and photons scattering in a sparse medium (but not, for instance, photons concentrated in a laser beam).

In particular, we are interested in the behavior of perfect fluids, which are isotropic in their rest frame, and have 
negligible particle--particle interactions (except some scattering to keep them isotropic).

## Intuition
