---
title:  "How do fundamental particles work?"
date: 2024-11-08
---

*Inspired by some exposition from Bea Noether, [Mark Srednicki](https://web.physics.ucsb.edu/~mark/qft.html), and [Quevedo & Schachner](https://arxiv.org/abs/2409.09211).*

## What are particles?

Our modern understanding of particle physics comes from the intersection of two revolutionary theories in the early 20th century:
quantum mechanics and special relativity.

In quantum mechanics, the state of the universe is described by some vector living in an **abstract vector space**. 
This vector space is the set of all possible configurations the universe can be in, *not* our literal, physical spacetime. 
Inside it, the universe’s state vector evolves over time via a special kind of transformation that preserves its total probability in measurements ---
i.e. the probabilities of measuring the universe to be “X” or “not X” always adds up to one. 
This is just like how rotational transformations preserve lengths in ordinary physical space. 
We call probability-conserving transformations of quantum state vectors **unitary transformations**.

In special relativity, every event in the universe is described by 4 coordinates: 
3 for space and 1 for time, constituting a vector in physical spacetime. 
These event vectors can be observed by different people in different places traveling at different velocities --- different **reference frames**. 
Special relativity relates how events are seen in different
reference frames by transformations which preserve the speed of light, which is experimentally measured to always be ~300 million m/s no matter where you are or how fast you’re moving. 
In fact, the constancy of the speed of light is necessary to prevent backwards time travel, i.e. preserve the notion of cause-and-effect. 
We call these causality-conserving transformations of spacetime events **Poincaré transformations**.

These two descriptions of the universe: quantum mechanics and special relativity, are not only extremely mathematically elegant, 
but have by now also withstood a century of experimental tests. 
It then stands to reason that a complete description of fundamental physics should be a unification of the description from both theories: 
States of the universe must transform both *unitarily* in the abstract quantum vector space, and *Poincaré*-ly in physical spacetime. 
Fundamental particles --- the building blocks of the universe --- must then be the irreducible building blocks of these states.
(In math-speak, they are unitary, irreducible **representations** of the Poincaré group of transformations.)

## How to particles interact?

Mathematicall, there are actually many redundancies in this description of fundamental particles: 
a single physical state in the universe (comprised of a collection of particles) can correspond to an infinite set of different representations, 
each differing from each other by a transformation in the abstract vector space. 
These **gauge transformations** do not change any observable phenomenon, 
but are inherent to this mathematical formalism. 
However, weirdly enough, they do have physical meaning!

Firstly, just as rotations conserve lengths, gauge transformations can conserve some quantity, 
called its corresponding **Noether charge**. 
Particles that transform under a gauge are said to be **charged** under that gauge --- 
meaning they *interact* according to the strengths of these Noether charges.
Why is this? Well, the equations which govern the motions of the universe can all be derived from one function (the **Lagrangian**), 
which takes in its state and spits out a single number for each point in spacetime. 
But even when these equations are invariant under gauge transformations, the Lagrangian may change, since the state itself isn’t invariant either. 
To make the Lagrangian itself invariant, we are forced to introduce additional particles (**gauge particles**) whose terms in the Lagrangian 
transform exactly opposite to the terms from other particles that are charged under the same transformation. 
For this exact cancellation to occur, these gauge particles must know about these charged particles’ existence and reflect that in their behaviors. 
In other words, they must interact with their corresponding charged particles.

The exact nature of these interactions depend on the kind of gauge transformation they are derived from. 
For example, the electron’s representation can change by an arbitrary local 1-parameter unitary transformation, $U(1)$, 
without changing the equations of motion in electrodynamics. 
This gives rise to a $U(1)$ charge: the electric charge; 
as well as a gauge particle that transforms like a regular vector in spacetime: the photon.

In this model, charged particles can only directly interact with their partner gauge particles. 
When we see charged particles attracting or repelling each other, they are actually exchanging gauge particles who serve to transmit the interaction. 
We therefore call gauge particles **force carriers**. 
Continuing with the previous example, photons are said to be the carriers of the electromagnetic force.

But the electromagnetic force cannot explain one very important phenomena: 
the existence of atomic nuclei, 
composed of a bunch of tightly-packed protons and neutrons, 
each carrying either zero or positive electric charge. 
Without the existence of another fundamental force ---
equivalently, another gauge particle, Noether charge, or gauge transformation ---
atoms would fly apart. 
Which would be kinda bad.

This is the “strong” force, 
mediated by a gauge particle called the “gluon”, interacting with a “color charge”, 
conserved under 3-parameter, continuous unitary transformations $SU(3)$. 
It is attractive with particles carrying color charge, like gravity is with “mass charge”, 
only much stronger (as its name implies) at the relevant subatomic distances.
Protons and neutrons are actually both net neutral in color charge, 
but they are made of fundamental “quark” particles that are charged under $SU(3)$, 
and the attraction between their quarks binds the nucleus together.

## What are antiparticles?

Alongside these continuous gauge transformations,
there are also discrete transformations that we can imagine performing on the universe.
One is **charge conjugation**, which changes a particle to its antiparticle, and vise versa. 

But what are antiparticles, and how do we know they exist in the first place? 
It turns out without them, our model of particle physics would violate causality, a principle both from special relativity, and from common sense.

This all boils down to the idea of **commuting operations**: 
If two operations $X$ and $Y$ are performed on the universe’s state, 
does it matter if $X$ happens before $Y$, or $Y$ before $X$? 
If it doesn’t, 
then $X$ and $Y$ do not have a fixed ordering in time. 
One does not cause the other. We say $X$ and $Y$ commute, 
or that they are **causally disconnected**. 
Xccording to special relativity, 
then if $X$ and $Y$ are performed in two places sufficiently far apart that they cannot possibily have communicated prior to the operations, 
then they must be causally disconnected.

But suppose our universe contains a particle species “$A$”, as well as two operations “$a^+$” and “$a^-$”.
The $a^+$ operation creates an $A$ particle in the universe, while $a^-$ removes an $A$ from it.
Do $a^+$ and $a^-$ commute, even when performed with arbitrary distance separating them?
Generally not! Suppose our universe began with zero $A$ particles.
If we perform $a^+$, then $a^-$, we would be left with zero $A$ particles still.
But if we perform $a^-$first, there’s nothing to destroy.
Performing $a^-$, then $a^+$ would leave one $A$ particle in the universe.
So the ordering of $a^+$ and $a^-$matters- there is cause-and-effect here.
Somehow $a^-$ "knows" if the other side of the universe contains a particle $A$ or not, violating causality!

The resolution to this dilemma is to redefine how $a^-$ acts on an empty universe.
Instead of needing to know whether there is or is not a nonzero amount of a particle $A$ in the entire universe, 
we say $a^-$creates an antiparticle $A^\*$, 
which travels at or below light speed, and can annihilate a normal particle $A$ if it ever finds one.
What charge must this antiparticle have? 
Because charge is conserved (again, due to the invariance of the universe under gauge transformations), 
the antiparticle-particle pair must have a net charge equal to the zero charge left in the universe after their annihilation. 
Antiparticles must then carry *exactly the opposite charge* as their regular particle partners, 
be it electric or color. (But not mass, for reasons I do not know.)
