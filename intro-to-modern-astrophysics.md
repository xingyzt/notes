---
layout: default
title: Xing’s notes on Carroll & Ostlie, Introduction to Modern Astrophysics
---

{% raw %}

# Carroll & Ostlie, Introduction to Modern Astrophysics

Xing’s notes

## 1. The Celestial Sphere

### Altitude–Azimuth Coordinate System

* Local to observer location and day of year
* **Altitude** $h$ measured from horizon towards zenith.
* **Zenith** $z = 90^\circ - h$
* **Azimuth** $A$ measured from north on the horizon, clockwise, first towards
  east on the horizon

### Equatorial Coordinate System

* **Celestial equator**: plane passing through Earth’s equator
* **Ecliptic**: path of sun across seasons
* **Vernal equinox**: When ecliptic/sun crosses the celestial equator northward
* **Declination** $\delta$ measured in degrees N/S of celestial equator
* **Right ascension** $\alpha$ measured in hour-minute-seconds from vernal equinox
  eastward to hour circle of object
* **Meridian**: great circle intersecting observer and the poles
* **Hour angle** $H$ angle between object and meridian of observer
* **Local sidereal time**: time elapsed since vernal equinox last traversed
  the meridian; hour angle of that intersection

### Precession

* ~25,770 year cycle of wobbling
* Necessitates use of epochs in the coordinate system
* **J2000.0**: Standard epoch, zero at noon 2000−01−01 UT
  (universal time, measured from Greenwich)
* Uses Julian Calendar (365.25 d/yr), without Gregorian corrections

### Time standards

* **Julian Date** (JD): Days since Noon 4713BC−01−01
* **Heliocentric Julian Date** (HJD): From center of Sun
* **Terrestrial Time** (TT): From surface of Earth factoring in relativity

### Motions

* **Radial velocity** $v_r$ (positive is away)
* **Proper motion** $\mu$ (arcseconds/year)
* Transverse/**tangential velocity** $v_\theta = \mu \cdot r$

### Spherical Trigonometry

* Sides with arclengths $a, b, c$ meet at angles $A, B, C$
* **Law of sines**: ${ \sin a \over \sin A} = {\sin b \over \sin B} = {\sin c \over \sin C}$
* **Law of cosines for sides**: $\cos a = \cos b \cos c + \sin b \sin c \cos A$
* **Law of cosines for angles**: $\cos A = − \cos B \cos C + \sin B \sin C \cos a$
* For a proper motion with position angle $\phi$ (direction measured) and
  arclength $\Delta\theta$,
  * $\Delta\alpha = {\Delta\theta \sin \phi \over \cos \delta}$
  * $\Delta\delta = \Delta\theta \cos \phi$
  * $(\Delta\theta)^2 = (\Delta\alpha \cos \delta)^2 + (\Delta\delta)^2$

## 2. Celestial Mechanics

### Ellipses

* **Semi−major axis** $a$; **Semi−minor axis** $b$
* **Eccentricity** $e : ae = {d\over2}$, where $d$ is distance between foci
* $b^2 = a^2 (1 − e^2)$
* Polar form $r = {a (1 − e^2) \over 1 + e \cos \theta}$
* Area $A = \pi ab$

### Kepler’s Laws

1. Ellpitical orbits: $r = {L^2/\mu^2 \over GM(1 + e \cos \theta)}$
2. Equal area over equal time: ${dA \over dt} = {L \over 2\mu}$ (constant)
3. Period–axis relation: $P^2 = ka^3$, $k = {4\pi^2 \over G(M+m)} = 1 \text{ year}^2\text{ AU}^{-3}$

### Orbital equations

* **Virial theorem**: $\langle E \rangle = {\langle U \rangle \over 2}$, $\langle U \rangle = −2\langle K\rangle$
* **Vis−viva**: $v^2 = G(M+m)\left({2\over r} − {1 \over a}\right)$

## 3. Continuous Spectrum of Light

### Parallax and parsecs (pc)

* Parallax angle $p''$ (arcseconds), subtended by quarter orbit of Earth (1 AU)
* Distance $d_\text{AU} = {206,265 \over p''}$, $d_\text{pc} = {1 \over p″}$

### Apparent magnitude

* Originally compiled by Hipparchus: 1 for brightest, 6 for dimmest
* Modern definition: ${F_2 \over F_1} = 100^{(m_1 − m_2)/5}$
  * Where radiant flux is density of emitted light: $F = {L\over 4\pi r^2}$

### Absolute magnitude

* $M$ for stars, apparent magnitude measured as if it’s located at 10 pc
* $H$ for planets, measured if located at 1 AU
* **Distance modulus** $d = 10^{(m − M + 5)/5}$,
  $m − M = 5 \log_{10}d − 5$
* For equally distant stars, ratio of radiant fluxes equals
  ratio of luminosities, thus ${L_2 \over L_1} = 100^{(M_1 − M_2)/5}$
* $M_\odot = +4.74, L_\odot = 3.839 \times 10^{26} \text{ W}$

### Light

* **Wave speed** $c = \lambda\nu$
* **Poyinting vector** $\vec S = {\vec E \times \vec B \over \mu_0}, \langle S \rangle = {E_0B_0 \over 2\mu_0}$
* **Radiation pressure** $P_\text{absorption} = {\langle S \rangle A \cos \theta \over c}$,
  $P_\text{reflection} = {2\langle S \rangle A \cos^2\theta \over c}$

### Blackbody radiation

* **Wien’s displacement**: $\max(\lambda) = {0.002898 \text{ m K} \over T}$
* **Stefan–Boltzmann**: $L = 4\pi R^2 \sigma T^4$
  * Where $T$ is effective temperature of the star’s surface,
    and $\sigma = 5.670 \times 10^{-8} \text{ W}\;\text{m}^{-2}\;\text{K}^{-4}$
* **Planck’s function**: $B(T) = {2hc^2 / \lambda^5 \over e^{hc / \lambda kT} − 1}$

### Color index

* **Bolumetric magnitudes**: measured over all wavelengths
* In practice, detectors only measure certain ranges, to varying
  degrees of sensitivity
* Standard filters:

|     | name   | center | bandwidth |
| --- | ------ | ------ | --------- |
| U   | UV     | 365 nm | 68 nm     |
| B   | blue   | 440 nm | 98 nm     |
| V   | visual | 550 nm | 89 nm     |

* Isolated letters are apparent magnitudes; $M_x$ are absolute
* U−B color index: $U − B = M_U − M_B$
* **B−V color index**: $B − V = M_B − M_V$ (smaller is bluer)
* **Bolumetric correction**: $BC = m_\text{bol} − V = M_\text{bol} − M_V$
* **Color−color diagrams** graph U−B against B−V,
  showing stars are non−ideal blackbodies

## 4. Special Relativity

* _Notation: prime (′) means observed, unprimed means rest_

### Einstein’s postulates:

1. Principle of Relativity: Laws of physics are the same in all inertial
   reference frames
2. Constancy of the Speed of Light: In vacuum, $c$ is independent of the motion
   of the light source

### Lorentz transformations:

* **Lorentz factor** $\gamma \equiv {1 \over \sqrt{1 − {u^2 \over c^2}}}$
* **Time dilation** $t' = \gamma t$
* **Length contraction** $x' = {x\over\gamma}$
* **Proper time** $\tau$ and **proper length** measured at rest w.r.t. events ($\gamma = 1$)

### Relativistic doppler shifts

* **Doppler shift** $\nu' = {\nu \over \gamma (1 + {u \cos \theta \over c}) }$
  ($\theta = 0^\circ$ away from observer, $180^\circ$ towards observer)
* **Redshift** $z = {\lambda'\over\lambda} − 1 = \sqrt{ 1 + {u \cos \theta \over c} \over 1 − {u \cos \theta \over c} }$
* Time dilation from redshift ${t' \over t} = z+1$

### Relativistic momentum and energy

* **Rest mass** $m_0$; **Rest energy** $E_0 = m_0c^2$
* **Relativistic momentum** $p = \gamma m_0v$
* **Relativistic kinetic energy** $K = m_0c^2(\gamma − 1)$
* **Relativistic total energy** $E = \gamma m_0 c^2 = \sqrt{ p^2c^2 + m_0^2c^4 }$

## 5. Interaction of Light and Matter

* **Fraunhofer lines**: absorption lines in the Solar spectrum

### Kirchoff’s Laws

1. Hot, dense gas/solid  produces continuous, lineless spectrum
2. Hot, diffuse gas produces emission lines
3. Cold, diffuse gas in front of a continuous spectrum source produces
   absorption lines in the spectrum

### Spectrographs

* Produced via diffraction grating: $d \sin \theta = n\lambda$,
  $d$ is grating separation; $n$ is spectral order
* Resolving power ${\lambda\over\Delta\lambda} = nN$,
  $\Delta\lambda$ is smallest resolvable wavelength; $N$ is # of gratings illuminated
* Low−speed Doppler shift approximation: ${\Delta\lambda\over\lambda} = {u \cos \theta\over c}$

### Photoelectric effect

* Shining light onto metal ejects electron with energy proportional to frequency, not intensity
* Max kinetic energy of an ejected electron $K = h\nu − \phi$
* Work function $\phi$ is the metal’s min electron binding energy
* Demonstrates quantization of light into packets of energy $E = h\nu$

### Compton scattering

* Relativistic scattering of electrons by high−energy photons
* Photon loses wavelength $\Delta\lambda = {h(1 − \cos \theta)\over m_ec}$, $m_e = 9.109^{-31} \text{ kg}$

### Bohr’s semiclassical atom

* Assumes quantization of angular momentum $L = \mu v r = n \hbar$,
  $n$ is principle quantum number
* However, still puts electrons in classical circular orbits modeled by
  Coulombic attraction
* Kinetic energy of atom $K = {\mu v^2 \over 2} = {e^2 \over 8 \pi \epsilon_0 r} = {n \hbar^2 \over 2\mu r^2}$
* Allowed orbital radii $r = a_0n^2$,
  **Bohr radius** $a_0 = {4\pi \epsilon_0 \hbar^2 \over \mu e^2} \approx 5.292^{-11} \text{ m}$
* Allowed energy levels: $E = -{E_0 \over n^2}$, $E_0 \approx 13.6 \text{ eV}$

### De Broglie waves

* Generalizes photon wavelength-frequency-momentum relation to all partilces
* $\nu = {E \over h}, \lambda = {h \over p}$

### Heisenberg’s Uncertainty Principle

* $\Delta x \Delta p \geq {\hbar\over2}$, estimated as $\Delta x \Delta p \approx \hbar$ , equivalent to $\Delta E \Delta t \approx \hbar$
* **Reduced Planck constant**: $\hbar = {h \over 2\pi} = 1.054 \times 10^{-34} \text{ J}\cdot\text{s} = 6.582 \times 10^{-16} \text{ eV}\cdot\text{s}$
* Allows for quantum tunnelling when barrier is within a few wavelengths wide
  * e.x. total internal reflection can be disrupted by placing a another prism sufficiently close to the boundary

### Schrödinger's quantum atom

* Schrödinger equation models particles as probability waves
* Analytically solvable for H with same result as the Bohr atom
* Quantum orbitals are probability density clouds
* Adds additional quantum numbers $l \in [0, n-1]$ and $m_l \in [−l, +l]$
* Allowed angular momentum magnitudes $L = \hbar \sqrt{ l (l + 1) }$
* Allowed angular momentum $z$-components $L_z = \hbar m_l$
* Different orbitals with same $n$ are degenerate
  (i.e. have the same energy) *unless* under a magnetic field

### Normal Zeeman effect

* Previously degenerate orbitals split into slightly different energies
* 3 frequencies: $\nu \in \{ \nu_0, v_0 \pm {eB \over 4\pi\mu} \}$
* Allows detection of magnetic field strength around e.g. starspots
* Even if field is too weak to cleanly split,
  different polarizations are still detectable

### Anomalous Zeeman effect and Spin

* Anomalous Zeeman effect: more complex splitting patterns
* Resulted in discovery of the spin quantum number
* Allowed spin angular momentum magnitudes $S = \hbar{\sqrt3\over2}$
* Allowed spin angular momentum z-component $S_z = \hbar m_s$

### Fermions and bosons

* Dirac's relativistic Schrödinger's equation divided particles
  into fermions and bosons
  * Also predicted antiparticles: opposite electric changes and magnetic moments
* Fermions have spins of $\text{odd} \times {\hbar\over2}$; Bosons have spins of $\text{integer} \times \hbar$
* Pauli exclusion principle: No two fermions can share the same set
  of 4 quantum numbers $\{n, l, m_l, m_s\}$

### Complex spectra

* Quantum state transitions go from one set of quantum numbers to another
* Selection rules significantly lower the probabilities of certain transitions
* Restrictions: $\Delta l = \pm1, \Delta m_l \in \{ 0, \pm1 \}, m_l$ of orbitals cannot both be 0
* Allowed transitions triggered by collisions or occur spontaneously $\approx 10^{-8} \text{ s}$
* Some **forbidden transitions** occur in very low gas densities
  * e.x. diffuse interstellar medium; outer stellar atmospheres

## 6. Telescopes

### Basic optics

* **Index of refraction** $n(\lambda) \equiv {c \over \nu}$
* **Snell's law**: $n_1 \sin \theta_1 = n_2 \sin \theta_2$
* Lensmaker's formula: $P = {1\over f} = (n - 1)\left( {1\over R_1} + {1\over R_2} \right)
  $, 
  $P$ is optical power, $f$ is focal length, $R_n$ are radii of curvature
* Focal length is wavelength-independent for mirrors
* A lens projects angular separation onto a plate (focal plane)
* Linear separation of point sources on focal plane increases with focal length:
  ${d\theta \over dy} = {1\over f}$

### Resolution

* Single-slit diffraction of wavefronts through aperture limits resolving power
* **Airy disk**: bright spot around central maximum
* Destructive interference creates minimums: $\sin \theta = {m\lambda\over D}, m = 1, 2, 3, ...$
* Images are “unresolved” when one's central maximum overlaps with
  another's 1st minimum
* **Rayleigh criterion** defines arbitrary min resolvable distance
  of a circular aperture $\theta \equiv {1.22\lambda \over D}$, $D$ is aperture diameter
  * Possible to differentiate sources within Rayleigh criterion through careful
    analysis of diffraction patterns
* **Seeing**: ctual ground-based optical resolution worse than ideal due to
  refraction from atmospheric turbulence, unless corected by **adaptive optics**
* Only lenses suffer from **chromatic aberration**;
  somewhat reduceable with correcting lenses
* Spherical lenses and mirrors suffer from **spherical aberration**;
  mitigated by paraboloids, which are harder to produce
* Paraboloid lenses and mirrors suffer from **coma**: elongation of off-axis images
* **Astigmatism**: different parts of a lens or mirror converge an image
  at different locations; correction may result in curvature of field issues

### Brightness

* Illumination $J \propto {1 \over F^2}$, the amount of light energy per second focused onto a unit
  area of the resolved image
* Focal ratio: $F ≡ {f \over D}$
  * Larger aperture increases resolution and illumination
  * Longer focal length increases image size but decreases illumination
  * For fixed focal ratio, larger telescope diameter increases resolution but
    not illumination

### Optical telescopes

* **Refracting telescope**:
  * Objective lens focus as much light as possible onto the focal plane
  * Sensor is placed on focal plane or eyepiece is placed at its focal length away
  * **Angular magnification** $m = {f_\text{objective} \over f_\text{eye}}$
  * Issues:
    * Lens can only be supported near the edge, so gravity deforms heavy ones
  * Entire volume must be fabricated precisely
  * Slow thermal response and thermal expansion
  * Chromatic aberrations
* **Reflecting telescope**:
  * Replaces refractor’s objective lens with a mirror
  * Minimizes weight, support, and distortion issues,
    esp. with active support system
  * Main issue: Prime focus of the mirror is in the path of the incoming light
    * **Newtonian telescope**: flat diagonal secondary mirror;
      suffers from faraway eyepiece that introduces torque
    * **Cassegrain telescope**: parabolic (or hyperbolic for **Ritchey-Chrétien telescope**) primary mirror;
      secondary mirror (usually convex to increase focal length) reflects back
      through a hole in the primary’s center
    * **Coudé telescope**: mirror system directs light to an instruments room;
      very long focal length
    * **Schmidt telescope**: spheroidal primary mirror minimizing coma;
      correcting lens to remove spherical aberrations;
      provides wide field of view (degrees vs arcminutes) with low distortion
* Mounts:
  * Equatorial: polar axis, easy to adjust, hard to build for massive telecopes
  * Altitude-azimuth: easy to build, generally needs computer adjustments
* Adaptive optics: small deformable mirror with many piezoelectric actuators
  counteract atmospheric distortions, constantly calibrated via a guide star
* **Charge-coupled device** (CCD) detects ~100% of incident photons with a linear
  response across wide wavelength and intensity (dynamic) ranges

### Radio telescopes

* Strength measured is spectral flux density $S(\nu)$
* Typical radio source has $S \approx 1 \text{ jansky (Jy)} = 10^{-26} \text{ W}\text{ m}^{-2}\text{ Hz}^{-1}$;
  weak sources ~mJy
* Integrate $S$ over area and bandwidth for total power received
* Long wavelengths require large apertures, but less manufacturing precision
* Addition of other telescopes enables interferometry;
  reduces side lobes and narrows main lobe of antenna sensitivity pattern
* Interferometry determines angle of source from phase difference between antennae
* Pointing angle $\theta | \sin \theta = {L\over d}$,
  $L$ is difference in wavefront’s distances to the antennae;
  $d$ is baseline distance between antennae
* Very long baseline interferometry can span multiple continents

### IR/UV/X-ray/gamma-ray astronomy

* Water vapor is primary contributor to IR absorption; thus low-humidity
  mountain peaks, balloons and aircraft observatories are used
* IR observation also requires very cold telescopes and detectors
* Glass is opaque to UV; UV telescopes need very precise reflecting surfaces
* X-ray and gamma-ray photons penetrate traditional mirrors;
  instead imaged by **graze-incidence** (near 90° incidence) reflections, or by
  **Bragg scattering** through crystal lattice inteference

# Part II. The Nature of Stars

## 7. Binary Systems and Stellar Parameters

* In many cases, allows the calculation of stellar masses, highlighting a
  well-defined mass-luminosity relation for the large majority of stars

### Classifications

* **Optical double**: Fake binaries that just look close
* (The rest of the classes are not mutually exclusive)
* **Visual binary**: Both stars are resolved
* **Astrometric binary**: Visual oscillatory motion of one star implies companion
* **Eclipsing binary: Periodic variation in light curves reveal two stars
* **Spectrum binary**: Different spectral classes / Doppler effect reveal superimposed spectra
* **Spectroscopic binary**: Doppler effect in spectra reveal oscillation of
  radial velocity curves

### Mass determination of visual binaries

* Mass ratio from ratio of subtended angles: ${m_1 \over m_2} = {r_2 \over r_1} = {\alpha_2 \over \alpha_1}$
* If distance to system and inclination are known, masses can be determined:
  $m_1 + m_2 = {4\pi^2 \over G} {d^3 \over \cos^3 i} {\alpha^3 \over P ^2}$
  * Where i is the angle of inclination between the orbit and the sky’s planes

### Mass determination of spectroscopic binarie

* Mainly detected from binaries with high inclinations
* “**Double-line**” if both spectra are visible; “**single-line**” otherwise
* Radial velocity curves cross at the radial velocity of the center of mass
* Radial velocity curve amplitudes scaled uniformly by $\sin i$
* Eccentricity skews the curves, but close binaries tend to quickly
  circularize from tidal forces
* Mass ratio: ${m_1 \over m_2} = {v_2 \over v_1}$
  * Where v are radial velocities
* If inclination is known, masses can be determined:
  $m_1 + m_2 = {P \over 2 \pi G} \left(v_1 + v_2 \over \sin i\right)^3$,
* If only one radial velocity is known (single-line binary),
  then mass ratio cannot be found,
  but the mass function still provides rough constraints:
  ${(m_2 \sin i)^3\over(m_1 + m_2)^2} = {P \over 2\pi G} v_1^3$

### Temperature ratio and radii determination of eclipsing binaries

* Two dips in the light curve:
  * Primary eclipse: hotter star passes behind cooler one
  * Secondary eclipse: cooler star passes behind warmer one
  * Size doesn’t matter
* If smaller star $T_1$ is hotter, temperature ratio from brightness is
  ${B_0 − B_1 \over B_0 − B_2} = \left( T_1 \over T_2 \right)^4$
* Unless extremely close binary, inclination $i \approx 90^\circ$,
  thus allowing accurate mass and velocity determination
* Radius of star $r = {v\Delta t\over2}$
  * Where $v$ is relative velocity,
    and $\Delta t$ is transit time for the smaller star, or eclipse time for the larger

### Computer modeling

* May incorporate tidal deformations, surface tempreature variations,
  flux distributions, etc
* Generates synthetic light curves to compare against observational data

### Extrasolar planets (exoplanets)

* Direct observation very difficult due to vast luminosity differences
* Indirect methods: spectral radial velocities, astrometric wobbles, and eclipses
* First exoplanet around a Sun-like star discovered in 1995
* Many rapidly followed due to dramatic advances in detector technology,
  large-aperture telescopes, and long-term observing campaigns

## 8. Classification of Stellar Spectra

* Stars emit a continuous spectra (the continuum)
  with absorption lines in certain wavelengths

### Harvard spectral types

* Pickering and Fleming initially labeled alphabetically by strength of H absorption lines
* Different strengths due to different temperatures causing different electron
  ionization levels and orbitals
  * $\text{H I}$ (neutral)’s visible spectral (Balmer) lines strongest at A0
    (effective T = 9520 K)
  * $\text{He I}$ (neutral)’s visible spectral lines strongest at B2
    (effective T = 22,000 K)
  * $\text{Ca II}$ (singly ionized)’s visible spectral lines strongest at K0
    (effective T = 5250 K)
* Cannon rearranged by temperature (OBAFGKM, O hottest)
  and added decimal subdivisions: A0 hotter than A9
* Additional spectral types for very cool stars and brown dwarfs: **OBAFGKMLT**

### Spectral physics

* Statistical mechanics studies macroscopic behavior of stellar atmospheres
* **Maxwell−Boltzmann velocity distribution** $\rho(v) = \left(m \over 2\pi k T \right)^{3/2} e^{−mv^2/2kT} 4\pi v^2$
  * Restricted to gasses in thermodynamic equilibrium with $\rho \lesssim 1\text{ kg}\text{ m}^{-3}$
  * Most probable speed $v_\text{mode} = \sqrt{2kT\over m}$
  * Root-mean-square speed $v_\text{RMS} = \sqrt{3kT \over m}$ due to right-skewed exponential tail
* Electrons’ orbital energies are affected by their atoms’ kinetic energies
  through collisions: higher orbitals are less likely to be occupied
* Probability ratio of electron states: ${P(s_1) \over P(s_2)} = e^{−(E_1−E_2) / kT}$
  * As thermal energy $kT \to 0$, $P(\text{ higher energy state }) \to 0$; confined to lowest state
  * As $kT \to \infty$, ${P(s_1) \over P(s_2)} \to 1$; all states are equally likely
* Statistical weights of some energy levels increased by degenerate orbitals
* **Boltzmann equation**, the probability ratio of electron energies:
  ${P(E_1) \over P(E_2)} = {g_1 \over g_2} e^{−(E_1−E_2) / kT}$
  * Where $g$ is the statistical weight, or number of states with that energy
* Balmer lines require exitations from the first excited state $N_2$
* Greater fraction of H I are in the excited state $N_2$ at higher temperatures;
  However, at the same time, significant fractions of $\text{H I}$ are ionized to $\text{H II}$
* **Partition function**, the number of possible electron states weighted by energy,
  varies by ionization state of the atom, thus affecting the probability of
  atoms at different ionization states: $Z = \sum_j g_j e^{−(E_j−E_1) / kT}$
  * For H I, the ground state (with two orbitals $s=\pm{1\over2}$) dominates for most temperatures: $Z \approx 2$
  * For H II, no electrons means only one possible configuration: $Z = 1$
* **Saha equation**, the probability ratio of ionization states:
  ${N_+ \over N} = {2Z_+ \over n_e Z} \left( 2\pi m_e k T \over h^2 \right)^{3/2} e^{−\chi / kT}$
  * Very sensitive to ionization energy $\chi$, as $kT$ only ranges around 0.5−2 eV
  * Also effected by free electron number density $n$, increased by presence of ionized He II&III
* Combining the Boltzmann and Saha equations shows varying narrow partial ionization
  zones for different atoms and ionization states
  * For $\text{H}$, there is significant fraction of electrons in the first excited state for
    temperatures between 8300−11300 K, matching the peaking of Balmer lines ≅ 9900 K

### Hertzsprung–Russell diagram

* Hertzsprung discovered type G and cooler stars had a range of magnitudes;
  termed the brighter ones giants; Russell termed the dimmer ones dwarfs
* H–R diagram plots absolute magnitude (brighter upward) against
  spectral type (warmer leftward)
* 80–90% of stars lie as dwarfs in a main sequence from upper left to lower right
  * Width of main sequence due to varying ages and compositions
* Initial theory was stars evolved from upper left to lower right;
  terminology of “earlier” (warmer) and “later” (cooler) stars remain
* Simple luminosity–temperature relation reveals fundamental dependence on mass
* Iso-radius lines run diagonally roughly parallel to main sequence
  * Warmer, more massive stars have a lower average density
* Supergiants such as Betelgeuse occupy extreme upper right

### Morgan−Keenan luminosity classes

* Maury noted subtle line width variations amongst stars with similar effective
  temperatures and different luminosities; found by Hertzsprung to differentiate
  main-sequence stars and giants
* Morgan and Keenan published atlas appending Roman numeral luminosity classes
  to Harvard spectral types
  * Ia/Ib for supergiants
  * V for main requence
  * VI for metal-deficient subdwarfs
  * Excludes white dwarfs, class D
* Luminosity classes roughly correlates with absolute magnitudes:
  enables placement of star on H−R diagram entirely from its spectrum
* Spectroscopic parallax: calculating distance modulus from the
  spectrally determined absolute magnitude

## 9. Stellar Atmospheres

* Line blanketing: dense series of metallic absorption lines significantly
  decrease intensity relative to ideal blackbody at certain wavelengths

### Radiation field

* Radiation energy originates from thermonuclear reactions,
  gravitational contraction, and core cooling
* But observed light comes from gasses in atmosphere, not the opaque interior
* Temperature, density, and composition of atmospheric layers determine
  spectral features
* Specific ray intensity:
  $I(\lambda) \;d\lambda \equiv {\partial I\over\partial\lambda} = {{\partial E \over \partial\lambda} d\lambda  \over d\lambda \;dt\;dA \cos \theta \;d\Omega }$,
  power transmitted at a certain wavelength within a certain solid angle
* Mean intensity $\langle I(\lambda) \rangle = I(\lambda)$ for isotropic field; $= B(\lambda)$ for blackbody
* Specific energy density $u(\lambda) \;d\lambda = {4\pi\over c} \langle I(\lambda) \rangle \;d\lambda$
  * Total energy density for blackbody $u = {4\sigma T^4 \over c}$
  * Blackbody radiation pressure $P = {u\over3}$
* Specific radiative flux
  $F(\lambda) \;d\lambda = \int I(\lambda) \;d\lambda \cos \theta \;d\Omega$,
  net energy at a certain wavelength flowing towards one single direction
* For a resolved source, specific ray intensity is measured,
  thus detector distance does not affect measured intensity,
  though the angular size decreases
* For an unresolved source, specific radiative flux is measured,
  falling off with ${1\over r^2}$

### Temperature

* **Mean free path** $l = {1\over n\sigma}$,
  * Average distance traveled by particles and photonsbetween collisions
  * Where $n$ is atomic number density and $\sigma \equiv \pi(2a_0)^2$ is collision cross section
* Multiple descriptions of temperature:
  * **Effective temperature**: obtained from the Stefan–Boltzmann law
  * **Excitation temperature**: defined by the Boltzmann equation
  * **Ionization temperature**: defined by the Saha equation
  * **Kinetic temperature**: defined by the Maxwell–Boltzmann distribution
  * **Color temperature**: defined by fitting the Planck function on the spectrum
* Effective temperature is a global descriptor;
  rest apply locally, varying by gas location and other conditions
* All local temperature definitions result in the same value
  at thermodynamic equilibrium
* **Local thermodynamic equilibrium** (LTE) occurs when temperature change is
  very gradual compared to the mean free path length: $H \ll L$
  * Temperature scale height $H(T) \equiv {T \over \left|dT \over dr\right|}$

### Opacity

* Both scattering and pure absorption reduce the intensity of directed light:
  $\;dI = −\kappa\rho I \;ds$
  * Where absorption coefficient (opacity) $\kappa$ is dependent on
    ray wavelength and gas composition, density, and temperature
* Characteristic distance $l = {1\over k\rho}$ is equivalent to main free path for photons
* Pure absorption decreases intensity exponentially: $I = I_0 e^{−\tau}$
  * Where $\tau$ is **optical depth** with $d\tau = −\kappa\rho \;ds$
* Optical depth is equivalent to the number of mean free path lengths along the
  path of a ray; gas is “**optically thick**” if $\tau \gg 1$, “**optically thin**” if $\tau \ll 1$
* **Balmer jump**: opacity of stellar material suddenly increases for wavelengths
  below the ionization energy of H I in the first excited state
* Sources of opacity classified by initial and final quantum states of the
  interacting electrons:
  * **Bound–bound transitions** (excitation and de-excitation)
    * Small except for wavelengths corresponding to specific excitation energies;
    * If de-excitation returns to the initial state through emission,
      the photon is essentially scattered
    * If multiple photons are emitted from one absorption, the average photon
      energy is reduced
    * If de-excitation is instead triggered by collision, then the energy is
      converted into kinetic/thermal
  * **Bound-free absorption** (photoionization) and **free–bound emission** (recombination)
    * Wavelength-dependent cross section comparable to that of collisions
    * Adds to continuum opacity, as any photon above an electron’s
      ionization energy may be absorbed
    * Recombination emits one or more photons,
      again reducing average photon energy and scattering
  * **Free–free absorption** and **bremsstrahlung** (“braking radiation”) emission
    * Free electron gains speed from absorption, or loses speed from emission
    * Must take place next to an ion to conserve both energy and momentum
    * Adds to continuum opacity
  * **Free electron** (Thomson) scattering
    * Photon is scattered when the electron is made to oscillate in its EM field
    * Very tiny wavelength-independent cross section
    * Only dominates in stellar interiors and the hottest atmospheres,
      where near-total ionization eliminates bound-electron processes
  * **Loosely-bound electron** (Compton/Rayleigh) scattering
    * Compton if $\lambda \ll a_0$;
      very small change in scattered photon wavelength, much like Thomson
    * Rayleigh if $\lambda \gg a_0$;
      proportional to ${1 \over \lambda^4}$,
      so only significant in UV for supergiant stars’ extended envelopes,
      cool main-sequence stars, and planet atmospheres
* Electron scattering and photoionization of $\text{He}$ are primary sources of continuum opacity
  for type O
* Photoionization of $\text{H}$ and free–free absorption primary sources for types B–A
* Photoionization of $\text{H}^-$ is primary source for type F0 and cooler
  * Very low binding energy of 0.754 eV ($\lambda = 1640 \text{ nm}$) for bound–free opacity
  * Also contributes to free–free absorption at longer wavelengths
* Molecules survive in planetary and cooler stellar atmospheres
  * High opacity from large number of discrete absorption lines
  * May break apart during absorption through photodissociation
* Rosseland mean opacity: weighted harmonic mean of opacity across all wavelengths,
  accounting for rate of change of the blackbody spectrum with temperature;
  ${1\over\kappa} \equiv {\int {1\over\kappa_\nu} {\partial B_\nu \over \partial T} \;d\nu \over \int {\partial B_\nu \over \partial T} \;d\nu }$
  * Obeys **Kramer’s opacity law**: $\kappa \propto {\rho\over T^3.5}$
  * No analytic solution, but approximations exist for bound–free and free–free opacities
  * $\kappa_\text{bf} \propto {g_\text{bf} \over t} {Z (1+X) \rho \over T^3.5}, k_\text{ff} \propto g_\text{ff} {(1−Z)(1+X) \rho \over T^3.5}$
    * $X$, $Y$, and $Z$ are **mass fractions** of $\text{H}$, $\text{He}$, and metals respectively
    * $g_\text{bf}$ and $g_\text{ff}$ are **Gaunt factors** correcting for quantum effects; ~1 for visible and UV
    * t is the **guillotine factor**, typically ranging 1–100

### Radiative transfer

* No net energy change occurs in any layer of a star that is in steady-state equilibrium
* Emission processes complement each of the primary absorption processes,
  resulting in randomly-directed scattering
* Specific intensity: $dI = −\kappa\rho I \;ds + j\rho \;ds$
  * Where $\kappa$ is the absorption and $j$ is the emission coefficient; both are $\lambda$-dependent
* Equation of radiative transfer: $−{1\over\kappa\rho} {dI\over ds} = {dI\over d\tau} = I − S$
  * Where source function: $S \equiv {j\over\kappa}$
  * Expresses how light beam composition tends to resemble the local source of photons
* In local thermodynamic equilibrium, $S = B$, the blackbody Planck function
  * Integrating over all wavelengths, $S = B = {\sigma T^4\over\pi}$
  * If $\tau \gg 1, I = B$ as well
* Radiation differential “driving” net radiative flux of photons flowing to the surface:
  ${dP\over d\tau} = −{F\over c}$
* Radiative flux throughout star is constant: $F = \sigma T_\text{eff}^4$
* In a random walk with $N$ steps of average size $I$, $d = l\sqrt N$
* On average, a photon from optical depth $\tau$ needs $\tau^2$ steps to reach the surface
  * Photons from $\tau \approx 1$ may escape without scattering
* **Eddington approximation**: $T^4 \approx {3\over4} T_\text{eff}^4 \left(\tau + {2\over3}\right)$
  * Star has its effective temperature at optical depth $\tau \approx {2\over3}$
  * Average observed photon is emitted from $\tau \approx {2\over3}$, independent of angle
* Higher opacity corresponds to shorter pathlength for the same optical depth,
  thus absorption lines must come from outer, cooler layers
* Limb darkening: line of sight reaches $\tau = {2\over3}$ in cooler, dimmer layers when
  observing closer to the edge of the disk

### Spectral line profile

* Core of line is formed at higher and cooler layer;
  formation descends down to continuum region at wings of the line
* Line is **spectrally thin** if radiant flux $F(\lambda)$ is never 0
* **Equivalent width** $W \equiv \int 1 − {F(\lambda) \over F_0} d\lambda$
  * Where the integrand is the depth of the line at $\lambda$
* Broadening processes:
  * **Natural broadening**: uncertainty principle means momentary occupancy of
    an excited state for small $\Delta t$ amplifies uncertainty of orbital energy
    $\Delta E \approx {\hbar \over \Delta t}$, resulting in uncertainty of wavelength
    $\Delta\lambda \approx {\lambda^2 \over 2\pi c} \left( {1 \over \Delta t_0} + {1 \over \Delta t_1} \right)$, where $\Delta t$ are lifetimes in the two states
  * **Doppler broadening**: random thermodynamic motion results in nonrelativistic Doppler shift
    $\Delta\lambda \approx {2\lambda \over c} \sqrt{2kT \over m}$
    * in giant and supergiant stars, random large-scale turbulence increases to
      $\Delta\lambda \approx {2\lambda \over c} \sqrt{2kT\over m + v^2}$, where $v$ is the most probable turbulence speed
    * coherent mass motions such as rotation, pulsation and mass loss also
      substantial factors
  * **Pressure broadening**: collisions with neutral atoms and pressure from nearby
    ion electrical fields may perturb the orbitals with
    $\Delta\lambda = {\lambda^2 \over \pi c \Delta t_0}$, where $\Delta t_0 \approx {1\over n\sigma\sqrt{2kT \over m} }$ is the average collision time
    * Dependence on atomic number density explains narrower lines of sparser
      giant and supergiant extended atmospheres used for the
      Morgan–Keenan classification
* Damping/**Lorentz profile**: shape of lines from natural and pressure broadening,
  characteristic of radiation from electric charge in damped harmonic motion
* **Voigt profile**: total line profile from both Doppler and damping profiles,
  with Doppler dominating at the core and damping dominating at wings
* **Schuster–Schwarzchild model**: assumes photosphere is a blackbody
  and atoms above it create absorption lines
* Number of atoms involved in absorption per surface area: $fN$
  * Where oscillator strength $f$ is the relative likelihood of each transition,
    and column density $N$ is the area density of absorbing atoms in a column from the surface to the observer
* **Curve of growth** for line width as a function of column density:
  * Initially optically thin core: $W \propto N$
  * Then, saturated core with optically thin wings: $W \propto \sqrt{\ln N}$
  * For high $N$, pressure-broadening profile dominates wings: $W \propto \sqrt N$
* Applying Boltzmann and Saha equations to curve of growth finds total number of atoms above continuum layer
  * Applying Boltzmann equation also finds excitation temperature
  * Applying Saha equation also finds either $\text{e}^-$ pressure or ionization temperature from the other

### Computer modeling

* Construction of an atmospheric model with extensive physics fine-tuned to observations
  can provide information on line profile, chemical composition, effective temperature,
  surface gravity, etc

## 10. The Interiors of Stars

* Direct observation only possible with neutrinos and the occasional supernovae
* Study requires physically accurate computer models that match observable surface features

### Hydrostatic equilibrium

* Equilibrium condition requires pressure gradient to counteract gravity: ${dP \over dr} = −\rho g$
  * **Local gravity** $g \equiv {G M_r \rho \over r^2}$, where $M_r$ is mass enclosed within $r$
* Mass conservation ${dM_r \over dr} = 4\pi r^2 \rho$
* Total pressure $P = {\rho kT \over \mu} + {aT^4\over3}$
  * Where first term is ideal gas law and second term is radiation pressure
  * With as mean molecular weight

### Kelvin–Helmholtz mechanism

* Contraction releases gravitational energy $\Delta E \approx {3\over10} {GM^2 \over R}$
* Gravitational energy contribution $\epsilon = −{dQ \over dt} = −T {dS \over dt}$
  * Where specific entropy: $dS \equiv {dQ \over T}$
  * Contraction produces heat and decreases entropy;
    vise versa for expansion
* **Kelvin–Helmholtz timescale** of star $t = {\Delta E\over L}$, $t_\odot \approx 10^7 \text{ years}$

### Nuclear fusion

* Releases difference in strong nuclear force binding energy:
  * $E = \Delta m_0c^2 = ( Z m_\text{p} + (A−Z) m_\text{n} − m_\text{nucleus} ) c^2$
* Timescale $t_\odot \approx 10^{10}$
* Must overcome Coulombic repulsion between positive nucleons: ~10⁶ eV at $r \approx 10^{-15} \text{ m}$
* Thermal energy of gas: $E = {\mu v^2 \over 2} = {p^2 \over 2\mu}$
  * Where $\mu$ is reduced mass and $v$ is average relative velocity
* Requires impossibly high temperatures in classical physics:
  $E = {3\over2} kT = {q_1q_2 \over 4\pi r \epsilon_0} \to T \approx 10^{10} \text{ K}$ for \text{H-H}$ fusion
* With quantum tunneling (de Broglie $\lambda = {h \over p} \approx r$):
  $E = {\left(h\over\lambda\right)^2\over2\mu} = {q_1q_2 \over 4\pi\lambda\epsilon_0} \to T \approx 10^7 \text{ K}$, consistent with solar core temperature estimates
* Reaction rate depends on velocity (wavelength) distribution, densities, and
  cross-sectional area of particles 1 and 2: $r = \int n_1 n_2 \sigma(E) v(E) { n(E)\over n} dE$
  * Power-law approximation: $r = r_0 X_1 X_2 \rho^{\alpha'} T^\beta$,
    * Where $X$ are mass fractions, $\alpha' \approx 2$, and $\beta$ ranging from ~1 to ≥40
* Velocity distribution follows Maxwell–Boltzmann
* Cross section $\sigma(E) = {S \over E} e^{−b/\sqrt E}$
  * Where de Broglie area is $\pi\lambda^2 \propto {1\over p^2} \propto {1\over E}$,
    and tunneling across barrier of energy $U$ is $e^{‐2\pi^2 U/E} \propto e^{−b/\sqrt E}$,
    with $b \propto q_1q_2\sqrt\mu$
  * $S$ may be slow-varying function of $E$, or may have sharp peaks from resonance of
    specific energy levels within the nucleus
* **Electron screening** from sea of free $\text{e}^-$ partially hides nuclei,
  reducing effective charge and Coulomb barrier,
  sometimes enhancing $\text{He}$ production by 10–50%
* Likelihood of nuclear reaction: $e^{−b/\sqrt E} e^{−E/kT}$
  * Product of velocity distribution’s high-energy tail and the quantum tunneling terms
  * Most likely energy ("**Gamow peak**") $E_0 = \left(bkT\over2\right)^{2\over3}$
* **Luminosity gradient** ${dL_r \over dr} = 4\pi r^2\rho$,
  * Where $L_r$ is luminosity enclosed in $r$,
    and $\epsilon$ is specific power (power released per mass)

### Nucleosynthesis

* *Notation: $^A_Z\text{He}$, where $A$ is atomic mass number and $Z$ is proton number*
* Simultaneous 4-body collision $4 {\;}^1_1\text{H} + ? \to {\;}^4_2\text{He} + ?$ extremely unlikely
* Reaction chain of 2-body interactions more probable
* Interactions must obey conservation laws:
  electric charge, nucleon number, and lepton number
* **Proton–proton chain** ($\text{H}$ burning): $4 {\;}^1_1\text{H} \to {\;}^4_2\text{He} + 2\text{e}^+ + 2\nu + 2\gamma$
  * Three branches (see Fig. 10.8)
  * In Sun, 69% PP I, 31% PP II, 0.3% PP III
  * Near $T = 1.5 \times 10^7 \text{ K}$, $\epsilon \approx \epsilon_0 \rho X_\text{H}^2 \left(T \over 10^6 \text{ K}\right)^4$
    * Where $\epsilon_0 = 1.08 \times 10^{-12} \text{ W}\;\text{m}^3\;\text{kg}^{-2}$
* **CNO cycle** ($\text{H}$ burning): $\text{C}$, $\text{N}$, and $\text{O}$ are catalysts
  * Two main branches, with the second only occuring ~0.04% of the time
  * Near $T = 1.5 \times 10^7 \text{ K}$, $\epsilon \approx \epsilon_0 \rho X_\text{H} X_\text{CNO} \left(T \over 10^6 \text{ K}\right)^{20}$
    * Where $\epsilon_0 = 8.24 \times 10^{-31} \text{ W}\;\text{m}^3\;\text{kg}^{-2}$
    * Much more temperature-dependent, only dominating in more massive stars
* **Triple-alpha process** ($\text{He}$ burning): $2 {\;}^4_2\text{H} \leftrightarrow {\;}^8_4\text{Be}$, ${\;}^8_4\text{Be} + {\;}^4_2\text{He} \to {\;}^{12}_6\text{C} + \gamma$
  * First step produces an extremely unstable $^8_4\text{Be}$,
    thus combined is essentially a 3-body interaction: $r \propto (\rho Y)^3$
  * $\epsilon \approx \epsilon_0 \rho^2 Y^3 \left(T \over 10^8 \text{ K}\right)^{41}$
    * Incredibly strong temperature dependence
* **Alpha process** ($\text{C}$ & $\text{O}$ burning): ${\;}^{12}_6\text{C} + {\;}^4_2\text{He} → {\;}^{16}_8\text{O} + \gamma, {\;}^{16}_8\text{O} + {\;}^4_2\text{He} \to {\;}^{20}_{10}\text{Ne} + \gamma$
  * Alpha capture becomes prohibitive at higher $Z$ due to higher Coulomb barrier
  * $\text{C}–\text{C}$ burning near $6 \times 10^8 \text{ K}$ and $\text{O}–\text{O}$ burning near $10^9 \text{ K}$ can produce
    $\{\text{ Na, Mg, Si, P, S }\}$; some are endothermic but are normally less likely
* Binding energy per nucleon ${E\over A}$
  * Relative to atomic number $A$, local maxima are very stable
  * Magic nuclei: some elements ($^4_2\text{He}, {\;}^{16}_8\text{O}$) have unusually high ${E\over A}$
  * Broad peak around $^{56}_{26}\text{Fe}$, the most stable nuclei

### Energy transport and thermodynamics

* Three mechanisms: radiation of photons (affected by opacity),
  convection, and conduction (generally insignificant)
* **Radiative temperature gradient** ${dT \over dr} = −{3\over4ac} {\kappa\rho \over T^3} {L\over4\pi r^2}$
* If temperature gradient becomes too steep, convection takes hold
* Convection is much more complicated than radiation
  * Strongly coupled to the star’s dynamic behavior
  * 3D Navier–Stokes with turbulence is hard compute
  * **Pressure scale height**, convection’s characteristic length scale,
    is big: $H \equiv −P {dr \over dP} \approx {P\over\rho g} \approx {R_*\over 10}$
* **First law of thermodynamics**: $dU = dQ − dW = dQ − P \;dV$
  * $U$ is a state function; $Q$ and $W$ are not — $dQ$ and $dW$ are inexact differentials
* **Specific heat capacity** $C \equiv {\partial Q \over \partial T}$, $C_P = C_V + nR$
  * $C_P$ is at constant pressure; $C_V$ is at constant volume
  * Heat capacity ratio / **adiabatic index** $\gamma \equiv {C_P \over C_V}$
  * $\gamma = {5\over3}$ for a monoatomic gas;
    approaches 1 in a partial ionization zone as both specific heats increase
* Isochoric process ($dV = 0$): $dU = dQ = C_V \;dT$
* Adiabatic process ($dQ = 0$): $dU = −P \;dV$
  * Gas law: $P \propto {1 \over V^\gamma} \propto \rho^\gamma \propto T^{\gamma/(\gamma−1)}$
  * **Speed of sound** $v = \sqrt{B\over\rho} = \sqrt{\gamma P \over \rho}$
    * Where bulk modulus $B \equiv −V {\partial P \over \partial V} \;\{dQ = 0\} $
* **Adiabatic temperature gradient** ${dT \over dr} = −\left(1 − {1\over\gamma}\right) {\mu \over k} {GM \over r^2} = −{g \over C_P}$
* If the surrounding temperature gradient is steeper than the bubble’s,
  even just slightly, the condition becomes **superadiabatic**,
  and nearly all luminosity is transferred outwards adiabatically,
  via convection instead of radiation
  * Equivalent criterion: ${d(\ln P) \over d(\ln T)} < {\gamma \over \gamma − 1}$, = 2.5 for ideal monoatomic gas
* In general, convection occurs if a region
  1. has high opacity (surrounding ${dT \over dr} \propto \kappa$),
  2. is ionizing and raising the specific heat capacity (bubble ${dT \over dr} \propto {1 \over C_P}$), or
  3. has a highly temperature-dependent fusion process
* First two conditions can occur simultaneously;
  third only occurs deep in the interior with the CNO or triple-alpha processes
* Convective flux under mixing-length theory:
  $F = \rho C_P {k\over\mu} \left({T \over g} \delta\left(dT \over dr\right)\right)^{3/2} \alpha^2 \sqrt\beta$
  * Where $0.5 \lesssim \alpha \lesssim 3$ and $0 \lesssim \beta \lesssim 1$ are free parameters
    * $\alpha \equiv {l\over H}$, the ratio of the mixing length
      (distance traveled by bubble before thermaliing with surrounding)
      and the pressure scale height
    * $\beta : \beta v^2$ is the average kinetic energy of the bubble as it travels over $l$

### Stellar model building

* Basic stellar models need constructive relations for $P$, $\bar\kappa$, and $\epsilon$:
  expressing them in terms of density, temperature, and composition
  * $P$ (pressure) generally modelable with ideal gas law and radiation pressure,
    but is more complex in certain stars’ deep interiors
  * $\bar\kappa$ (mean opacity) calculated explicitly from tabular data or fitting functions
  * $\epsilon$ (specific power) calculated analytically from reaction networks
* Boundary conditions:
  * As $r \to 0$, $\{M_r, L_r\} \to 0$
  * And ignoring extended atmospheres and mass loss:
    as $r \to R_*$, $\{T, P, \rho\} \to 0$
* **Vogt–Russell theorem**:
  due to a star’s dependence on nuclear burning,
  its mass and internal composition uniquely determine
  its radius, luminosity, internal structure, and subsequent evolution
  * Ignores smaller influences such as magnetic fields and rotation
* General modeling numeric integrates shell-by-shell
  with the system of differential equations,
  often from from a transition point towards both the surface and the center
* **Polytropes**: simplified stellar models in which $P(\rho) \propto \rho^\gamma$
* **Lane–Emden equation**: ${1\over\xi^2} {d \over d\xi} \left(\xi^2 {d(D_n) \over d\xi}\right) = −D_n^n$
  * Where the **polytropic index** $n : \gamma \equiv {n+1 \over n}$
  * TODO: too complicated

### Main sequence

* Vast majority of stars have $\text{H}$ mass fraction $X \approx 0.7$ and metal mass fraction $0 \lesssim Z \lesssim 0.03$
* Changes in core composition affects observed surface features
* Very light stars ($M \lesssim 0.08 M_\odot$) are not hot enough to let fusion
  stabilize against gravitational contraction
  * Highly opaque and fully convective
* Lower-initiation-energy PP chain dominates for low-mass stars
  * Shallow core thermal gradient leads to radiative core
  * High shell opacity leads to convective shell
* Highly-temperature-dependent CNO cycle dominates for high-mass stars ($M \gtrsim 1.2 M_\odot$)
  * Steep core thermal gradient leads to convective core
  * Low shell opacity leads to radiative shell
* Very massive stars ($M \gtrsim 90 M_\odot$) have rapid core thermal oscillations
  affecting fusion rates
* Very massive stars may also have radiation pressure exceed gas pressure at
  outer layers, with maximum stable luminosity given by **Eddington limit**:
  $L ≤ {4\pi GcM\over\bar\kappa}$
* Main-sequence lifetimes decrease with increasing luminosity

## 11. The Sun

* Loosely treated as two parts:
  an optically thin atmosphere and an optically thick core
* Transition is ~600 km thick

### Solar interior

* Hydrogen burning below $R \lesssim 0.3 R_\odot$; convection above $R \gtrsim 0.7 R_\odot$
* Heterogeneous composition due to nucleosynthesis, convection, and elmental diffusion
* $^4_2\text{He}$ is more abundant than $^1_1\text{H}$ below $R \lesssim 0.1 R_\odot$
* $^3_2\text{He}$ abundance peaks at the top of the hydrogen-burning region,
  where cooler temperatures slow the $^3_2\text{He}–^3_2\text{He}$ reaction
* Convection zone turbulence creates homogeneous composition
* Peak energy generation region is shell around $r \approx 0.1 R_\odot$
  * ${dL\over dr}$ affected by shell volume and fuel availability,
    both smaller at the very center, as well as temperature and pressure
* Surface \text{Li} abundance somewhat less than expected for the current solar model
* **Mikheyev–Smirnov–Wolfenstein effect**: neutrino oscillation between the
  $\text{e}^-$, $\mu^-$, and $\tau^-$ flavors explains the solar neutrino problem of only ⅓ as many
  $\text{e}^-$ neutrinos as expected being detected from the solar core

### Photosphere

* Where optical photons originate
* Starts 100 km below where $\tau_\text{500 nm} = 1$, extending ~600 km
* Temperature drops from ~9400 K at base to \gtrsim4400 K at top
* Continuum opacity partly due to $\text{H}^-$ near base of photosphere,
  as the far more abundant neutral $\text{H I}$ cannot contribute much to the continuum
* Absorption lines produced in higher, cooler, more opaque regions of the photosphere
* **Granulation**: patches of bright and dark regions (~700 km) at base of photosphere due to
  underlying convection zone; Doppler shifts (~0.4 km/s) cause wiggles in absorption lines
  * Characteristic lifetime (~10 minutes) corresponds to the time a
    convection eddy takes to rise and fall 1 mixing length
* **Differential rotation**: Doppler shifts at solar limb and solar oscillations
  show the solar rotation varies by latitude and by radius
  * Period of 25 days at equator lengths to 36 days at poles
  * **Tachocline**: base of convection zone ($\approx 0.65 R_\odot$) where differing rotation rates converge,
    resulting in strong shear theorized to create plasma which
    generate the solar magnetic field

### Chromosphere

* Where temperature begins rising again, from 4400 to ~10,000 K
* Starts 525 km above where $\tau_\text{500 nm} = 1$, extending ~1600 km
* Density and intensity $\approx 10^{-4}$ of the photosphere
* Low density and high temperature produce certain absorption and emission lines,
* As the blackbody continuum emission peaks ~500 nm, visible-spectrum emission lines are
  only clearly seen in a flash spectrum of the limb near total eclipse
* Strength of $\text{H}\alpha$ emission line allows filters to selectively observe the chromosphere structure
* **Supergranulation**: patches of ~30,000 km wide from underlying convection
* **Spicules**: vertical gas filaments extending upwards for ~10,000 km,
  ejecting mass at ~15 km/s

### Transition region

* Starts ~2100 km above $\tau_\text{500 nm} = 1$, extending to the corona
* Temperature rises rapidly to 10⁵ K in 100 km, then slowly to >10⁶ K
* Selectively observed in various UV bands (e.g. Lyman $\alpha$ at 20,000 K)

### Corona

* Faint region ($\lesssim 10^{15} \text{ particles/m}^3$) with vaguely defined outer boundary
* High-temperature, high-thermal-conductivity, approximately isothermal plasma
* Quiet corona near sunspot minimum (low solar activity)
  * More extended at equator than poles, consistent with nearly dipole magnetic field
* Active corona near sunspot maximum
  * More complex magnetic field shape and structure
* Essentially transparent to most EM radiation
* Not in local thermodynamic equilibrium thus no strictly definable temperature,
  but estimated to be $\gtrsim 2 \times 10^6 \text{ K}$
* **Parker wind model**: not in hydrostatic equilibrium as pressure does not vanish at infinity,
  implying solar wind
* Kontinuierlich/continuous **K corona**:
  * From free electron scattering of photospheric light primary between $1 \sim 2.3 R_\odot$
  * Spectral lines essentially blended into continuum from high-velocity Doppler broadening
* Fraunhofer **F corona**:
  * From  dust scattering of photospheric light beyond 2.3 R_\odot
  * Slower dust grains have less Doppler broadening and leave detectable Fraunhofer lines
  * Merges with zodiacal light from interplanetary dust
* Emission **E corona**:
  * From highly ionized coronal atoms, very rich in emission lines in the X-ray spectrum
  * Low number density enable forbidden transitions between metastable energy levels,
    as well as low-energy free-free transitions that produce radio waves
  * Radio waves also produced by relativistic electrons’ synchrotron radiation

### Solar wind

* Stream of escaped ions and electrons
* Deflects comets’ ion tails differently than pure radiation pressure on their dust tails
* Trapped into Van Allen belts by planetary magnetospheres and
  create aurorae upon reaching atmosphere
* Fast solar wind: ~750 km/s, produced from open magnetic fields
  * Associated with coronal holes — darker, cooler coronal regions
* Slow solar wind: ~300 km/s, produced by coronal streamers around closed magnetic fields;
  * Associated with X-ray bright spots from trapped spiralling charges
* **Heliopause**: outer limit of the Sun’s EM influence;
  where solar wind meets the interstellar medium and produces a **termination shock**
* **Heliosheath**: beyond the termination shock,
  particles slow, magnetic field strengthens, and density increases

### Magnetohydrodynamics (MHD)

* Longitudinal pressure waves propagate outward from top of convection zone
* Drastic drop in medium density turn waves supersonic, creating shock fronts
  that drastically heat up chromospheric gas as they dissipate
* Magnetic energy density and pressure $u = P = {B^2\over2\mu_0}$
* Transverse **Alfvén waves** propagate from oscillations in magnetic field lines
* Resistive Joule heating from electrical currents in Alfvén waves also
  contributes to temperature rise, in particular the steep gradient of the transition region
* Open magnetic field lines are dragged by stellar rotation across interplanetary space,
  slowing stars down significantly over their lifetimes

### Sunspots

* Zeeman splitting indicates strong magnetic fields,
  which inhibit convective motion below and create dark spots
* Umbra: darkest portion with vertical field lines; diameter \lesssim 30,000 km
* Penumbra: border with filament-like structure with field lines becoming horizontal
* Generally located in groups,
  with one dominant sunspot leading several in the direction of rotation
  * Lead sunspots in the same hemisphere always have the same polarity
    during an 11-year cycle
  * Opposite polarity in the other hemisphere
* Frequency follows 11-year cycle, starting from sunspot minimum
* **Butterfly diagram**: average latitude starts near ±40^\circ, migrating to the equator over cycle
* Solar polarity reverses during sunspot minimum, thus technically 22-year cycle
* Long-term variations include the Maunder minimum spanning 1645–1715

### Plages

* Chromospheric regions of higher density and bright H\alpha emission
* Located near sunspots, forming before they appear and vanishing after they disappear
* Caused by magnetic fields

### Flares

* Eruptive events releasing $10^{17}–10^{25} \text{ J}$ of energy over timespans from milliseconds to over an hour
* May reach 100,000 km in length
* Develop in sunspot groups with intense stored magnetic energy
* Reconnection of magnetic field lines creates a sheet of current in the plasma,
  Joule heating the gas up to 10⁷ K
* Charged particles are accelerated away from the reconnection point
  * $\text{H}\alpha$ line becomes locally in emission from ejected particles
    recombining by the base of the field lines
  * Solar cosmic rays from particles ejected towards outer space
* Nonthermal radio waves from synchrotron radiation around field lines
* Soft X-rays from high temperatures in loop below reconnection point
* Hard X- and gamma-rays from surface nuclear reactions,
  including **spallation** of heavy nuclei:
  $^1_1\text{H} + {\;}^{16}_8\text{O} \to ^{12}_6\text{C} + {\;}^4_2\text{He} + {\;}^1_1\text{H} + \gamma$

### Prominences

* **Quiescent prominence**:
  * Curtains of cooler (~8000 K) ionized gas collected along magnetic field lines
  * May be stable for weeks
  * Appear as dark filaments against the disk
* Eruptive/**active prominence**:
  * Suddenly destabilized magnetic field configuration
  * Energy converted into lifting prominence away from the Sun

### Coronal mass ejections (CMEs)

* Ejects $5 \times 10^{12} \sim 5 \times 10^{13} \text{ kg}$ of material at speeds of 400–1000 km/s
* ~1/day, more frequent near sunspot maximum
* 70% correlated with eruptive prominences; 40% correlated with flares

### Magnetic activity in other stars

* **Flare stars**: M-type main-sequence stars whose occasional rapid fluctuation in brightness
  may be due to large flares on their relatively dimr surface
* Starspots may be used to measure stellar rotation
* Magnetic field lines measured from Zeeman broadening correlate with luminosity variations

## 12. The Interstellar Medium and Star Formation

### Interstellar medium (ISM)

* Gas and dust between the stars
* 70% $\text{H}$, ~30% $\text{He}$, rest metals ($\text{C}$, $\text{Si}$, etc)
* Used in star formation; returned by stellar winds and explosive events

### Dust

* ~1% of molecular clouds by mass, but significant in light absorption and cloud chemistry
* Likely facilitate formation of many molecules besides $\text{H}_2$,
  including solid $\text{CO}$, $\text{H}_2\text{O}$, etc that give the grains icy mantles
* Formed in dense envelopes of very cool stars, from supernovae and stellar winds
* Unexpected abundance of large grains an area of active research
* Responsible for **interstellar extinction**, a $\lambda$-dependent increase in apparent magnitude
  from scattering and absorption of starlight: $A = \Delta m$
  * Approximately equal to the optical depth: $A \approx 1.086\tau$
  * Visual band extinction commonly used as reference: $A_V$
* **Mie theory** applies to IR–visible range
  * Shorter wavelengths are scattered more (though spectral lines go unaltered),
    causing interstellar reddening along direct line of sight,
    and blue tones in reflection nebulae
  * Extinction coefficient: $Q(\lambda) \equiv {\sigma(\lambda)\over\sigma_g}$
  * With dust grain radius: $r$ such that grain cross section: $\sigma_g = \pi r^2$
  * Where photon cross section approaches 0 for $\lambda \gg r$ and a constant for $\lambda \ll r$:
    $\sigma(\lambda) \propto {r^3\over\lambda}$ for $\lambda \gtrsim r$, and $\sigma(\lambda) \propto r^2$ for $\lambda \ll r$
* Color excesses relative to Mie prediction grow for shorter wavelengths
* Bump at 217.5 nm corresponding to resonance of graphite,
  along with IR emission bands corresponding to C–C and C–H bond vibrations,
  indicate presence of polycyclic aromatic hydrocarbons  (PAHs)
* Near-IR absorption bands indicate presence of silicate grains
* Slight, $\lambda$-dependent polarization due to anisotropic, somewhat aligned dust grains;
  alignment most likely due to interactions with a weak magnetic field

### Hydrogen

* **21-cm line**: anti-aligned spins ↑↓ has slightly less energy than aligned spins ↑↑
  * Very stable, “forbidden” transition that takes millions of years
    on average to occur per atom
  * Only possible (although still rare) for $\text{H I}$ in low-density diffuse ISM
  * Optically thin, thus optical depth $\propto H$ column density
* $\text{H}_2$ hard to directly observe except for rovibrational bands;
  tracer molecules and their isotopomers used instead
* Optically thick dust shields $\text{H}_2$ from UV photodissociation,
  as well as enhancing $\text{H}_2$ formation
* $\text{H I}$ generally proportional in number to dust when $A_V \lesssim 1$;
  displaced by $\text{H}_2$ when dust becomes optically thick
* Shells of $\text{H I}$ clouds surround molecular clouds of $\text{H}_2$

### Interstellar cloud structure

* Most diffuse ISM are hydrogen clouds of ground state H I,
  only absorbing UV photons and emitting in 21-cm line
* Diffuse\overtranslucent molecular clouds:
  irregularly shaped; similar in conditions to H I clouds;
  primarily atomic H with regional concentrations of H_2
* Giant molecular clouds (GMCs): enormous complexes with clumpy structure
  * Contains slightly denser dark cloud complexes, denser clumps and dense cores,
    and very dense, star-forming hot cores
* Bok globules:
  almost-spherical clouds located outside of larger molecular complexes;
  possibly dense cores stripped of surrounding gas by stellar winds

| thing                   | $T\;(\text{K})$ | $M\;(M_\odot)$ | $A_V$        | $n\;(\text{m}^{-3})$                     | $D\;(\text{pc})$ |
| ----------------------- | --------------- | -------------- | ------------ | ---------------------------------------- | ---------------- |
| diffuse molecular cloud | $15\sim50$      | $3\sim100$     | $1\sim5$     | $5 \times 10^8 \sim 5 \times 10^9$       | $1\sim10$        |
| giant molecular cloud   | $\sim15$        | $10^5\sim10^6$ | $\gtrsim1$   | $1 \times 10^8 \sim 3 \times 10^8$       | $\sim50$         |
| · dark cloud complex    | $\sim10$        | $\sim10^4$     | $\sim5$      | $\sim5 \times 10^8$                      | $\sim10$         |
| · clump                 | $\sim10$        | $\sim30$       | $\sim10$     | $\sim 1 \times 10^9$                     | $1\sim5$         |
| · dense core            | $\sim10$        | $\sim10$       | $\gtrsim10$  | $\sim 1 \times 10^{10}$                  | $\sim0.1$        |
| · hot core              | $100\sim300$    | $10\sim1000$   | $50\sim1000$ | $1 \times 10^{13} \sim 1 \times 10^{15}$ | $0.05\sim0.1$    |
| Bok globule             | $10\sim$        | $1\sim1000$    | $\sim10$     | $\gtrsim 1 \times 10^{10}$               | $\lesssim1$      |

### Interstellar cloud heating and cooling

* Primarily heated by cosmic ray protons
  * $E$ ranges $10–10^{14} \text{ MeV}$, with $10^3–10^8 \text{ MeV}$ common
  * Ionizes $\text{H}$ and $\text{H}_2$, ejecting $\text{e}^-$’s that distribute thermal energy
    throughout ISM via collisions with molecules
* Also heated by UV starlight ionizing $\text{C}$ and dust grains,
  X-ray starlight ionizing $\text{H}$, dust grains absorbing starlight,
  and shocks from supernovae and strong stellar winds
* Primarily cooled by IR radiation from post-collision de-excitations,
  particularly those of $\text{C}^+$ and $\text{CO}$ after colliding with $\text{H}$ and $\text{H}_2$, respectively

### Protostar formation

* Protostars: pre-nuclear-burning objects formed from molecular clouds
* Deviation from hydrostatic equilibrium means imbalance in the virial theorem:
  $2K < |U|$ ⇒ gravitational collapse
* Jeans criterion, minimum necessary for spontaneous collapse:
  * **Jeans mass** $M \approx \sqrt{ \left(5kT \over G\mu\right)^3 {3\over4\pi\rho_0} }$
  * **Jeans length** $R \approx \sqrt{ 15kT \over 4\pi G\mu\rho_0}$
  * Neglecting rotation, turbulence, magnetic fields,
    and pressure of surrounding gas
* Bonnor–Ebert criterion accounts for external pressure $P_0$
  * **Bonnor–Ebert mass**: $M = { c_BE v^4 \over \sqrt{P_0 G^3} }$
  * Where isothermal sound speed $v \equiv \sqrt{kT\over\mu}$,
    and constant $c_\text{BE} \approx 1.18$
* Low-pressure, optically thin first phase of collapse essentially isothermal and in free-fall
  * Free-fall timescale $t = \sqrt{3\pi \over 32G\rho_0}$
* **Homologous collapse**:
  Initially uniform spherical cloud would uniformly contract and condense
* Inside-out collapse:
  Somewhat centrally condensed cloud will have
  shorter free-fall time and thus condense faster near center
* Significant rotation results in disk-like structure
* Magnetic pressure keeps clouds balanced; eddy braking slows collapse
  * Critical mass with magnetic pressure: $M = {c_B \pi R^2 B \over \sqrt G}$
  * Supercritical mass achieved by merger of subcritical clouds,
    or (more commonly) regional lessening of magnetic field
* **Ambipolar diffusion**: neutral particles drift gravitationally,
  but are slowed by collisions with magnetically-frozen charged particles
* **Fragmentation**: increasing density with constant $T$ reduces Jeans mass,
  causing density inhomogeneities to collapse locally,
  forming large numbers of small protostars in place of a single massive one
* **Initial mass function** (IMF):
  number of stars formed per mass interval strongly mass-dependent,
  with abundance of low-mass stars
* Increasing density purely adiabatically increases Jeans mass
* Fragmentation stops as center becomes optically thick when $\sim10^{-10} \text{ kg}\;\text{m}^{-3}$,
  trapping radiation and making collapse more adiabatic
* Central region reaches near hydrostatic equilibrium and becomes protostar
  * Larger dust photosphere has radius and effective temperature where $\tau = {2\over3}$,
    putting star on H–R diagram
  * ~10 AU wide for collapse from $1 M_\odot$ cloud

### Protostar evolution

* Surrounding material still in free-fall until meeting protostellar core,
  releasing significant kinetic energy at shock front as they slow from
  supersonic speeds
  * Emits blackbody radiation predominantly in IR
  * Luminosity increases with temperature for several thousand years
* Dust begins to vaporize and opacity drops ~1000 K
  * Substatially reduces photosphere radius to nearly the surface of the
    hydrostatic core
  * Increases effective temperature with little change in luminosity
* $\text{H}_2$ dissociates ~2000 K,
  absorbing thermal energy otherwise used to maintain hydrostatic equilibrium,
  triggering second collapse
* Core radius decreases to $\sim1.3 R_\odot$ for collapse from $1 M_\odot$ cloud, while
  core mass $\ll 1 M_\odot$, indicating ongoing accretion at a second shock front
* Hot core begins burning deuterium ($^2_1\text{H}$), producing ~60% of luminosity,
  which stays roughly constant
* Deuterium burn-out leads to sharp drop in luminosity with small decrease
  in effective temperature, leading to pre-main-sequence star

### Protostar observation

* Direct observation shielded by thick dust of molecular cloud,
  and made rare by relatively short free-fall time
* Small IR sources embedded in dense cores or Bok globules indicate collapse
* Infalling material around embedded IR objects show twin spectral splitting
  from Doppler effect

### Pre-main-sequence evolution

* Increasing temperatures lead to $\text{H}^-$ opacity dominating in outer layers,
  creating envelope to become deeply convective
* **Hayashi track**: convection constrains H–R path to near vertical line
  * Luminosity decreases with slight increase in temperature
    as protostar slows collapse and reaches hydrostatic equilibrium
* Formation of $\sim1 M_\odot$ stars:
  * Completely convective from high $\text{H}^-$ opacity for the first ~10⁶ years of collapse,
    slowed little by scarce deuterium burning
  * Rising central temperature decreases opacity and creates expanding radiative core
  * Core begins generating energy from first steps in the
      PP I chain ($2 {\;}^1_1\text{H} \to {\;}^3_2\text{He}$) and CNO cycle ($^{12}_6\text{C} + {\;}^1_1\text{H} \to ^{14}_7\text{N}$),
  * Surface luminosity increases with temperature again and moves star off Hayashi track
  * Core burning comes to increasingly dominate luminosity generation
  * CNO reactions give core steep temperature gradient and some convection
  * Core expansion from intense nuclear energy production expends gravitational energy,
    lowering luminosity and effective temperature to main-sequence values
  * After exhaustion of $^{12}_6\text{C}$ for CNO, core reaches temperature
    for hydrostatic hydrogen burning via the full PP I chain
* Formation of $\lesssim 0.5 M_\odot$ stars:
  * Central temperature never reaches efficient $^{12}_6\text{C}$-burning temperatures
    for upward branch between Hayashi track and main sequence
  * Fully convection: temperature stays sufficiently low and opacity sufficiently high
    to never develop radiative core
* Formation of brown dwarfs ($0.013 \sim 0.072 M_\odot$):
  * Low nuclear burning rate cannot form main-sequence star
  * $\text{Li}$ burning $\gtrsim 0.06 M_\odot$
  * Deuterium burning $\gtrsim 0.013 M_\odot$ (~13 Jupiter masses)
  * Spectral types L and T
* Formation of massive stars:
  * High central temperatures leads to early, high-luminosity departure from Hayashi track
  * Evolve nearly horizontally across H–R diagram
  * Full CNO cycle becomes dominant H-burning mechanism,
    with steep thermal gradient keeping core convective
* **Zero-age main sequence** (ZAMS): when stars begin equilibrium hydrogen burning

### Stellar formation effects on medium

* **OB associations**: star groups dominated by O and B main-sequence stars
* Massive O and B protostars will first vaporize surrounding dust,
  then dissociate molecules, and finally ionize immediate surroundings
  into an $\text{H II}$ region within an $\text{H I}$ region
* **$\text{H II}$ regions**: $\text{H I}$ ionized into $\text{H II}$ by protostar UV radiation fluoresce in visible light
  during recombination, creating emission nebulae
  in the visible spectrum during recombination
  * **Strömgren radius** of $\text{H II}$ region: $r \approx \sqrt[3]{3N \over 4\pi\alpha n^2}$
  * Where $N$ is rate of ionizing photon production,
    $\alpha$ is likelihood of recombination, and n is number density of protons and electrons
* Radiation pressure from cluster of highly luminous O and B stars
  drives significant mass loss, disperses surrounding cloud (halting star formation),
  and weakens gravitational bound of star group (generally unbinding them)
* **Circumstellar disks**
  * Formed by spin-up of cloud transferring angular momentum away from collapsing star,
    possibly via stellar winds coupled to magnetic fields from within the convection zones
  * Continuous spectrum from reflection of protostar light
  * May be accretion or debris disks, possibly forming protoplanets
  * **Proplyds**: protoplanetary disks $\gtrsim 10^{25} \text{ kg}$
* **Herbig–Haro objects**:
  narrow beams of supersonic gas jets ejected from poles of young protostars
  * Collision with ISM results in excitations, producing bright emission lines
* **T Tauri stars**:
  low-mass ($0.5\sim2 M_\odot$) pre-main-sequence objects in intermediate phase
  between IR source and main-sequence star
  * Large irregular variations in luminosity with timescales ~ days
  * Strong H, Ca II, and Fe emission and Li absorption
  * Often forbidden [O I] and [S II] indicative of extremely low gas densities
  * Some lines have a P Cygni profile indicative of significant mass loss $\approx 10^{-8} M_odot/\text{year}$:
    blueshifted absorption trough before an emission peak and
    a redshifted emission tail
  * P Cygni profile sometimes inverts, indicating significant mass accretion in
    a very unstable environment
* **FU Orionis stars**:
  T Tauri stars undergoing extreme mass accretion ($\approx 10^{-4} M_odot/\text{year}$)
  and increase in luminosity ($\gtrsim4$ magnitudes)
  * Circumstellar disk instabilities dump $\approx 0.01 M_\odot$ of material
    onto central star over ~100 years
  * Inner disk outshines centural star by $100–1000\times$
  * May occur to a T Tauri star several times over its lifetime
* **Herbig Ae/Be stars**: type A and B pre-main-sequence stars with strong emission lines

### Modifications to classical model

* Classical model neglects rotation, turbulence, magnetic fields,
  initial inhomogeneities, strong stellar winds, ionizing stellar radiation,
  pressure-free protostellar collapse, likely smaller initial radii,
  and upper limits to massive star accretion,
* Birth line: evolutionary theories with smaller initial radii
  place upper limit on observed protostar luminosity
* Some observations suggest massive starts $\gtrsim 10 M_\odot$ may
  instead form from mergers or with accretion disk
  due to limiting feedback mechanism on non-rotational accretion
* Classical model predicts inverse relation of collapse time and mass,
  implying massive stars would disperse surrounding cloud before
  any low-mass stars can form

## 13. Main Sequence and Post-Main-Sequence Stellar Evolution

* Generally dominated by nuclear reaction timescale
* Kelvin–Helmholtz timescale relevant when transitioning between nuclear sources

### Main sequence (MS)

* Stars $0.3\sim1.2 M_\odot$:
  * Initially, core burns H, predominantly via PP chain
    * Radiative due to low temperature dependence
    * Core burning radius, density, and temperature rise due to fusion’s increase of
      mean molecular weight causing contraction and release of gravitational energy
    * Surface luminosity, radius, and effective temperature consequently rise over lifetime
  * Temperature high enough to burn thick H shell around core immediately after core H depletion
    * Core now isothermal and predominantly He, growing as He ash rains down from shell
    * Shell H burning generates more power than core H burning, but some of it is
      converted into gravitational potential through slow expansion of envelope
    * Effective temperature decreases slightly
  * **Schönberg–Chandrasekhar limit**, max mass fraction of an isothermal core
    in hydrostatic equilibrium supported by ideal gas pressure:
    ${M' \over M} \lesssim 0.37 \left(\mu\over\mu'\right)^2$
    * Where M′ and \mu′ are mass and mean molecular weights of the core
    * Derived from hydrostatic pressure: ${dP \over dM} = –{Gm \over 4\pi r^4}$
  * After reaching SC limit, core contracts on Kelvin–Helmholtz timescale
    unless mass is low enough to be supported by temperature-independent
    electron degeneracy pressure: $P \propto \rho^{5\over3}$
* Stars $\gtrsim 1.2 M_\odot$:
  * Initially, core nearly homogeneous due to convective mixing
  * Convection zone decreases in mass during H burning
    * Shrinks faster for more massive stars,
      disappearing entirely before H exhaustion for stars \gtrsim 10 M_\odot
    * Leaves slight composition gradient
  * After near-exhaustion of core H ($X \approx 0.05$ for $M \approx 5 M_\odot$), entire star contracts
    on Kelvin–Helmholtz timescale while increasing luminosity and effective temperature

### Subgiant branch (SGB)

* Low-mass stars:
  * Rapid core contraction releases gravitational energy,
    expanding envelope and decreasing effective temperature
  * Raised H shell temperature and density rapidly increases power generation,
    again expanding envelope
* Massive stars:
  * Rising temperatures eventually trigger rapid ignition of thick H shell around core,
    forcing slight expansion of star envelope,
    momentarily decreasing luminosity and temperature

### Red giant branch (RGB)

* Low- and intermediate-mass stars:
  * Near-surface convection zone created for
    due to $\text{H}^-$ ions increasing photospheric opacity
* All stars:
  * Star rises along Hayashi track as convection zone dominates star interior
  * First dredge-up: convection sinks elements from surface ($\text{Li}$, $^{12}_6\text{C}$, etc)
    and surfaces products of nuclear processes ($^3_2\text{He}$, $^{14}_7\text{N}$, etc)

### Red giant tip / start of helium fusion

* Stars $\lesssim 1.8 M_\odot$:
  * He core becomes strongly electron-degenerate
  * Significant neutrino losses create negative temperature gradient near center
  * Helium core flash: explosive initiation of triple-alpha process in He core in seconds
    * Core is rapidly lifted out of degeneracy outside-in, with
      strong temperature dependence of triple-alpha driving extreme thermal runaway
    * Briefly generates $\approx 10^{11} L_\odot$ of power, but most is absorbed in thermalizing core
* All stars;
  * Core expands from He burning, pushing H-burning shell outward
  * Luminosity abruptly decreases due to cooling of H-burning shell
  * Effective temperature begins to increase due to envelope contraction

### Horizontal branch (HB) / blue loop

* Stars $\lesssim 15 M_\odot$:
  * Envelope contraction raises effective temperature
    and compresses H-burning shell, increasing power generation
  * Deep convection zone rises towards surface while triple-\alpha makes core convective as well
  * Many develop instabilities in outer envelope, leading to periodic pulsations
  * Increase in mean molecular mass eventually causes core contraction,
    expanding envelope and lowering effective temperature
  * Depletion of core He rapidly accelerates core contraction
    * Core now mostly C and O, becoming very degenerate and
      cooling through significant neutrino production
  * Contraction and heating of surrounding He shell initiates He shell burning,
    in turn forcing expansion and cooling of surrounding H shell,
    temporarily halting H shell burning
* Stars $\gtrsim 15 M_\odot$ do not experience the blue loop

### Early asymptotic giant branch (E-AGB)

* All stars:
  * He-burning shell generates most power while H-burning shell is nearly inactive
  * Convective envelope absorbs energy and expands, lowering effective temperature
  * **Second dredge-up**: convection deepens again,
    bringing He and N from interior onto H-rich envelope
  * Star back on Hayashi track, approaching previous RGB path asymptotically from the left
  * Rapid mass loss due to low surface gravity

### Thermally-pulsating AGB (TP-AGB)

* All stars:
  * As He-burning shell nears exhaustion, H shell reignites
  * He ash from H-burning shell builds up and makes base of He shell slightly degenerate,
    eventually triggering helium shell flash
  * Convection zone established between He- and H-burning shells
  * He-burning expands and cools H-burning shell, gradually turning it off
  * Process repeats with growing amplitude after every pulse
  * Cycle evident from abrupt changes in surface luminosity
  * Period ranges from $\sim10^3 \text{ years}$ for stars $\sim5 M_\odot$ to $10^5 \text{ years}$ for stars $\sim0.6 M_\odot$
    * Long-Period Variables (LPVs) such as Mira have 100-700-day periods
* Stars $\gtrsim 2 M_\odot$
  * **Third dredge-up**: convection envelope merges with inter-shell convection zone
    and surfaces material from C-synthesizing region
  * Inverting number ratio of O and C in stellar atmosphere creates carbon stars,
    with spectral type C overlapping the traditional K and M
    * Spectral type S inbetween M and C has approximately equal O and C abundances
  * Mass loss in cool temperatures (~3000 K) expels into ISM
    silicate grains from O-rich atmospheres and graphite grains from C-rich atmospheres
  * **S-process nucleosynthesis**:
    nuclei in deep interior capture neutrons produced by nuclear burning,
    at a slow enough rate to radioactively decay before their next capture
  * S and C-type stars dredge up elements with no stable isotopes (e.g. Tc) into atmosphere,
    indicating active s-process

### Late and post-AGB

* Stars $\lesssim 8 M_\odot$:
  * He-burning increases mass of CO core, causing contraction until
    electron degeneracy pressure dominates
  * Mass loss accelerates with decreasing mass and increasing radius,
    developing $\sim10^{-4} M_\odot/\text{year}$ superwind
    * Precise mechanism unknown
    * Energizes shroud of optically thick dust clouds into OH\overIR sources:
      the $\text{OH}$ molecules emit IR as masers
  * Mass loss prevents catastrophic core collapse from Chandrasekhar limit,
    allowing stars between $4\sim8 M_\odot$ to additionally synthesize $\text{Ne}$ and $\text{Mg}$ in cores
  * Surrounding dust cloud eventually becomes optically thin from expansion,
    revealing an F or G supergiant
  * AGB phase ends as envelope is expelled, moving star nearly horizontally blueward
  * Luminosity drops rapidly as H and He-burning shells lose pressure and extinguish,
    revealing hot, degenerate C–O or ONeMg core: a white dwarf
  * Planetary nebula: expanding shell of gas around white dwarf
    * Emitting visible spectrum due to UV from central star remnant
    * Complex morphologies due to preferentially equatorial ejection
    * Receding at 10~30 km/s, with character length scales of ~0.3 pc
* Stars $\gtrsim 8 M_\odot$: see Ch. 15

### Stellar populations

* Stars give ISM back material enriched with more heavy elements than given
* **Population III**: stars formed soon after the Big Bang with virtually no metals ($Z = 0$)
  * Found in extreme deep field observations to primordial galaxies
* **Population II**: succeeding generations of metal-poor stars ($Z \gtrsim 0$)
  * Found outside galactic disk and in globular clusters
* **Population I**: current generations of metal-rich stars ($Z \approx 0.03$)
  * Found inside galactic disk and in open clusters

### Stellar clusters

* Stars formed from the same cloud, with similar compositions and birth times
  * Globular clusters: Population II, larger and older
  * Galactic \over open clusters: Population I, smaller and younger
* Graphed onto color-magnitude diagrams:
  H–R diagrams with B–V indices rather than effective temperatures,
  and apparent instead of absolute magnitudes if distance is not known
* Distance generally calculated with spectroscopic parallax
* Isochrone: curve fitting all stars of a cluster at a certain age
* Main-sequence turn-off point:
  color and magnitude where stars are currently leaving main sequence
  * Reflects age of cluster as point becomes redder and less luminous over time,
* Hertzsprung gap: absence of stars between late main sequence and red giant regions
  due to rapid Kelvin–Helmholtz-timescale intermediate evolution
* Blue stragglers: group of stars found above turn-off point, possibly
  due to mass exchange with binary companion or collision between stars

## 14. Stellar Pulsation

### Observations

* Pulsating stars dim and brighten as radius and temperature change
* **Long-Period Variables** (LPVs; prototype: $\omicron$ Ceti / Mira):
  thermally-pulsating asymptotic giants with 100–700-day periods and somewhat irregular light curves
  * May pulsate in either fundamental or first overtone mode
* **Classical Cepheids** (prototype: $\delta$ Cephei):
  supergiant Ib stars with 1–50-day periods proportional to their average luminosities
  * **Period–luminosity relation** (Leavitt’s law):
    $V = –2.81 \log_{10} P – 1.43$,
    $H = –3.234 \log_{10} P + 16.079$
    * Where $P$ is period in days
    * Measuring magnitude in the IR H-band mitigates some interstellar extinction
  * **Period–luminosity–color relation**:
    $H = –3.428 \log_{10} P + 1.54(J – K_s) + 15.637$
    * Where J–K_s is an IR color index
    * Adding color term improves data fit
  * Used as “**standard candles**” for measuring intergalactic distances
  * Luminosity variation primarily due to ~1000 K variation in temperature,
    with a phase lag of maximum luminosity occuring behind minimum radius
  * Vast majority pulsate in fundamental mode
* **W Virginis stars**:
  Population II Cepheids, around $4 \times$ less luminous than classical Cepheids of the same period
  * Vast majority pulsate in fundamental mode
* **RR Lyrae stars**:
  horizontal-branch Population II stars with 1.5–24-hour periods,
  all having nearly the same luminosity
  * May pulsate in either fundamental or first overtone mode
* **$\delta$ Scuti variables**:
  evolved F stars near the main sequence, with both radial and nonradial oscillations
  on 1–3-hour periods
* **ZZ Ceti stars**: pulsating white dwarfs with 100–1000-second periods
* **Instability strip**: narrow, nearly vertical region on H–R diagram,
  right of the main sequence, where most pulsating stars are found,
  including the above types
* **$\beta$ Cephei stars**: luminous (class III–V) blue variables with 3–7-hour periods;
  found in H–R’s upper left, outside the instability strip

### Radial pulsation mechanisms

* Oscillations result from standing sound waves in interior
  * Sound waves in medium with adiabatic index $\gamma$ travel at $v = \sqrt{\gamma P\over\rho}$
  * No displacement at nodes; maximum displacement at antinodes
* **Period–mean density relation**: $\Pi \approx \sqrt{3\pi\over2\gamma G\rho}$
  * Denser stars (e.g. white dwarfs) pulsate faster (than e.g. supergiants)
* Modes given by conical harmonics
  * Fundamental mode: node at star center, antinode at surface
  * Each overtone adds a node between center and surface
* Surface amplitudes decreases with overtone:
   ${\delta r\over R} \approx 0.07$ for fundamental, $\lesssim 0.01$ for first, $\approx 0$ for second
* Eddington modeled stars as heat engines:
  layers doing positive net work on surroundings drive oscillations;
  those doing negative net work dampen them
  * Positive work done if layer absorbs heat around max compression;
    releases heat and reaches max pressure during expansion
* **Nuclear $\epsilon$-mechanism**:
  compressing center raises temperature and density, increasing power generation
  * Amplitude usually too small to drive pulsation
  * May contribute to preventing formation of $\gtrsim 90 M_\odot$ stars
* Opacity $\kappa$ and $\gamma$ mechanisms:
  compressing layer increases opacity and traps heat, driving expansion
  * For most layers, opacity decreases with increased temperature from compression
  * **$\kappa$-mechanism**: in partial ionization zones however,
    some compression energy goes into further ionization instead of direct heating,
    letting opacity increase with the higher density
  * **$\gamma$-mechanism**: heat prefers to flow into partial ionization zone in compression
    due to its lower relative temperature, reinforcing the $\kappa$-mechanism
* Most stars have two main ionization zones
  * **Hydrogen partial ionization zone**: broader and closer to surface;
    ionizes $\text{H I} \to \text{H II}$ and $\text{He I} \to \text{He II}$ at characteristic temperature of $1–1.5 \times 10^4 \text{ K}$;
    changing depths during pulsation accounts for phase lag in classical Cepheids and RR Lyrae
  * **$\text{He II}$ partial ionization zone**: narrower and deeper;
    ionizes $\text{He II} \to \text{He III}$ at characterstic temperature of $4 \times 10^4 \text{ K}$;
    primarily responsible for driving oscillations within instability strip
* Temperature-dependent depths of ionization zones determine pulsation properties:
  * Blue edge of instability strip (~7500 K):
    high temperatures puts ionization zone too close to low-density surface,
    where there is insufficient mass to effectively drive oscillations
  * First overtone may be excited ~6500 K
  * Fundamental mode takes hold ~5500 K
  * Red edge of instability strip (~5000 K):
    heat transfer by convection bypasses high-opacity ionization zones
    and quenches pulsation
* $\beta$ Cephei stars pulsate due to iron ionization zone
  * Effective temperature (20,000~30,000 K) too high for $\text{H}$ and $\text{He}$ ionization zones
  * $\kappa$ and $\gamma$ mechanisms rely on $\text{Fe}$ opacity bump near $10^5 \text{ K}$

### Pulsation model

* Newton’s second law must be used instead of hydrostatic equilibrium model
  to account for oscillation of mass shells:
  $\rho {d^2r \over dt^2} = – {GM\rho \over r^2} – {dP \over dr}$
* Nonlinear evaluation can model complex, nonsinusoidal behavior of large-amplitude pulsations,
  but is very computationally expensive
* Linearizable by approximating with small amplitudes,
  but results in sinusoidal oscillations with no amplitude information
* Adiabatic approximation also used to minimize complexity,
  but nonlinear, nonadiabatic models are necessary for some variable stars
* One-zone linear, adiabatic model: $\Pi = {2\pi \over \sqrt{{4\pi \over 3} G\rho(3\gamma–4) }}$
* **Dynamic stability**: star collapses if $\gamma < {4\over3}$

### Nonradial pulsation mechanisms

* Some regions of surface expand while others contract
* Oscillations result from standing sound waves with latitudinal nodal circles
  and traveling sound waves with longitudinal nodal circles
  * No displacement at nodal circles
* **Angular modes** given by real parts of spherical harmonic functions: $Y^m_l(\theta, \phi)$
  * There are $l$ nodal circles: $|m|$ longitudes and $l – |m|$ latitudes
  * Where $l \in Z^+$ and $m \in [–l, l]$
  * Traveling waves take $|m|\cdot\Pi$ long to travel around star,
    with direction dependent on sign of $m$
* Pressure waves: **p-modes**
  * Confined to low depths, revealing conditions in stellar surface layers
  * Both radial and angular nodes
  * Acoustic frequency $S = \sqrt{\gamma P \over \rho} {\sqrt{l(l+1)} \over r}$
  * Frequencies split by prograde and retrograde wave motion
    proportional to star rotation rate: $\Delta S \propto m\Omega$
* Internal gravity waves: **g-modes**
  * Reveals movement of stellar material in deep interior
  * Only have angular nodes
  * **Brunt–Väisälä bouyancy frequency**: $N = \sqrt{–Ag}$
  * Confined to radiative zones, where $A < 0$: $A \equiv {1\over\gamma P} {dP \over dr} – {1\over\rho} {d\rho \over dr}$
* Surface gravity waves: **f-modes**
  * Frequency inbetween p- and g-modes

### Helioseismology and asteroseismology

* All observed solar oscillations are in the p-mode, with 3–8-minute periods and
  very short (high $l$) horizontal wavelengths
* Likely driven by turbulent energy of convection zone
* Latitide-dependent rotation rate revealed by m-dependent frequency splitting
* Depth-dependent rotation rate revealed by $l$-dependent attenuation below convection zone
* Thick convection zone prevents surface observation of of g-modes
* $\delta$ Scuti stars tend to pulsate in low-overtone radial modes, low-order p-modes,
  and possibly g-modes
* **Rapidly oscillating peculiar A stars** (roAp) primarily pulsate in higher-order
  p-modes, with the pulsation axis aligned with their strong magnetic fields
  instead of the rotation axis

## 15. The Fate of Massive Stars

### Post-main-sequence evolution:

| $M\;(M_\odot)$ | path to supernova            |
| -------------- | ---------------------------- |
| $\gtrsim85$    | O →  Of → LBV → WN → WC → SN |
| $40\sim85$     | O →  Of →  WN → WC → SN      |
| $25\sim40$     | O → RSG →  WN → WC → SN      |
| $20\sim25$     | o → RSG →  WN → SN           |
| $10\sim20$     | O → RSG → BSG → SN           |

* **Luminous blue variables** (LBVs; prototype: S Doradus):
  massive stars whose brightness suddenly erupt from time to time
  * Effective temperatures between 15,000~30,000 K, masses $\gtrsim 85 M_\odot$,
    and luminosities $\gtrsim10^6 L_\odot$ approaching the Eddington limit
  * Possible mechanisms behind behavior include
    envelope mass loss from temporarily exceeding Eddington limit,
    atmospheric pulsation instabilities, and binary companions
* **Wolf–Rayet stars** (WR):
  very hot, rapidly rotating stars with unusually strong broad emission lines and high mass loss
  * Effective temperatures between 25,000~100,000 K, masses $\gtrsim 20 M_\odot$,
    mass loss rates $\gtrsim 10^{-5} M_\odot/\text{year}$, and equatorial rotation speed ~300 km/s
  * Atypical spectral composition due to mass loss progressive stripping away outermost layers:
    WNs emission lines dominate in He and N; WCs in He and C; WOs in O
* Blue supergiants (BSG)
* Red supergiants (RSG)
  * Humphreys–Davidson luminosity limit: the most massive stars never evolve to RSG portion
    due to max luminosity cutoff
* Of stars: O supergiants with pronounced emission lines

### Supernovae

* Recorded Milky Way supernovae:
  SN 1006, SN 1054 with Crab supernova remnant,
  Tycho’s supernova SN 1572, and Kepler’s supernova SN 1604
* SN 1987A: first nearby supernova in age of modern astronomy
  * Occurred in Large Magellanic Clouds, 50 kpc from Earth
    *
* Classifications:
  * Type I: no H lines
    * Ia: Si II lines
    * Ib/Ic: no Si II lines
      * Ib: He lines
      * Ic: no He lines
  * Type II: H lines
    * II-P: plateau
    * II-L: linear

### Core-collapse supernovae

### Gamma-ray bursts

### Cosmic rays

## 16. The Degenerate Remnants of Stars

## 17. General Relativity and Black Holes

## 18. Closed Binary Star Systems

# III. The Solar System

# IV. Galaxies and the Universe

To be continued...

{% endraw %}
