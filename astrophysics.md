# Carroll & Ostlie, Introduction to Modern Astrophysics: Xing’s notes

## 1. The Celestial Sphere

### Altitude−Azimuth Coordinate System

* Local to observer location and day of year
* Altitude: h, measured from horizon towards zenith.
* Zenith: z = 90° − h
* Azimuth: A, measured from north on the horizon, clockwise, first towards
east on the horizon

### Equatorial Coordinate System

* Celestial equator: plane passing through Earth’s equator
* Ecliptic: path of sun across seasons
* Vernal equinox: When ecliptic/sun crosses the celestial equator northward
* Declination δ: Measured in degrees N/S of celestial equator
* Right ascension α: Measured in hour-minute-seconds from vernal equinox
  eastward to hour circle of object
* Meridian: great circle intersecting observer and the poles
* Hour angle H: angle between object and meridian of observer
* Local sidereal time: time elapsed since vernal equinox last traversed
  the meridian; hour angle of that intersection

### Precession

* ~25,770 year cycle of wobbling
* Necessitates use of epochs in the coordinate system
* J2000.0: Standard epoch Noon 2000−01−01 UT
  (universal time, measured from Greenwich)
* Uses Julian Calendar (365.25 d/yr), without Gregorian corrections

### Time standards

* Julian Date (JD): Days since Noon 4713BC−01−01
* Heliocentric Julian Date (HJD): From center of Sun
* Terrestrial Time (TT): From surface of Earth factoring in relativity

### Motions

* Radial velocity: v_r, positive is away
* Proper motion: μ, arcseconds/year
* Transverse/tangential velocity: v_θ = μ * r

### Spherical Trigonometry

* Sides: a, b, c
* Angles: A, B, C
* Law of sines: (sin a / sin A) = (sin b / sin B) = (sin c / sin C)
* Law of cosines for sides: cos a = cos b cos c + sin b sin c cos A
* Law of cosines for angles: cos A = − cos B cos C + sin B sin C cos a
* For a proper motion with position angle Φ (direction measured) and
  arclength Δθ,
  * Δα = Δθ sin Φ / cos δ
  * Δδ = Δθ cos Φ
  * (Δθ)² = (Δα cos δ)² + (Δδ)²

## 2. Celestial Mechanics

### Ellipses

* Semi−major axis: a
* Semi−minor axis: b
* Eccentricity: e | ae = d/2, where d = distance between foci
* b² = a² (1 − e²)
* Polar equation: r = a (1 − e²) / (1 + e cos θ)
* Area: A = πab

### Kepler’s Laws

1. r = (L²/μ²)/GM(1 + e cos θ)
2. dA/dt = L/2μ, constant
3. P² = ka³, k = 4π²/G(M+m), k = 1 using years and AU

### Orbital equations

* Virial theorem: 〈E〉 = 〈U〉/2, 〈U〉 = −2〈K〉
* Vis−viva: v² = G(M+m)[(2/r) − (1/a)]

## 3. Continuous Spectrum of Light

### Parallax and parsecs

* Parallax angle: p″ (arcseconds),
  subtended by quarter orbit of Earth (1 AU)
* Distance in AU: d = 206,265 / p″
* Distance in parsecs (pc): d = 1 / p″

### Apparent magnitude

* Originally compiled by Hipparchus: 1 for brightest, 6 for dimmest
* Modern definition: F₂/F₁ = 100^[(m₁ − m₂)/5]
* Where radiant flux is density of emitted light: F = L/4πr²

### Absolute magnitude

* M, apparent magnitude of star if it’s located at 10 pc
* H for planets, measured if located at 1 AU
* Distance modulus: d = 10^[(m − M + 5)/5],
m − M = 5 log₁₀(d) − 5
* For equally distant stars, ratio of radiant fluxes equals
  ratio of luminosities, thus: L₂/L₁ = 100^[(M₁ − M₂)/5],
* M_sun = +4.74, L_sun = 3.839E26 W

### Light

* c = λν
* Poyinting vector: S = E×B / μ₀, 〈S〉= E₀B₀ / 2μ₀
* Radiation pressure: P = 〈S〉cos(θ) / c, absorption;
P = (2〈S〉A cos²θ)/c, reflection

### Blackbody radiation

* Wien’s displacement: λT = 0.002898 m·K, for peak wavelength
* Stefan−Boltzmann: L = 4πR²σT⁴,
  where T is effective temperature of the star’s surface,
  and σ = 5.670 × 10⁻⁸ W/m²K⁴
* Planck’s function: B(T) = (2hc²/λ⁵)/(e^(hc/λkT) − 1)

### Color index

* Bolumetric magnitudes: measured over all wavelengths
* In practice, detectors only measure certain ranges, to varying
  degrees of sensitivity
* Standard filters:

|   | name | center | bandwidth |
| - | ---- | ------ | --------- |
| U |  UV  | 365 nm |   68 nm   |
| B | blue | 440 nm |   98 nm   |
| V |visual| 550 nm |   89 nm   |

* Isolated letters are apparent; M_letter is absolute
* U−B color index: U − B = M_U − M_B
* B−V color index: B − V = M_B − M_V, smaller is bluer
* Bolumetric correction: BC = m_bol − V = M_bol − M_V
* Color−color diagrams graph U−B against B−V,
  showing stars are non−ideal blackbodies

## 4. Special Relativity

* (Notation: prime (′) means observed, unprimed means rest)

### Einstein’s postulates:

1. Principle of Relativity: Laws of physics are the same in all inertial
  reference frames
2. Constancy of the Speed of Light: In vacuum, c is independent of the motion
  of the light source

### Lorentz transformations:

* Lorentz factor: γ ≡ 1 / √(1 − u²/c²)
* Time dilation: t′ = γt
* Length contraction: x′ = x/γ
* Proper time (τ) and proper length measured at rest w.r.t. events, γ = 1

### Relativistic doppler shifts

* Doppler shift: ν′ = ν / [ γ (1 + (u cos θ / c)) ],
  θ = 0° away from observer, 180° towards observer
* Redshift: z = (λ′/λ) − 1 = √[ (1 + (u cos θ / c)) / (1 − (u cos θ / c)) ]
* Time dilation from redshift: z + 1 = t′/t

### Relativistic momentum and energy

* Rest mass: m; Rest energy: mc²
* Relativistic momentum: p = γmv
* Relativistic kinetic energy: K = mc²(γ − 1)
* Relativistic total energy: E = γmc² = √( p²c² + m²c⁴ )

## 5. Interaction of Light and Matter

* Fraunhofer lines: absorption lines in the Solar spectrum

### Kirchoff’s Laws

1. Hot, dense gas/solid  produces continuous, lineless spectrum
2. Hot, diffuse gas produces emission lines
3. Cold, diffuse gas in front of a continuous spectrum source produces
   absorption lines in the spectrum

### Spectrographs

* Produced via diffraction grating: d sin θ = nλ,
  d is grating separation; n is spectral order
* Resolving power: λ/Δλ = nN,
  Δλ is smallest resolvable wavelength; N is # of gratings illuminated
* Low−speed Doppler shift approximation: Δλ/λ = (u cos θ)/c

### Photoelectric effect

* Shining light onto metal ejects e− with energy ∝ frequency, not intensity
* Max kinetic energy of an ejected e−: K = hν − Φ
* Work function Φ is the metal’s min e− binding energy
* Demonstrates quantization of light into packets of energy hν

### Compton scattering

* Relativistic scattering of electrons by high−energy photons
* Photon loses wavelength Δλ = h(1 − cos θ)/mc, e− mass m = 9.109⁻³¹ kg

### Bohr’s semiclassical atom

* Assumes quantization of angular momentum: L = μvr = nħ,
  n is principle quantum number
* However, still puts electrons in classical circular orbits modeled by
  coulombic attraction
* Kinetic energy of atom: K = μv²/2 = e²/8πε₀r = nħ²/2μr²
* Allowed orbital radii: r = a₀n²,
  Bohr radius a₀ = 4πε₀ħ²/μe² ≈ 5.292⁻¹¹ m
* Allowed energy levels: E = −E₀/n², E₀ ≈ 13.6 eV

### deBroglie waves

* Generalizes photon wavelength-frequency-momentum relation to all partilces
* ν = E/h, λ = h/p

### Heisenberg’s Uncertainty Principle

* Δx Δp ≥ ħ/2, estimated as Δx Δp ≅ ħ, equivalent to ΔE Δt ≈ ħ
* Reduced Planck constant: ħ = h/2π = 1.054 × 10⁻³⁴ J·s = 6.582 × 10⁻¹⁶ eV·s
* Allows for quantum tunnelling when barrier is within a few wavelengths wide
  * e.x. total internal reflection can be disrupted by placing a another    prism sufficiently close to the boundary

### Schrödinger's quantum atom

* Schrödinger equation models particles as probability waves
* Analytically solvable for H with same result as the Bohr atom
* Quantum orbitals are probability density clouds
* Adds additional quantum numbers l (0..n-1) and m_l (−l..+l)
* Allowed angular momentum magnitudes: L = ħ √[ l (l + 1) ]
* Allowed angular momentum z-components: L_z = ħ m_l
* Different orbitals with same n are degenerate
  (i.e. have the same energy) *unless* under a magnetic field

### Normal Zeeman effect

* Previously degenerate orbitals split into slightly different energies
* 3 frequencies: ν ∈ { ν₀, v₀ ± eB/4πμ }
* Allows detection of magnetic field strength around e.g. starspots
* Even if field is too weak to cleanly split,
  different polarizations are still detectable

### Anomalous Zeeman effect and Spin

* Anomalous Zeeman effect: more complex splitting patterns
* Resulted in discovery of the spin quantum number
* Allowed spin angular momentum magnitudes: S = ħ √3/2
* Allowed spin angular momentum z-component: S_z = ħ m_s

### Fermions and bosons

* Dirac's relativistic Schrödinger's equation divided particles
  into fermions and bosons
  * Also predicted antiparticles: opposite electric changes and magnetic moments
* Fermions have spins of odd * ħ/2; Bosons have spins of integer * ħ
* Pauli exclusion principle: No two fermions can share the same set
  of 4 quantum numbers (n, l, m_l, m_s)

### Complex spectra

* Quantum state transitions go from one set of quantum numbers to another
* Selection rules significantly lower the probabilities of certain transitions
* Restrictions: Δl = ±1, Δm_l ∈ { 0, ±1 }, m_l of orbitals cannot both be 0
* Allowed transitions triggered by collisions or occur spontaneously ~10⁻⁸ s
* Some forbidden transitions occur in very low gas densities
  * e.x. diffuse interstellar medium; outer stellar atmospheres

## 6. Telescopes

### Basic optics

* Index of refraction: n ≡ c/ν, wavelength dependent
* Snell's law: n₁ sin θ₁ = n₂ sin θ₂
* Lensmaker's formula: P = 1/f = (n - 1)[ (1/R₁) + (1/R₂) ],
  P is optical power, f is focal length, R are radii of curvature
* Focal length is wavelength-independent for mirrors
* A lens projects angular separation onto a plate (focal plane)
* Linear separation of point sources on focal plane increases with focal length:
  dθ/dy = 1/f

### Resolution

* Single-slit diffraction of wavefronts through aperture limits resolving power
* Airy disk: bright spot around central maximum
* Destructive interference creates minimums: sin θ = mλ/D, m = 1, 2, 3, ...
* Images are “unresolved” when one's central maximum overlaps with
  another's 1st minimum
* Rayleigh criterion defines arbitrary min resolvable distance
  of a circular aperture: θ = 1.22λ/D, D is aperture diameter
  * Possible to differentiate sources within Rayleigh criterion through careful
         analysis of diffraction patterns
* Actual ground-based optical resolution worse than ideal due to
  refraction from atmospheric turbulence, unless corected by adaptive optics
* Only lenses suffer from chromatic aberration;
  somewhat reduceable with correcting lenses
* Spherical lenses and mirrors suffer from spherical aberration;
  mitigated by paraboloids, which are harder to produce
* Paraboloid lenses and mirrors suffer from coma: elongation of off-axis images
* Astigmatism: different parts of a lens or mirror converge an image
  at different locations; correction may result in curvature of field issues

### Brightness

* Illumination: J ∝ 1/F², the amount of light energy per second focused onto a unit
  area of the resolved image
* Focal ratio: F ≡ f/D
  * Larger aperture increases resolution and illumination
  * Longer focal length increases image size but decreases illumination
  * For fixed focal ratio, larger telescope diameter increases resolution but
         not illumination

### Optical telescopes

* Refracting:
  * Objective lens focus as much light as possible onto the focal plane
  * Sensor is placed on focal plane or eyepiece is placed at its focal length away
  * Angular magnification: m = f_objective / f_eye
  * Issues:
    * Lens can only be supported near the edge, so gravity deforms heavy ones
         * Entire volume must be fabricated precisely
         * Slow thermal response and thermal expansion
         * Chromatic aberrations
* Reflecting:
  * Replaces refractor’s objective lens with a mirror
  * Minimizes weight, support, and distortion issues,
    esp. with active support system
  * Main issue: Prime focus of the mirror is in the path of the incoming light
    * Newtonian: flat diagonal secondary mirror;
           suffers from faraway eyepiece that introduces torque
    * Cassegrain: parabolic (or hyperbolic for Ritchey-Chrétien) primary mirror;
           secondary mirror (usually convex to increase focal length) reflects back
                through a hole in the primary’s center
    * Coudé: mirror system directs light to an instruments room;
           very long focal length
    * Schmidt: spheroidal primary mirror minimizing coma;
           correcting lens to remove spherical aberrations;
                provides wide field of view (degrees vs arcminutes) with low distortion
* Mounts:
  * Equatorial: polar axis, easy to adjust, hard to build for massive telecopes
  * Altitude-azimuth: easy to build, generally needs computer adjustments
* Adaptive optics: small deformable mirror with many piezoelectric actuators
         counteract atmospheric distortions, constantly calibrated via a guide star
* Charge-coupled device (CCD) detects ~100% of incident photons with a linear
  response across wide wavelength and intensity (dynamic) ranges

### Radio telescopes

* Strength measured is spectral flux density: S(ν)
* Typical radio source has S ≈ 1 jansky (Jy) = 10⁻²⁶ W/(m²·Hz);
  weak sources ≈ mJy
* Integrate S over area and bandwidth for total power received
* Long wavelengths require large apertures, but less manufacturing precision
* Addition of other telescopes enables interferometry;
  reduces side lobes and narrows main lobe of antenna sensitivity pattern
* Interferometry determines angle of source from phase difference between antennae
* Pointing angle: θ | sin θ = L/d,
  L is difference in wavefront’s distances to the antennae;
  d is baseline distance between antennae
* Very long baseline interferometry can span multiple continents

### IR/UV/X/γ astronomy

* Water vapor is primary contributor to IR absorption; thus low-humidity
         mountain peaks, balloons and aircraft observatories are used
* IR observation also requires very cold telescopes and detectors
* Glass is opaque to UV; UV telescopes need very precise reflecting surfaces
* X-ray and gamma-ray photons penetrate traditional mirrors;
  instead imaged by graze-incidence (near 90° incidence) reflections, or by
  Bragg scattering through crystal lattice inteference

# Part II. The Nature of Stars

## 7. Binary Systems and Stellar Parameters

* In many cases, allows the calculation of stellar masses, highlighting a
  well-defined mass-luminosity relation for the large majority of stars

### Classifications

* Optical double: Fake binaries that just look close
* (The rest of the classes are not mutually exclusive)
* Visual: Both stars are resolved
* Astrometric: Visual oscillatory motion of one star implies companion
* Eclipsing: Periodic variation in light curves reveal two stars
* Spectrum: Different spectral classes / Doppler effect reveal superimposed spectra
* Spectroscopic: Doppler effect in spectra reveal oscillation of
  radial velocity curves

### Mass determination of visual binaries

* Mass ratio from ratio of subtended angles: m₁/m₂ = r₂/r₁ = α₂/α₁
* If distance to system and inclination are known, masses can be determined:
  m₁ + m₂ = (4π²/G) (d / cos i)³ (α³/P²),
  where i is the angle of inclination between the orbit and the sky’s planes

### Mass determination of spectroscopic binaries

* Mainly detected from binaries with high inclinations
* “Double-line” if both spectra are visible; “single-line” otherwise
* Radial velocity curves cross at the radial velocity of the center of mass
* Radial velocity curve amplitudes scaled uniformly by sin i
* Eccentricity skews the curves, but close binaries tend to quickly
  circularize from tidal forces
* Mass ratio: m₁/m₂ = v₂/v₁, where v are radial velocities
* If inclination is known, masses can be determined:
  m₁ + m₂ = (P/2πG) (v₁ + v₂)³/(sin³ i),
* If only one radial velocity is known (single-line binary),
  then mass ratio cannot be found,
  but the mass function still provides rough constraints:
  (m₂ sin i)³/(m₁ + m₂)² = (P/2πG) v₁³

### Temperature ratio and radii determination of eclipsing binaries

* Two dips in the light curve:
  * Primary eclipse: hotter star passes behind cooler one
  * Secondary eclipse: cooler star passes behind warmer one
  * Size doesn’t matter
* If smaller star (T₁) is hotter, temperature ratio from brightness is:
  (B₀ − B₁)/(B₀ − B₂) = (T₁/T₂)⁴
* Unless extremely close binary, inclination must be ~90°,
  thus allowing accurate mass and velocity determination
* Radius of star: r = vΔt/2,
  where v is relative velocity,
  and Δt is transit time for the smaller star or eclipse time for the larger

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

### Harvard spectral types

* Pickering and Fleming initially labeled alphabetically by strength of H absorption lines
* Different strengths due to different temperatures causing different electron
  ionization levels and orbitals
  * H I (neutral H)’s visible spectral (Balmer) lines strongest at A0
    (effective T = 9520 K)
  * He I (neutral He)’s visible spectral lines strongest at B2
    (effective T = 22,000 K)
  * Ca II (singly ionized Ca)’s visible spectral lines strongest at K0
    (effective T = 5250 K)
* Cannon rearranged by temperature (OBAFGJM, O hottest)
  and added decimal subdivisions: A0 hotter than A9
* Additional spectral types for very cool stars and brown dwarfs: OBAFGKM + LT

### Spectral physics

* Statistical mechanics studies macroscopic behavior of stellar atmospheres
* Maxwell−Boltzmann velocity distribution: ρ(v) = √(m/2πkT)³ e^(−mv²/2kT) 4πv²
  * Restricted to gasses in thermodynamic equilibrium with density less than ≈1 kg/m³
  * Most probable speed: v = √(2kT/m)
  * Root-mean-square speed: v = √(3kT/m), due to right-skewed exponential tail
* Electrons’ orbital energies are affected by their atoms’ kinetic energies
  through collisions: higher orbitals are less likely to be occupied
* Probability ratio of electron states: P(s₁)/P(s₂) = e^[−(E₁−E₂)/kT]
  * As thermal energy kT → 0, P(higher energy state) → 0; confined to lowest state
  * As kT → ∞, probability ratios of states → 1; all states are equally likely
* Statistical weights of some energy levels increased by degenerate orbitals
* Boltzmann equation, the probability ratio of electron energies:
  P(E₁)/P(E₂) = (g₁/g₂) e^[−(E₁−E₂)/kT],
  where g is the statistical weight, or number of states with that energy
* Balmer lines require exitations from the first excited state N₂
* Greater fraction of H I are in the excited state N₂ at higher temperatures;
  However, at the same time, significant fractions of H I are ionized to H II
* Partition function, the number of possible electron states weighted by energy,
  varies by ionization state of the atom, thus affecting the probability of
  atoms at different ionization states: Z = Σ_j[ g_j e^[−(E_j−E₁)/kT] ]
  * For H I, the ground state (with two orbitals s=±½) dominates for most temperatures: Z ≅ 2
  * For H II, no electrons means only one possible configuration: Z = 1
* Saha equation, the probability ratio of ionization states:
  N₊/N = (2Z₊/nZ) √(2πmkT/h²)³ e^(−χ/kT), where m is e− mass
  * Very sensitive to ionization energy χ, as kT only ranges around 0.5−2 eV
  * Also effected by free e− density n: presence of ionized He II&III increases free e− density,
* Combining the Boltzmann and Saha equations shows varying narrow partial ionization
  zones for different atoms and ionization states
  * For H, there is significant fraction of electrons in the first excited state for
    temperatures between 8300−11300 K, matching the peaking of Balmer lines ≅ 9900 K

### Hertzsprung−Russell diagram

* Hertzsprung discovered type G and cooler stars had a range of magnitudes;
  termed the brighter ones giants; Russell termed the dimmer ones dwarfs
* H−R diagram plots absolute magnitude (brighter upward) against
  spectral type (warmer leftward)
* 80−90% of stars lie as dwarfs in a main sequence from upper left to lower right
  * Width of main sequence due to varying ages and compositions
* Simple luminosity-temperature relation reveals fundamental dependence on mass
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
