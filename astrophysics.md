# Carroll & Ostlie, Introduction to Modern Astrophysics: Xing’s notes

* _Notation: (a/bc)d = (a/(b·c))·d_

## 1. The Celestial Sphere

### Altitude–Azimuth Coordinate System

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
* M☉ = +4.74, L☉ = 3.839 × 10²⁶ W

### Light

* c = λν
* Poyinting vector: S = E×B / μ₀, 〈S〉= E₀B₀ / 2μ₀
* Radiation pressure: P = 〈S〉cos(θ) / c, absorption;
P = (2〈S〉A cos²θ)/c, reflection

### Blackbody radiation

* Wien’s displacement: λT = 0.002898 m·K, for peak wavelength
* Stefan–Boltzmann: L = 4πR²σT⁴
  * Where T is effective temperature of the star’s surface,
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

* _Notation: prime (′) means observed, unprimed means rest_

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

* Shining light onto metal ejects e⁻ with energy ∝ frequency, not intensity
* Max kinetic energy of an ejected e⁻: K = hν − Φ
* Work function Φ is the metal’s min e⁻ binding energy
* Demonstrates quantization of light into packets of energy hν

### Compton scattering

* Relativistic scattering of electrons by high−energy photons
* Photon loses wavelength Δλ = h(1 − cos θ)/mc, e⁻ mass m = 9.109⁻³¹ kg

### Bohr’s semiclassical atom

* Assumes quantization of angular momentum: L = μvr = nħ,
  n is principle quantum number
* However, still puts electrons in classical circular orbits modeled by
  Coulombic attraction
* Kinetic energy of atom: K = μv²/2 = e²/8πϵ₀r = nħ²/2μr²
* Allowed orbital radii: r = a₀n²,
  Bohr radius a₀ = 4πϵ₀ħ²/μe² ≈ 5.292⁻¹¹ m
* Allowed energy levels: E = −E₀/n², E₀ ≈ 13.6 eV

### De Broglie waves

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
  weak sources ~mJy
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
  m₁ + m₂ = (4π²/G) (d / cos i)³ (α³/P²)
  * Where i is the angle of inclination between the orbit and the sky’s planes

### Mass determination of spectroscopic binaries

* Mainly detected from binaries with high inclinations
* “Double-line” if both spectra are visible; “single-line” otherwise
* Radial velocity curves cross at the radial velocity of the center of mass
* Radial velocity curve amplitudes scaled uniformly by sin i
* Eccentricity skews the curves, but close binaries tend to quickly
  circularize from tidal forces
* Mass ratio: m₁/m₂ = v₂/v₁
  * Where v are radial velocities
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
* Radius of star: r = vΔt/2
  * Where v is relative velocity,
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

* Stars emit a continuous spectra (the continuum)
  with absorption lines in certain wavelengths

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
* Cannon rearranged by temperature (OBAFGKM, O hottest)
  and added decimal subdivisions: A0 hotter than A9
* Additional spectral types for very cool stars and brown dwarfs: OBAFGKM + LT

### Spectral physics

* Statistical mechanics studies macroscopic behavior of stellar atmospheres
* Maxwell−Boltzmann velocity distribution: ρ(v) = √(m/2πkT)³ e^(−mv²/2kT) 4πv²
  * Restricted to gasses in thermodynamic equilibrium with density ≲ 1 kg/m³
  * Most probable speed: v = √(2kT/m)
  * Root-mean-square speed: v = √(3kT/m), due to right-skewed exponential tail
* Electrons’ orbital energies are affected by their atoms’ kinetic energies
  through collisions: higher orbitals are less likely to be occupied
* Probability ratio of electron states: P(s₁)/P(s₂) = e^[−(E₁−E₂)/kT]
  * As thermal energy kT → 0, P(higher energy state) → 0; confined to lowest state
  * As kT → ∞, probability ratios of states → 1; all states are equally likely
* Statistical weights of some energy levels increased by degenerate orbitals
* Boltzmann equation, the probability ratio of electron energies:
  P(E₁)/P(E₂) = (g₁/g₂) e^[−(E₁−E₂)/kT]
  * Where g is the statistical weight, or number of states with that energy
* Balmer lines require exitations from the first excited state N₂
* Greater fraction of H I are in the excited state N₂ at higher temperatures;
  However, at the same time, significant fractions of H I are ionized to H II
* Partition function, the number of possible electron states weighted by energy,
  varies by ionization state of the atom, thus affecting the probability of
  atoms at different ionization states: Z = Σ_j[ g_j e^[−(E_j−E₁)/kT] ]
  * For H I, the ground state (with two orbitals s=±½) dominates for most temperatures: Z ≅ 2
  * For H II, no electrons means only one possible configuration: Z = 1
* Saha equation, the probability ratio of ionization states:
  N₊/N = (2Z₊/nZ) √(2πmkT/h²)³ e^(−χ/kT)
  * Where m is e⁻ mass
  * Very sensitive to ionization energy χ, as kT only ranges around 0.5−2 eV
  * Also effected by free e⁻ number density n: presence of ionized He II&III increases free e⁻ number density,
* Combining the Boltzmann and Saha equations shows varying narrow partial ionization
  zones for different atoms and ionization states
  * For H, there is significant fraction of electrons in the first excited state for
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
  I(λ) ≡ ∂I/∂λ ≡ [ (∂E/∂λ) dλ ] / [ dλ dt dA cos θ dΩ ],
  power transmitted at a certain wavelength within a certain solid angle
* Mean intensity:〈I(λ)〉= I(λ) for isotropic field; = B(λ) for blackbody
* Specific energy density: u(λ) dλ = (4π/c)〈I(λ)〉dλ
  * Total energy density for blackbody: u = 4σT⁴/c
  * Blackbody radiation pressure is u/3
* Specific radiative flux:
  F(λ) dλ = ∫ I(λ) dλ cos θ dΩ,
  net energy at a certain wavelength flowing towards one single direction
* For a resolved source, specific ray intensity is measured,
  thus detector distance does not affect measured intensity,
  though the angular size decreases
* For an unresolved source, specific radiative flux is measured,
  falling off with 1/r²

### Temperature

* Mean free path: average distance traveled by particles and photons
  between collisions
* Multiple descriptions of temperature:
  * Effective temperature: obtained from the Stefan–Boltzmann law
  * Excitation temperature: defined by the Boltzmann equation
  * Ionization temperature: defined by the Saha equation
  * Kinetic temperature: defined by the Maxwell–Boltzmann distribution
  * Color temperature: defined by fitting the Planck function on the spectrum
* Effective temperature is a global descriptor;
  rest apply locally, varying by gas location and other conditions
* All local temperature definitions result in the same value
  at thermodynamic equilibrium
* Local thermodynamic equilibrium (LTE) occurs when temperature change is
  very gradual compared to the mean free path length: H « L
  * Temperature scale height: H(T) ≡ T/|dT/dr|
  * Mean free path: l = 1/nσ,
    * Where n is atomic number density and σ ≡ π(2a₀)² is collision cross section

### Opacity

* Both scattering and pure absorption reduce the intensity of directed light:
  dI = −κρI ds
  * Where κ is a absorption coefficient (opacity) dependent on
    ray wavelength and gas composition, density, and temperature
* Characteristic distance l = 1/kρ is equivalent to main free path for photons
* Pure absorption decreases intensity exponentially: I = I₀ e^(−τ)
  * Where τ is optical depth with dτ = −κρ ds
* Optical depth is equivalent to the number of mean free path lengths along the
  path of a ray; gas is “optically thick” if τ » 1, “optically thin” if t « 1
* Balmer jump: opacity of stellar material suddenly increases for wavelengths
  below the ionization energy of H I in the first excited state
* Sources of opacity classified by initial and final quantum states of the
  interacting electrons:
  * Bound–bound transitions (excitation and de-excitation)
    * Small except for wavelengths corresponding to specific excitation energies;
    * If de-excitation returns to the initial state through emission,
      the photon is essentially scattered
    * If multiple photons are emitted from one absorption, the average photon
      energy is reduced
    * If de-excitation is instead triggered by collision, then the energy is
      converted into kinetic/thermal
  * Bound-free absorption (photoionization) and free–bound emission (recombination)
    * Wavelength-dependent cross section comparable to that of collisions
    * Adds to continuum opacity, as any photon above an electron’s
      ionization energy may be absorbed
    * Recombination emits one or more photons,
      again reducing average photon energy and scattering
  * Free–free absorption and bremsstrahlung (“braking radiation”) emission
    * Free electron gains speed from absorption, or loses speed from emission
    * Must take place next to an ion to conserve both energy and momentum
    * Adds to continuum opacity
  * Free electron (Thomson) scattering
    * Photon is scattered when the electron is made to oscillate in its EM field
    * Very tiny wavelength-independent cross section
    * Only dominates in stellar interiors and the hottest atmospheres,
      where near-total ionization eliminates bound-electron processes
  * Loosely-bound electron (Compton/Rayleigh) scattering
    * Compton if λ « a₀;
      very small change in scattered photon wavelength, much like Thomson
    * Rayleigh if λ » a₀;
      proportional to 1/λ⁴,
      so only significant in UV for supergiant stars’ extended envelopes,
      cool main-sequence stars, and planet atmospheres
* Electron scattering and photoionization of He are primary sources of continuum opacity
  for type O
* Photoionization of H and free–free absorption primary sources for types B–A
* Photoionization of H⁻ is primary source for type F0 and cooler
  * Very low binding energy of 0.754 eV (λ = 1640 nm) for bound–free opacity
  * Also contributes to free–free absorption at longer wavelengths
* Molecules survive in planetary and cooler stellar atmospheres
  * High opacity from large number of discrete absorption lines
  * May break apart during absorption through photodissociation
* Rosseland mean opacity: weighted harmonic mean of opacity across all wavelengths,
  accounting for rate of change of the blackbody spectrum with temperature;
  1/κ ≡ [ ∫ (1/κ_ν) (∂B_ν/∂T) dν ] / [ ∫ (∂B_ν/∂T) dν ]
  * Obeys Kramer’s opacity law: κ ∝ ρ/(T^3.5)
  * No analytic solution, but approximations exist for bound–free and free–free opacities
  * κ_bf ∝ (g_bf/t) Z(1+X) ρ/(T^3.5), k_ff ∝ g_ff (1−Z)(1+X) ρ/(T^3.5)
    * X, Y, and Z are mass fractions of H, He, and metals respectively
    * g_bf and g_ff are Gaunt factors correcting for quantum effects; ~1 for visible and UV
    * t is the guillotine factor, typically ranging 1–100

### Radiative transfer

* No net energy change occurs in any layer of a star that is in steady-state equilibrium
* Emission processes complement each of the primary absorption processes,
  resulting in randomly-directed scattering
* Specific intensity: dI = −κρI ds + jρ ds
   * Where κ is the absorption and j is the emission coefficient; both are λ-dependent
* Equation of radiative transfer: −(1/κρ) dI/ds = dI/dτ = I − S
  * Where source function: S ≡ j/κ
  * Expresses how light beam composition tends to resemble the local source of photons
* In local thermodynamic equilibrium, S = B, the blackbody Planck function
  * Integrating over all wavelengths, S = B = σT⁴/π
  * If τ » 1, I = B as well
* Radiation differential “driving” net radiative flux of photons flowing to the surface:
  dP/dτ = −F/c
* Radiative flux throughout star is constant: F = σT_eff⁴
* In a random walk with N steps of average size l, d = l√N
* On average, a photon from optical depth τ needs τ² steps to reach the surface
  * Photons from τ ≈ 1 may escape without scattering
* Eddington approximation: T⁴ ≈ ¾ T_eff⁴ (τ + ⅔)
  * Star has its effective temperature at optical depth τ ≈ ⅔
  * Average observed photon is emitted from τ ≈ ⅔, independent of angle
* Higher opacity corresponds to shorter pathlength for the same optical depth,
  thus absorption lines must come from outer, cooler layers
* Limb darkening: line of sight reaches τ = ⅔ in cooler, dimmer layers when
  observing closer to the edge of the disk

### Spectral line profile

* Core of line is formed at higher and cooler layer;
  formation descends down to continuum region at wings of the line
* Line is spectrally thin if radiant flux F(λ) is never 0
* Equivalent width: W = ∫ 1 − (F(λ)/F₀) dλ
  * Where the integrand is the depth of the line at wavelength λ
* Broadening processes:
  * Natural: uncertainty principle means momentary occupancy of
    an excited state for small Δt amplifies uncertainty of orbital energy
    ΔE ≈ ħ/Δt, resulting in uncertainty of wavelength
    Δλ ≈ (λ²/2πc) [(1/Δt₀) + (1/Δt₁)], where Δt are lifetimes in the two states
  * Doppler: random thermodynamic motion results in nonrelativistic Doppler shift
    Δλ ≈ (2λ/c) √(2kT/m)
    * in giant and supergiant stars, random large-scale turbulence increases to

      Δλ ≈ (2λ/c) √(2kT/m + v²), where v is the most probable turbulence speed
    * coherent mass motions such as rotation, pulsation and mass loss also
      substantial factors
  * Pressure: collisions with neutral atoms and pressure from nearby
    ion electrical fields may perturb the orbitals with
    Δλ = λ²/πcΔt₀, where Δt₀ ≈ 1/[ nσ√(2kT/m) ] is the average collision time
    * Dependence on atomic number density explains narrower lines of sparser
      giant and supergiant extended atmospheres used for the
      Morgan–Keenan classification
* Damping (Lorentz) profile: shape of lines from natural and pressure broadening,
  characteristic of radiation from electric charge in damped harmonic motion
* Voigt profile: total line profile from both Doppler and damping profiles,
  with Doppler dominating at the core and damping dominating at wings
* Schuster–Schwarzchild model: assumes photosphere is a blackbody
  and atoms above it create absorption lines
* Number of atoms involved in absorption per surface area: fN
  * Where oscillator strength f is the relative likelihood of each transition,
    and column density N is the area density of absorbing atoms in a column from the surface to the observer
* Curve of growth for line width as a function of column density:
  * Initially optically thin core: W ∝ N
  * Then, saturated core with optically thin wings: W ∝ √(ln N)
  * For high N, pressure-broadening profile dominates wings: W ∝ √N
* Applying Boltzmann and Saha equations to curve of growth finds total number of atoms above continuum layer
  * Applying Boltzmann equation also finds excitation temperature
  * Applying Saha equation also finds either e⁻ pressure or ionization temperature from the other

### Computer modeling

* Construction of an atmospheric model with extensive physics fine-tuned to observations
  can provide information on line profile, chemical composition, effective temperature,
  surface gravity, etc

## 10. The Interiors of Stars

* Direct observation only possible with neutrinos and the occasional supernovae
* Study requires physically accurate computer models that match observable surface features

### Hydrostatic equilibrium

* Equilibrium condition requires pressure gradient to counteract gravity: dP/dr = −ρg
  * Local gravity: g ≡ GMρ/r², where M is mass enclosed within r
* Mass conservation equation: dM/dr = 4πr²ρ
* Total pressure: P = ρkT/μm_H + aT⁴/3
  * Where first term is ideal gas law and second term is radiation pressure
  * Mean relative molecular weight: μ ≡ m_avg/m_H

### Kelvin–Helmholtz mechanism

* Contraction releases gravitational energy: ΔE ≈ (3/10) GM²/R
* Gravitational energy contribution: ϵ = −dQ/dt = −T dS/dt
  * Where specific entropy: dS ≡ dQ/T
  * Contraction produces heat and decreases entropy;
    vise versa for expansion
* Kelvin–Helmholtz timescale of star: t = ΔE/L ≈ 10⁷ years for Sun

### Nuclear fusion

* Releases difference in strong nuclear force binding energy:
  * E = Δmc² = [ Z m_p + (A−Z) m_n − m_nucleus ] c²
* Timescale of ~10¹⁰ years for sun
* Must overcome Coulombic repulsion between + nucleons: ~10⁶ eV at r ≈ 10⁻¹⁵ m
* Thermal energy of gas: E = μv²/2 = p²/2μ
  * Where μ is reduced mass and v is average relative velocity
* Requires impossibly high temperatures in classical physics:
  E = (3/2) kT = q₁q₂/4πrϵ₀ → T ≈ 10¹⁰ K for H–H fusion
* With quantum tunneling (de Broglie λ = h/p ≈ r):
  E = (h/λ)²/2μ = q₁q₂/4πλϵ₀ → T ≈ 10⁷ K, consistent with solar core temperature estimates
* Reaction rate depends on velocity (wavelength) distribution, densities, and
  cross-sectional area of particles 1 and 2: r = ∫ n₁ n₂ σ(E) v(E) (n(E)/n) dE
  * Power-law approximation: r = r₀X₁X₂(ρ^α´)(T^β),
    * Where X are mass fractions, α´ ≈ 2, and β ranging from ~1 to ≥40
* Velocity distribution follows Maxwell–Boltzmann
* Cross section: σ(E) = (S/E) e^(−b/√E)
  * Where de Broglie area is πλ² ∝ 1/p² ∝ 1/E,
    and tunneling across barrier of energy U is e^(‐2π²U/E) ∝ e^(−b/√E),
    with b ∝ q₁q₂√μ
  * S may be slow-varying function of E, or may have sharp peaks from resonance of
    specific energy levels within the nucleus
* Electron screening from sea of free e⁻ partially hides nuclei,
  reducing effective charge and Coulomb barrier,
  sometimes enhancing He production by 10–50%
* Likelihood of nuclear reaction: e^(−b/√E) e^(−E/kT)
  * Product of velocity distribution’s high-energy tail and the quantum tunneling terms
  * Most likely energy is named Gamow peak: E₀ = (bkT/2)^(2/3)
* Luminosity gradient equation: dL/dr = 4πr²ρ,
  * Where L is luminosity enclosed in r,
    and ϵ is specific power (power released per mass)

### Nucleosynthesis

* _Notation: ⁴₂He, where 4 is mass number and 2 is proton number_
* Simultaneous 4-body collision 4 ¹₁H + ? → ⁴₂He + ? extremely unlikely
* Reaction chain of 2-body interactions more probable
* Interactions must obey conservation laws:
  electric charge, nucleon number, and lepton number
* Proton–proton chain (H-burning): 4 ¹₁H → ⁴₂He + 2e⁺ + 2ν_(e⁻) + 2γ
  * Three branches (see Fig. 10.8)
  * In Sun, 69% PP I, 31% PP II, 0.3% PP III
  * Near T = 1.5 × 10⁷ K: ϵ ≈ ϵ₀ ρ (X_H)² (T/10⁶)⁴
    * Where ϵ₀ = 1.08 × 10⁻¹² W·m³/kg²
* CNO cycle (H-burning): C, N, and O are catalysts
  * Two main branches, with the second only occuring ~0.04% of the time
  * Near T = 1.5 × 10⁷ K: ϵ ≈ ϵ₀ ρ (X_H) (X_CNO) (T/10⁶)²⁰
    * Where ϵ₀ = 8.24 × 10⁻³¹ W·m³/kg²
    * Much more temperature-dependent, only dominating in more massive stars
* Triple-α process (He burning): 2 ⁴₂He ↔ ⁸₄Be, ⁸₄Be + ⁴₂He → ¹²₆C + γ
  * First step produces an extremely unstable ⁸₄Be,
    thus combined is essentially a 3-body interaction: r ∝ (ρY)³
  * ϵ ≈ ϵ₀ ρ² Y³ (T/10⁸)⁴¹
    * Incredibly strong temperature dependence
* α process (C & O burning): ¹²₆C + ⁴₂He → ¹⁶₈O + γ, ¹⁶₈O + ⁴₂He → ²⁰₁₀Ne + γ
  * α capture becomes prohibitive at higher Z due to higher Coulomb barrier
  * C–C burning near 6 × 10⁸ K and O–O burning near 10⁹K can produce
    Na, Mg, Si, P, and S; some are endothermic but are normally less likely
* Binding energy per nucleon: E/A
  * Relative to atomic number A, local maxima are very stable
  * Magic nuclei: some elements (⁴₂He, ¹⁶₈O) have unusually high E/A
  * Broad peak around ⁵⁶₂₆Fe, the most stable nuclei

### Energy transport and thermodynamics

* Three mechanisms: radiation of photons (affected by opacity),
  convection, and conduction (generally insignificant)
* Radiative temperature gradient: dT/dr = −(3/4ac) (κρ/T³) L/4πr²
* If temperature gradient becomes too steep, convection takes hold
* Convection is much more complicated than radiation
  * Strongly coupled to the star’s dynamic behavior
  * 3D Navier–Stokes with turbulence is hard compute
  * Pressure scale height, convection’s characteristic length scale,
    is big: H ≡ −P dr/dP ≈ P/ρg ≈ R*/10
* First law of thermodynamics: dU = dQ − dW = dQ − P dV
  * U is a state function; Q and W are not — dQ and dW are inexact differentials
* Specific heat capacity: C ≡ ∂Q/∂T, C_P = C_V + nR
  * C_P is at constant pressure; C_V is at constant volume
  * Heat capacity ratio (adiabatic index): γ ≡ C_P / C_V
  * γ = 5/3 for a monoatomic gas;
    approaches 1 in a partial ionization zone as both specific heats increase
* Isochoric process (dV = 0): dU = dQ = C_V dT
* Adiabatic process (dQ = 0): dU = −P dV
  * Gas law: P ∝ 1 / V^γ ∝ ρ^γ ∝ T^[γ/(γ−1)]
  * Speed of sound: v = √(B/ρ) = √(γP/ρ)
    * Where bulk modulus: B ≡ −V ∂P/∂V | dQ = 0
* Adiabatic temperature gradient:
  dT/dr = −[1 − (1/γ)] (μ m_H / k) GM/r² = −g / C_P
* If the surrounding temperature gradient is steeper than the bubble’s,
  even just slightly, the condition becomes superadiabatic,
  and nearly all luminosity is transferred outwards adiabatically,
  via convection instead of radiation
* In general, convection occurs if a region
  1. has high opacity (surrounding dT/dr ∝ κ),
  2. is ionizing and raising the specific heat capacity (bubble dT/dr ∝ 1 / C_P), or
  3. has a highly temperature-dependent fusion process
* First two conditions can occur simultaneously;
  third only occurs deep in the interior with the CNO or triple-α processes
* Convective flux under mixing-length theory:
  F = ρ C_P (k / μ m_H) √[(T/g) δ(dT/dr)]³ √β α²
  * Where 0.5 ≲ α ≲ 3 and 0 ≲ β ≲ 1 are free parameters
    * α ≡ l/H, the ratio of the mixing length
      (distance traveled by bubble before thermaliing with surrounding)
      and the pressure scale height
    * β | βv² is the average kinetic energy of the bubble as it travels over l

### Stellar model building

* Basic stellar models need constructive relations for P, κ, and ϵ:
  expressing them in terms of density, temperature, and composition
  * P (pressure) generally modelable with ideal gas law and radiation pressure,
    but is more complex in certain stars’ deep interiors
  * κ (mean opacity) calculated explicitly from tabular data or fitting functions
  * ϵ (specific power) calculated analytically from reaction networks
* Boundary conditions:
  * As r → 0, enclosed M and L both → 0
  * And ignoring extended atmospheres and mass loss:
    as r → R*, T, P, and ρ all → 0
* Vogt–Russell theorem:
  due to a star’s dependence on nuclear burning,
  its mass and internal composition uniquely determine
  its radius, luminosity, internal structure, and subsequent evolution
  * Ignores smaller influences such as magnetic fields and rotation
* General modeling numeric integrates shell-by-shell
  with the system of differential equations,
  often from from a transition point towards both the surface and the center
* Polytropes: simplified stellar models in which P(ρ) ∝ ρ^γ
* Lane–Emden equation: (1/ξ²) (d/dξ)[ξ² d(D_n)/dξ] = −(D_n)^n
  * Where the polytropic index: n | γ ≡ (n+1)/n
  * TODO: too complicated

### Main sequence

* Vast majority of stars have H mass fraction X ≈ 0.7 and metal mass fraction 0 ≲ Z ≲ 0.03
* Changes in core composition affects observed surface features
* Very light stars (M ≲ 0.08 M☉) are not hot enough to let fusion
  stabilize against gravitational contraction
  * Highly opaque and fully convective
* Lower-initiation-energy p–p chain dominates for low-mass stars
  * Shallow core thermal gradient leads to radiative core
  * High shell opacity leads to convective shell
* Highly-temperature-dependent CNO cycle dominates for high-mass stars (M ≳ 1.2 M☉)
  * Steep core thermal gradient leads to convective core
  * Low shell opacity leads to radiative shell
* Very massive stars (M ≳ 90 M☉) have rapid core thermal oscillations
  affecting fusion rates
* Very massive stars may also have radiation pressure exceed gas pressure at
  outer layers, with maximum stable luminosity given by Eddington limit:
  L ≤ 4πGcM/κ
* Main-sequence lifetimes decrease with increasing luminosity

## 11. The Sun

*

## 12. The Interstellar Medium and Star Formation

## 13. Main Sequence and Post-Main-Sequence Stellar Evolution

## 14. Stellar Pulsation

## 15. The Fate of Massive Stars

## 16. The Degenerate Remnants of Stars

## 17. General Relativity and Black Holes

## 18. Closed Binary Star Systems

# III. The Solar System

# IV. Galaxies and the Universe

To be continued...
