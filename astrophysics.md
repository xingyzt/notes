# Carroll & Ostlie, Introduction to Modern Astrophysics: Xing's notes

## 1. The Celestial Sphere

### Altitude-Azimuth Coordinate System

* Local to observer location and day of year
* Altitude: h, measured from horizon towards zenith.
* Zenith: z = 90° - h
* Azimuth: A, measured from north on the horizon, clockwise, first towards
east on the horizon

### Equatorial Coordinate System

* Celestial equator: plane passing through Earth's equator
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
* J2000.0: Standard epoch Noon 2000-01-01 UT
  (universal time, measured from Greenwich)
* Uses Julian Calendar (365.25 d/yr), without Gregorian corrections

### Time standards

* Julian Date (JD): Days since Noon 4713BC-01-01
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
* Law of cosines for angles: cos A = - cos B cos C + sin B sin C cos a
* For a proper motion with position angle Φ (direction measured) and
  arclength Δθ,
  * Δα = Δθ sin Φ / cos δ
  * Δδ = Δθ cos Φ
  * (Δθ)² = (Δα cos δ)² + (Δδ)²

## 2. Celestial Mechanics

### Ellipses

* Semi-major axis: a
* Semi-minor axis: b
* Eccentricity: e | ae = d/2, where d = distance between foci
* b² = a² (1 - e²)
* Polar equation: r = a (1 - e²) / (1 + e cos θ)
* Area: A = πab

### Kepler's Laws

1. r = (L²/μ²)/[GM(1 + e cos θ)]
2. (dA/dt) = L/(2μ), constant
3. P² = ka³, k = 4π²/[G(M+m)], k = 1 using years and AU

### Orbital equations

* Virial theorem: 〈E〉 = 〈U〉/2, 〈U〉 = -2〈K〉
* Vis-viva: v² = G(M+m)[(2/r) - (1/a)]

## 3. Continuous Spectrum of Light

### Parallax and parsecs

* Parallax angle: p" (arcseconds),
  subtended by quarter orbit of Earth (1 AU)
* Distance in AU: d = 206,265 / p"
* Distance in parsecs (pc): d = 1 / p"

### Apparent magnitude

* Originally compiled by Hipparchus: 1 for brightest, 6 for dimmest
* Modern definition: F₂/F₁ = 100^[(m₁ - m₂)/5]
* Where radiant flux is density of emitted light: F = L/(4πr²)

### Absolute magnitude

* M, apparent magnitude of star if it's located at 10 pc
* H for planets, measured if located at 1 AU
* Distance modulus: d = 10^[(m - M + 5)/5],
m - M = 5 log₁₀(d) - 5
* For equally distant stars, ratio of radiant fluxes equals
  ratio of luminosities, thus: L₂/L₁ = 100^[(M₁ - M₂)/5],
* M_sun = +4.74, L_sun = 3.839E26 W

### Light

* c = λν
* Poyinting vector: S = (E × B)/μ₀, 〈S〉= E₀B₀/(2μ₀)
* Radiation pressure: P = (〈S〉cos θ) / c, absorption;
P = (2〈S〉A cos²θ)/c, reflection

### Blackbody Radiation

* Wien's displacement: λT = 0.002897755 (m·K), for peak wavelength
* Stefan-Boltzmann: L = 4πR²σT⁴,
  where T is effective temperature of the star's surface,
  and σ = 5.670400E-8 W/m²K⁴
* Planck's function: B(T) = (2hc²/λ⁵)/(e^(hc/λkT) - 1)

### Color Index

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
* U-B color index: U - B = M_U - M_B
* B-V color index: B - V = M_B - M_V, smaller is bluer
* Bolumetric correction: BC = m_bol - V = M_bol - M_V
* Color-color diagrams graph U-B against B-V,
  showing stars are non-ideal blackbodies

## 4. Special Relativity

* (Notation: prime (') means observed, unprimed means rest)

### Einstein's postulates:

1. Principle of Relativity: Laws of physics are the same in all inertial 
  reference frames
2. Constancy of the Speed of Light: In vacuum, c is independent of the motion
  of the light source

### Lorentz transformations:

* Lorentz factor: γ ≡ 1 / √(1 - u²/c²)
* Time dilation: t' = γt
* Length contraction: x' = x/γ
* Proper time (τ) and proper length measured at rest w.r.t. events, γ = 1

### Relativistic doppler shifts

* Doppler shift: ν' = ν / [ γ (1 + (u cos θ / c)) ],
  θ = 0° away from observer, 180° towards observer
* Redshift: z = (λ'/λ) - 1 = √[ (1 + (u cos θ / c)) / (1 - (u cos θ / c)) ]
* Time dilation from redshift: z + 1 = t'/t

### Relativistic momentum and energy

* Rest mass: m; Rest energy: mc²
* Relativistic momentum: p = γmv
* Relativistic kinetic energy: K = mc²(γ - 1)
* Relativistic total energy: E = γmc² = √( p²c² + m²c⁴ )
