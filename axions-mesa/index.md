---
title:  "Astrophysical axions, and using MESA simulate the lives of stars"
---

*I first met Prof. Wick Haxton as a student in his quantum mechanics class in spring 2024.
After the course ended, 
he kindly took me on to help with one of his ongoing research projects.
Below I briefly summarize our present paper [[Ha25]](#ref-Ha25).
then review how I performed stellar simulations
with the Modular Experiments for Stellar Astrophysics (MESA) program,
generously open-sourced by Bill Paxton and his collaborators at KITP.
Coming into this project without any Fortran experience,
I found MESA rather hard to learn on my own;
I hope others who are tinkering with MESA for their research will find this helpful.*

{% include_relative 1.md %}
{% include_relative 2.md %}
{% include_relative 3.md %}

## Acknowledgements

Our simulations were made possible by the generosity of the MESA collaborators,
and we have in turn open-sourced our modifications at [github.com/xingyzt/saltyaxions](https://github.com/xingyzt/saltyaxions).

I would like to thank my collaborators Anupam Ray and Wick Haxton for their patient mentorship on my first research project.
I am supported by the N3AS undergraduate research program at UC Berkeley,
ran by the wonderful Laura Fantone, Rebecca Singh,
and Kenneth McElvain.
Many thanks as well to Pablo Castaño and Winston Yin, who first introduced me to axions in their inspiring talks.
On stellar evolution theory, feedback from Stan Woosley, Luca Boccioli, and Liang Dai were invaluable.
Finally, I am forever indebted to Bea Noether for teaching me the particle physics necessary to understand these axions’ theoretical appeal.

## References

1. <span id="ref-Ha25"> W. C. Haxton et al. </span> “A Continuous Galactic Line Source of Axions: The Remarkable
Case of 23Na”. Submitted to Phys. Rev. Lett. arXiv: [2505.03038 [astro-ph.HE]](https://arxiv.org/abs/2505.03038).

1. <span id="ref-Ha91"> W. C. Haxton and K. Y. Lee. </span> “Red giant evolution, metallicity and new bounds on
hadronic axions”. In: Phys. Rev. Lett. 66 (1991), pp. 2557–2560. doi: [10.1103/PhysRevLett.66.2557](https://doi.org/10.1103/PhysRevLett.66.2557).

1. <span id="ref-An09"> S. Andriamonje et al. </span> “Search for 14.4-keV solar axions emitted in the M1-transition of
Fe-57 nuclei with CAST”. In: JCAP 12 (2009), p. 002. doi: [10.1088/1475-7516/2009/12/002](https://doi.org/10.1088/1475-7516/2009/12/002). arXiv: [0906.4488 [hep-ex]](https://arxiv.org/abs/0906.4488).

1. <span id="ref-Ta13"> Koh Takahashi, Takashi Yoshida, and Hideyuki Umeda.</span>
“Evolution of progenitors for
electron capture supernovae”. In: Astrophys. J. 771 (2013), p. 28. doi: [10.1088/0004-637X/771/1/28](https://doi.org/10.1088/0004-637X/771/1/28). arXiv: [1302.6402 [astro-ph.SR]](https://arxiv.org/abs/1302.6402).

1. <span id="ref-La76"> S. A. Lamb, I. Iben Jr., and W. M. Howard.</span>
“On the evolution of massive stars through
the core carbon-burning phase.” In: 207 (July 1976), pp. 209–232. doi: [10.1086/154486](https://doi.org/10.1086/154486)

1. <span id="ref-Li24"> Marco Limongi et al. </span> “Evolution and Final Fate of Solar Metallicity Stars in the Mass
Range 7–15 M . I. The Transition from Asymptotic Giant Branch to Super-AGB Stars,
Electron Capture, and Core-collapse Supernova Progenitors”. In: Astrophys. J. Suppl.
270.2 (2024), p. 29. doi: [10.3847/1538-4365/ad12c1](https://doi.org/10.3847/1538-4365/ad12c1). arXiv: [2312.00107 [astro-ph.SR]](https://arxiv.org/abs/2312.00107).

1. <span id="ref-Zh19"> Shuai Zha et al.</span>
“Evolution of ONeMg Core in Super-AGB Stars toward Electron-capture
Supernovae: Effects of Updated Electron-capture Rate”. In: 886.1, 22 (Nov. 2019), p. 22.
doi: [10.3847/1538-4357/ab4b4b](https://doi.org/10.3847/1538-4357/ab4b4b). arXiv: [1907.04184 [astro-ph.HE]](https://arxiv.org/abs/1907.04184).

1. <span id="ref-Ca24"> Francesca Calore et al. </span> “Uncovering axionlike particles in supernova gamma-ray spec-
tra”. In: Phys. Rev. D 109.4 (2024), p. 043010. doi: [10.1103/PhysRevD.109.043010](https://doi.org/10.1103/PhysRevD.109.043010).
arXiv: [2306.03925 [astro-ph.HE]](https://arxiv.org/abs/2306.03925).

1. <span id="ref-Pr09"> Dina Prialnik.</span> *An Introduction to the Theory of Stellar Structure and Evolution.* 2nd ed.
Cambridge University Press, 2009. isbn: 9780521866040

1. <span id="ref-Sc17"> Josiah Schwab, Lars Bildsten, and Eliot Quataert. </span> “The importance of Urca-process
cooling in accreting ONe white dwarfs”. In: 472.3 (Dec. 2017), pp. 3390–3406. doi:
[10.1093/mnras/stx2169](https://doi.org/10.1093/mnras/stx2169). arXiv: [1708.07514 [astro-ph.SR]](https://arxiv.org/abs/1708.07514).
