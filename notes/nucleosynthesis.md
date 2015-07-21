## Nucleosynthesis

Big Bang Nucleosynthesis is a process of generation of the light nuclei in the Early Universe
from the primordial protons and neutrons. Proton-neutron composition of the plasma and the
baryon-to-photon ratio are the main input parameters for the system of nuclear kinetic equations of
the BBN.

<center><img src="images/bbn_reactions.png" /></center>

### Deuterium bottleneck

Generation of the nuclei is controlled by the balance of the deuterium generation reaction.

\begin{equation}
n + p \leftrightarrow D + \gamma
\end{equation}

The bounding energy of the deuterium nucleus is $\Delta_D \sim 2.2 MeV$ - the minimal
energy requirement for deuterium generation. However, due to extremely small baryon-to-photon
ratio, the are plenty high-energetic photons that drive the inverse reaction.

The number of photons with $p > \Delta_D$ is given by:

\begin{align}
n_\gamma (E>\Delta_D) = \frac{1}{\pi^2} \int_{\Delta_D}^\infty \frac{p^2 dp}{e^\frac{p}{T}-1}
= \frac{T^3}{\pi^2} \int_{\Delta_D/T}^\infty \frac{x^2 dx}{e^x-1}
\end{align}

As the first approximation, for efficient deuterium generation to begin, this value has to be
comparable with $n_B$

\begin{equation}
n_\gamma (E>\Delta_D) \sim n_B = n_\gamma \eta \approx 0.244 T^3 \cdot \eta
\end{equation}

### Backreaction

Nucleosynthesis is a subleading process in the Early Universe as it governs mostly non-relativistic
particles in the radiation dominated epoch. Due to this one can neglect any backreaction of the
nuclear reactions on the expansion of the Universe and plasma.

### Computation

As all baryons are non-relativistic, their distribution functions are significantly spiked at the
zero momentum and it is sufficient to treat them using simplified Boltzmann equations as in
the case of sterile neutrinos. Originally nucleosynthesis was considered by Gamow as a system of
neutron capture/emission processes ([ $\alpha\beta\gamma$ paper](http://journals.aps.org/pr/pdf/10.1103/PhysRev.73.803)):

\begin{equation}
\frac{d n_i}{d t} = f(t)( \sigma_{i-1} n_{i-1} - \sigma_{i} n_i)
\end{equation}

($f(t)$ - expansion factor, $\sigma_i$ - neutron capture cross-sections).

More precisely, one can consider a coupled system of ordinary differential equations expressing not
only neutron capture processes, but all possible reactions between the nucleons and light nuclei.
In conditions of primordial nucleosynthesis it is sufficient to include only 2-particle scatterings
and nuclei up to $^7 Li, ^7 Be$. Heavier elements do not generate during the BBN because of absence
of stable nuclei with $A=5,8$ and more complicated reactions like $^4 He + ^4 He + ^4 He \to ^{12} C$
are negligible at the densities of the interest.
