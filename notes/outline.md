= Project outline =

== Sterile neutrinos generation and decoupling ==

Heavy sterile neutrinos ($N_{1,2}$) are produced thermally at the temperatures $\gtrsim 10^2 MeV$ \
(vaguely, depending on the mixing angle and the mass). Influence on the universe expansion at the\
later stages depends on the regime in which sterile neutrinos decouple from the plasma.

Decoupling happens when the interaction rate of the equilibrium-supporting reactions drops below\
the Hubble rate. In the radiation dominated epoch, Hubble rate can be conveniently approximated:

\begin{equation}
    H \approx \frac{T^2}{M_*}
\end{equation}

Typical interaction rate for 2-to-2 scattering ($N + \nu \to N + \nu$) of the sterile \
neutrino with the coupling constant $g$ is of the form

\begin{equation}
    \Gamma_N = \left< \sigma n_\nu v\right>
\end{equation}

where $\sigma$ is the cross-section of the reaction $\propto g^2$, and $v$ is the relative velocity.

Coupling constant for sterile neutrinos is $G_F |\theta|$, so cross-section from dimensional \
considerations looks like

\begin{equation}
    \sigma \sim [E]^{-2} = G_F^2 |\theta|^2 [E]^2
\end{equation}

Kinematically considered reaction contains few parameters: invariant mass and scattering angles\
The total cross-section does not depend on the angle, so $\sigma = \sigma(s)$.

\begin{equation}
    s \sim
    \begin{cases}
        m^2 & \text{for non-relativistic N} \\\\
        T^2 & \text{for ultra-relativistic}
    \end{cases}
\end{equation}

Then

\begin{equation}
    \Gamma_N \sim G_F^2 |\theta|^2 T^3
    \begin{cases}
        m^2 \\\\ T^2
    \end{cases}
    \sim \frac{T^2}{M_*}
\end{equation}

This gives us corresponding decoupling temperatures for sterile neutrinos:

\begin{equation}
    T^{-1} \sim
    \begin{cases}
        G_F^2 |\theta|^2 M_\* m^2
        \\\\ (G_F^2 |\theta|^2 M_\*)^\frac13
    \end{cases}
\end{equation}

For $|\theta|^2 \approx 10^{-8}$ and $m = 0.3 GeV$

\begin{equation}
    T^{-1} \sim
    \begin{cases}
        (300 GeV)^{-4} \cdot 10^{-8} \cdot 2.4 \cdot 10^{18} GeV \cdot (0.3 GeV)^2
        \\\\ \left((300 GeV)^{-4} \cdot 10^{-8} \cdot 2.4 \cdot 10^{18} GeV \right)^\frac13
    \end{cases}
    \sim
     \begin{cases}
        240 GeV^{-3} \cdot 10^2 \cdot 9 \cdot 10^{-2} GeV^2
        \\\\ (240 GeV^{-3})^\frac13
    \end{cases}
    \sim
     \begin{cases}
        22 GeV^{-1}
        \\\\ 6 GeV^{-1}
    \end{cases}
\end{equation}

\begin{equation}
    T \sim
     \begin{cases}
        45 MeV
        \\\\ 166 MeV
    \end{cases}
\end{equation}

Hence, the distribution function of sterile neutrinos at the later times highly depends on the\
following factors:

  * sterile-active mixing angle temperature dependence
  * mass of the sterile neutrino
  * regime of decoupling (relativistic/non-relativistic)

For example, if thermal corrections of $\theta$ are significant at the time of relativistic \
decoupling "in vacuum", sterile neutrinos might stay in the equilibrium longer and decouple being \
non-relativistic.

=== Mixing angle temperature dependence ===

A simple example of influence of the medium on sterile neutrino mixing (as well as active neutrinos\
oscillations) is the Miheev-Smirnov-Wolfenstein (MSW) effect. It is discussed in the theory of \
solar neutrinos oscillations and captures the influence of the solar $e^-$-$e^+$ plasma on the \
propagation of neutrinos. Basically, the reason for any interaction between sterile and active \
neutrinos and oscillation is the fact that mass operator of the Hamiltonian does not commute with\
the interaction operator. Hence, they don't have a shared set of eigenfunctions.

The easiest way to account for the thermal corrections to the mixing angle is the effective field\
theory. Let's consider the basic charged current scattering diagram of the active neutrinos on \
the electrons:

\begin{equation}
    \nu_e + e^+ \to \nu_e + e^+
\end{equation}

The corresponding term in the Lagrangian looks like:

\begin{equation}
    \Delta \mathcal{L} \sim \overline{\nu_e}(x) \gamma_\mu (1-\gamma^5) e(x)
        D_W^{\mu \nu}(x-y) \overline{e}(y) \gamma_\nu (1-\gamma^5) \nu_e(y)
\end{equation}

In the lowest order (i.e., Fermi theory) the propagator is simply:

\begin{equation}
    D_W^{\mu \nu}(x-y) \propto G_F g^{\mu \nu} \delta^4(x-y)
\end{equation}

If we average out the electron degrees of freedom, we will obtain the \
following effective Lagrangian term:

\begin{equation}
    <\Delta \mathcal{L}> \sim G_F \overline{\nu_e}(x) \gamma_\mu (1-\gamma^5)
        < e(x) \overline{e}(x) > \gamma_\nu (1-\gamma^5) \nu_e(x)
\end{equation}

Recalling the Dirac current:

\begin{align}
    < \overline{e}(x) \gamma^\mu e(x) > &= j^\mu(x)  \\\\
    < e(x) \overline{e}(x)> &= -j^0 (x) \gamma^0 = -(n_e(x)-n_\overline{e}(x)) \gamma^0
\end{align}

TODO: can't quite justify this part

Notice the difference from the MSW-effect, where the medium is absolutely anisotropic and one can\
neglect the density of the positrons.

Finally this gives us:

\begin{align}
    <\Delta \mathcal{L}>
        &\sim -G_F (n_e(x)-n_\overline{e}(x))
            \overline{\nu_e}(x) \gamma_\mu (1-\gamma^5) \gamma^0 \gamma_\nu (1-\gamma^5) \nu_e(x)
        \\\\ &= 8 G_F (n_e(x)-n_\overline{e}(x)) \, \overline{\nu_e}(x) \gamma^0 \nu_e(x)
\end{align}

Expanding this term in terms of chiral spinors, we see that it is basically an effective mass of
the active neutrino:

\begin{align}
    <\Delta \mathcal{L}>
        &= 8 G_F \Delta n_e \, {\nu_e}_L^+ {\nu_e}_L
\end{align}

From this we immediately see that this effect automatically cancels out in the lepton-symmetric\
medium.

Using the same technique, we can find a correction to $\theta$:

\begin{align}
    <\Delta \mathcal{L_\theta}> &\sim 8 \theta_0 G_F \Delta n_e \, \overline{\nu_e} \gamma^0 N
\end{align}

Where $\theta_0$ is the vacuum value of the mixing:

\begin{equation}
    \theta(T) = \frac{M_D(T)}{M_I} = \theta_0 \left(1 + \frac{G_F \Delta n_e}{M_I} + ...\right)
\end{equation}

This correction is insignificant when

\begin{equation}
    M_I \gg G_F \Delta n_e \approx G_F T^3 \eta_l
\end{equation}

Depending on the assumptions in the system, this correction might cancel out ($\eta_l = 0$). \
The maximal influence of this correction corresponds to the MSW situation: $\eta_l = 1$:

\begin{align}
    T &\ll \left(\frac{M_I}{G_F \eta_l}\right)^\frac13
        \sim \left(0.3 GeV (3 \cdot 10^2 GeV)^2 \right)^\frac13
        \\\\ &= \left(27 \cdot 10^3 \right)^\frac13 GeV = 30 GeV
\end{align}

So, as long as the temperature is below $1 GeV$ it is completely safe to neglect this effect.

However, due to the symmetry of the system, the subleading terms that give a rise to corrections\
proportional to $n_e + n_\overline{e}$ might become dominant.

=== Subleading effects ===

Returning to the interaction term of the Lagrangian:

\begin{equation}
    \Delta \mathcal{L} \sim \overline{\nu_e}(x) \gamma_\mu (1-\gamma^5) e(x)
        D_W^{\mu \nu}(x-y) \overline{e}(y) \gamma_\nu (1-\gamma^5) \nu_e(y)
\end{equation}

Expanding the W-boson propagator, we obtain new effective interactions of the kind:

\begin{align}
    \Delta \mathcal{L}
        & \sim G_F \overline{\nu_e} \gamma_\mu e \frac{\partial^2}{M_W^2} \overline{e} \gamma^\mu \nu_e +
        \\\\ &+ G_F \overline{\nu_e} \gamma_\mu e \frac{\partial^\mu \partial^\nu}{M_W^2} \overline{e} \gamma_\nu \nu_e
\end{align}

TODO: I don't see how to derive the $n+n$ effect in this framework. Looks like double-derivative\
    terms won't give the necessary relative sign between the electrons and positrons.

== Sterile neutrinos decay ==

After sterile neutrinos had decoupled from the plasma, their distribution as a function of \
conformal momentum stays constant. During this time sterile neutrinos do not scatter with the \
equilibrated plasma but they do decay.

Main decays channels for sterile neutrinos with masses $m_N < m_K$ range are:

\begin{align}
    N &\to \nu_\alpha + \nu_\beta + \overline{\nu_\beta} \\\\
    N &\to \nu_\alpha + l^+\_\beta + l^-\_\beta \\\\
    N &\to \nu_\alpha + l^+\_\alpha + l^-\_{\beta \neq \alpha} \\\\
    N &\to \pi^0 + \nu_\alpha \\\\
    N &\to \pi^+ + l^-_\alpha
\end{align}

Summed up decay rates of these processes gives us the total decay rate of the sterile neutrino. \
As the interaction rates of weak and electromagnetic processes are still sufficiently higher than\
the expansion rate of the Universe, all decay products of the sterile neutrino immediately\
equilibrate and affect only the temperature of the plasma.

Then the kinetic equation on the sterile neutrino reads:

\begin{equation}
    \frac{d n_N}{d t} = - (3 H + \Gamma_N) n_N
\end{equation}

Decay products function as a heating system of the plasma. This effect can be seen from the energy\
conservation law. For simplicity, let's consider a system consisting only of photons and \
non-relativistic sterile neutrinos:

\begin{align}
    &\dot{\rho} + 3 H ( \rho + p ) = 0
    \\\\ \rho_\gamma &= \frac{2 \pi^2}{30} T^4
        & p_\gamma &= \frac13 \rho_\gamma
        & \dot{\rho_\gamma} &= 4 \rho_\gamma \frac{\dot{T}}{T}
    \\\\ \rho_N &= m n_N
        & p_N &= 0
\end{align}

\begin{align}
    \dot{\rho} + 3 H (\rho + p) = \dot{\rho_\gamma} + \dot{\rho_N} + 3 H (\frac43 \rho_\gamma + m n_N) = 0
\end{align}

\begin{align}
    4 \rho_\gamma \frac{\dot{T}}{T} - (3 H + \Gamma_N) m n_N + 3 H (\frac43 \rho_\gamma + m n_N) &= 0
    \\\\ 4 \rho_\gamma \frac{\dot{T}}{T} - \Gamma_N m n_N + 4 H \rho_\gamma &= 0
\end{align}

\begin{equation}
    \frac{\dot{T}}{T} = \frac14 \Gamma_N \frac{\rho_N}{\rho_\gamma } - H
\end{equation}

This effect can be accounted for analytically without considering the precise dynamics of the \
decay products. From the equation follows that the amount of heating highly depends on the energy\
density ratio of the sterile neutrinos and photons, so non-relativistic particles won't affect the\
expansion rate significantly.

Depending on the decoupling temperature, one might also have to account for a transition in the \
plasma similar to the electron-positron annihilation at temperature $\sim m_e$ that changes the \
number of relativistic degrees of freedom - hence, the entropy conservation of the Universe will\
require another boost of temperature. Transition in this case is the QCD transition.

== Neutrinos decoupling and BBN setup ==


== Nucleosynthesis ==

Big Bang Nucleosynthesis is a process of generation of the light nuclei in the Early Universe\
from the primordial protons and neutrons. Proton-neutron composition of the plasma and the \
baryon-to-photon ratio are the main input parameters for the system of nuclear kinetic equations of\
the BBN.

<center><img src="images/bbn_reactions.png" /></center>

=== Deuterium bottleneck ===

Generation of the nuclei is controlled by the balance of the deuterium generation reaction.

\begin{equation}
n + p \leftrightarrow D + \gamma
\end{equation}

The bounding energy of the deuterium nucleus is $\Delta_D \sim 2.2 MeV$ - the minimal \
energy requirement for deuterium generation. However, due to extremely small baryon-to-photon\
ratio, the are plenty high-energetic photons that drive the inverse reaction.

The number of photons with $p > \Delta_D$ is given by:

\begin{align}
n_\gamma (E>\Delta_D) = \frac{1}{\pi^2} \int_{\Delta_D}^\infty \frac{p^2 dp}{e^\frac{p}{T}-1}
= \frac{T^3}{\pi^2} \int_{\Delta_D/T}^\infty \frac{x^2 dx}{e^x-1}
\end{align}

As the first approximation, for efficient deuterium generation to begin, this value has to be\
comparable with $n_B$

\begin{equation}
n_\gamma (E>\Delta_D) \sim n_B = n_\gamma \eta \approx 0.244 T^3 \cdot \eta
\end{equation}

=== Backreaction ===

Nucleosynthesis is a subleading process in the Early Universe as it governs mostly non-relativistic\
particles in the radiation dominated epoch. Due to this one can neglect any backreaction of the\
nuclear reactions on the expansion of the Universe and plasma.

=== Computation ===

As all baryons are non-relativistic, their distribution functions are significantly spiked at the\
zero momentum and it is sufficient to treat them using simplified Boltzmann equations as in \
the case of sterile neutrinos. Originally nucleosynthesis was considered by Gamow as a system of \
neutron capture/emission processes ([ $\alpha\beta\gamma$ paper](http://journals.aps.org/pr/pdf/10.1103/PhysRev.73.803)):

\begin{equation}
\frac{d n_i}{d t} = f(t)( \sigma_{i-1} n_{i-1} - \sigma_{i} n_i)
\end{equation}

($f(t)$ - expansion factor, $\sigma_i$ - neutron capture cross-sections).

More precisely, one can consider a coupled system of ordinary differential equations expressing not\
only neutron capture processes, but all possible reactions between the nucleons and light nuclei.\
In conditions of primordial nucleosynthesis it is sufficient to include only 2-particle scatterings\
and nuclei up to $^7 Li, ^7 Be$. Heavier elements do not generate during the BBN because of absence\
of stable nuclei with $A=5,8$ and more complicated reactions like $^4 He + ^4 He + ^4 He \to ^{12} C$\
are negligible at the densities of the interest.
