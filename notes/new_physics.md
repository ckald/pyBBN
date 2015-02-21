= Sterile neutrinos in the mass range $m_\pi \div m_K$ =

Interactions of sterile neutrinos mostly mimic the weak interactions of active neutrinos. However, \
non-zero mass of $N$ allows many channels previously forbidden by kinematics \
(e.g., $\nu \to \nu \nu \nu$ is forbidden, but $N \to \nu \nu \nu$ isn't).

In general, listing all leading order interactions of some particle requires building perturbation \
series for the corresponding Hamiltonian. This is a rather complicated task in case of the Standard\
Model, so it is convenient to use effective field theory ideas.

To account efficiently for all new interactions of sterile neutrinos with $m_N \ge m_\pi$, \
we need to list all possible effective interaction terms in the Lagrangian and take into \
consideration the thermodynamical state of the Universe at the moment of interest.

== Primordial plasma ==

We consider the Early Universe at temperatures around few $MeV$. Plasma mostly consists of photons,\
active neutrinos, electrons and positrons - and tiny additions of long-living heavier particles \
like muons, protons and neutrons.

Thermodynamical state of the system is characterized by equilibrium in the electromagnetic sector\
(maintained by high-rate EM-interactions), while weakly-only interacting particles begin to \
decouple. This means that all weak processes rate are small comparing to the Hubble rate and it is\
possible to consider only few lowest orders in the weak expansion.

Additionally, we are interested only in the influence of interactions on non-equilibrium species as\
fast EM interactions dominate. Here is the tentative hierarchy of interaction rates:

\begin{equation}
    \Gamma_{EM} \gg H \ge \Gamma_{Weak} > \Gamma_{Sterile}
\end{equation}

== Processes of interest ==

For our discussion it is sufficient to consider only processes where only a single sterile neutrino\
is present. This narrows our list to processes of 2 kinds:

 * $N \to 1 + 2 + \dotsc$
 * $N + 2 \to 3 + 4 + \dotsc$

=== Kinematics ===

In case of sterile neutrino decay, allowed processes are limited solely the the masses of the final\
states:

\begin{align}
    p_N^\mu &= p_1^\mu + p_2^\mu + \dotsc
    \\\\ p_N^2 = m_N^2 &= \sum_i m_i^2 + \sum_{ij} (p_i \cdot p_j)
\end{align}

Meaning that

\begin{equation}
    m_N^2 \ge \sum_i m_i^2
\end{equation}

Thus let's list all particles with masses below the Kaon mass:

| Particle                    | Mass          |
|-----------                  |-----------    |
| $\gamma$                    | 0             |
| $\nu_\alpha$                | 0             |
| $e^\pm$                     | 0.511 MeV     |
| $\mu^\pm$                   | 105.7 MeV     |
| $\pi^0$ (uu, dd)            | 135.0 MeV     |
| $\pi^\pm$ (ud)              | 139.6 MeV     |
| $K^\pm$ (us)                | 493.7 MeV     |
| $K^0, \overline{K}^0$ (ds)  | 497.6 MeV     |

The only states of interest are three pi-mesons composed of weakly interacting quark pairs. \
As they are relatively short-lived, abundant generation of them is improbable. Then, beside the \
previously considered processes, we need to anticipate the appearance of pions in the\
final states.

This leads to the following processes:

 *  Pionic decays:
    \begin{align}
        N &\to \nu_\alpha + \pi^0
        \\\\ N &\to l^- + \pi^+
    \end{align}
 *  Leptonic scatterings:
    \begin{align}
        N + l^+ \to \pi^+
    \end{align}
 *  Exotic:
    \begin{align}
        N + \mu^+ \to K^+
    \end{align}

Apparently, kaon reactions can appear already with $m_N > 384 MeV$, but we will see that\
interactions with kaons are suppressed by $|V_{us}|$ and the rareness of muons.

The backreaction of the pions on the plasma is dominated by few decay modes:

\begin{align}
    \pi^0 &\to \gamma + \gamma
    \\\\ \pi^+ &\to \mu^+ + \nu_\mu
    \\\\ \pi^+ &\to e^+ + \nu_e
\end{align}

== Pion interactions ==

According to the quark model, mesons are bound states formed by quark pairs and gluons. In practice\
this means that quark quantum states inside mesons are deformed and one cannot directly apply the\
Feynman calculus. However, in case of electroweak interactions, only the quarks carry the \
corresponding quantum numbers and the influence of gluons can be omitted.

=== Phenomenological model ===

According to experimental evidence, pi-mesons constitute a pseudoscalar triplet. Simplest possible\
interactions for a pion with fermionic current then have to be of the form

\begin{align}
    g J_\mu^{\pm(0)} F^\mu + g' \overline{f} \gamma^5 \pi f
\end{align}

The latter interaction is not realized in nature (proven by experiments and can be explained from\
theoretical point of view by Goldstone boson-like origin of the pions).

$F^\mu$ is referred to as the "pion form-factor". It's precise form can be reconstructed by \
recalling that for a spin-0 particle the only associated 4-vector is it's momentum. Moreover,\
multiplier of this vector has to be a Lorentz scalar, thus depending only on the $p^2 = m_\pi^2$.\
However, this multiplier can differ for neutral and charged pions.

\begin{equation}
    F^\mu = f_\pm p^\mu
\end{equation}

=== Gauge bosons interactions ===

According to the phenomenological model, pi-mesons is gradiently coupled to SM currents, thus we\
assume that this interaction has a form analogous to the regular current interactions. In \
particular, we can apply Feynman rules for quarks and fix the renormalization due to gluons using\
the observed lifetime of the meson.

<center><img src="pion_decay.png" width=200 /></center>

=== Feynman rules ===

\begin{align}
    \Delta \mathcal{L}\_{CC} = \frac{g}{2 \sqrt{2}}
    \begin{bmatrix} \overline{u} & \overline{c} & \overline{t} \end{bmatrix}
    \gamma^\lambda (1-\gamma^5) V_{CKM}
    \begin{bmatrix} d \\\\ s \\\\ b \end{bmatrix}
    W_\lambda^+ + h.c.
\end{align}

\begin{align}
    \Delta \mathcal{L}\_{NC} = \frac{g}{\cos \theta_W}
    \left(
          \epsilon_L^{(f)} \overline{f}_L \gamma^\lambda f_L
          + \epsilon_R^{(f)} \overline{f}_R \gamma^\lambda f_R
    \right) Z\_\lambda
\end{align}

where

\begin{align}
    \epsilon_L^{(f)} &= -\frac12 - Q_f \sin^2 \theta_W  && \text{for}  &f = e,\mu,\tau,d,s,b
    \\\\ \epsilon_L^{(f)} &= +\frac12 - Q_f \sin^2 \theta_W
    && \text{for}  &f = \nu_e,\nu_\mu,\nu_\tau,u,c,t
    \\\\ \epsilon_R^{(f)} &= - Q_f \sin^2 \theta_W
\end{align}

or

\begin{align}
    \Delta \mathcal{L}\_{NC} = \frac{g}{\cos \theta_W}
    \overline{f} (v_f - a_f \gamma^5) f \; Z\_\lambda
\end{align}

\begin{align}
    v_f = \pm_f \frac14 - Q_f \sin^2 \theta_W && a_f = \pm_f \frac14
\end{align}

\begin{align}
\begin{cases}
    \pm_f = -1 && \text{for} &f = e,\mu,\tau,d,s,b
    \\\\ \pm_f = +1 && \text{for} &f = \nu_e,\nu_\mu,\nu_\tau,u,c,t
\end{cases}
\end{align}

Direct interaction of the pion states with bosons is considered to be

\begin{align}
    \Delta \mathcal{L}\_{\pi W} &= \frac{g}{2\sqrt{2}} f_\pm |V_{ud}| (\partial_\mu \pi) {W^\dagger}^\mu + h.c.
    \\\\ \Delta \mathcal{L}\_{\pi Z} &= \frac{g}{2 \cos \theta_W} f_0 (\partial_\mu \pi^0) Z^\mu
\end{align}

Turns out that $f_\pm = \sqrt{2} f_0 = f_{\pi} \approx 130 MeV$

Relation between $f_\pm$ and $f_0$ can be related to the fact that precise quark composition of the\
$\pi^0$ is $\frac{1}{\sqrt{2}}(u \overline{u} + d \overline{d})$
