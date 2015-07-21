## Sterile neutrinos generation and decoupling

### Decoupling temperature and regime

Heavy sterile neutrinos ($N_{1,2}$) are produced thermally at the temperatures $\gtrsim 10^2 MeV$
(vaguely, depending on the mixing angle and the mass). Influence on the universe expansion at the
later stages depends on the regime in which sterile neutrinos decouple from the plasma.

Decoupling happens when the interaction rate of the equilibrium-supporting reactions drops below
the Hubble rate. In the radiation dominated epoch, Hubble rate can be conveniently approximated:

$$
    H \approx \frac{T^2}{M_*}
$$

Typical interaction rate for 2-to-2 scattering ($N + \nu \to N + \nu$) of the sterile
neutrino with the coupling constant $g$ is of the form

$$
    \Gamma_N = \left< \sigma n_\nu v\right>
$$

where $\sigma$ is the cross-section of the reaction $\propto g^2$, and $v$ is the relative velocity.

Coupling constant for sterile neutrinos is $G_F |\theta|$, so cross-section from dimensional
considerations looks like

$$
    \sigma \sim [E]^{-2} = G_F^2 |\theta|^2 \cdot [E]^2
$$

Kinematically considered reaction contains few parameters: invariant mass and scattering angles
The total cross-section does not depend on the angle, so $\sigma = \sigma(s)$.

$$
    s \sim
    \begin{cases}
        m^2 & \text{for non-relativistic N}
        \\ T^2 & \text{for ultra-relativistic}
    \end{cases}
$$

Then

$$
    \Gamma_N \sim
    \begin{cases}
        G_F^2 |\theta|^2 T^3 m^2
        \\ G_F^2 |\theta|^2 T^5
    \end{cases}
$$

$$
    \Gamma_N \sim \frac{T^2}{M_*}
$$

This gives us corresponding decoupling temperatures for sterile neutrinos:

$$
    T^{-1} \sim
    \begin{cases}
        G_F^2 |\theta|^2 M_* m^2
        \\ (G_F^2 |\theta|^2 M_*)^\frac13
    \end{cases}
$$

Hence, the distribution function of sterile neutrinos at the later times highly depends on the
following factors:

  * sterile-active mixing angle temperature dependence
  * mass of the sterile neutrino
  * regime of decoupling (relativistic/non-relativistic)

For example, if thermal corrections of $\theta$ are significant at the time of relativistic
decoupling "in vacuum", sterile neutrinos might stay in the equilibrium longer and decouple being
non-relativistic.

Using the analysis above, we assume no significant thermal corrections for the mixing angles and,
for given mass and temperature, find the decoupling regime:

For $ \theta \approx 10^{-4} $
<img src="images/decoupling_regimes_10^-4.png" >
For $ \theta \approx 10^{-2} $
<img src="images/decoupling_regimes_10^-2.png" >

The plotted regions correspond to non-equilibrium species in different regimes. In practice,
as the temperature drops, when particle with given mass enters any of the regions - it decouples in
the corresponding regime. For the sakes of certainty, the "intermediate" regime corresponds to the
particle specie with $\frac15 T < m < 5 T$.

Analysis of the decoupling regimes for the mass range of interest shows that for
$\theta \lesssim 10^{-4}$ decoupling regime is rather relativistic, while for bigger values it is
possible to get a fully-non-relativistic decoupling.

### Decoupling and the QCD transition

Tuning of the parameters allows to shift the sterile neutrinos decoupling around the QCD transition
temperature $~200 MeV$. The strong sector constituents of the plasma change enormously, requiring
separate treatment of the quark-gluon and hadron phases and, possibly, the transition itself.

#### Decoupling above the QCD scale

At the temperatures above $\Lambda_{QCD} \sim 200 MeV$, plasma is full of quarks and gluons in
equilibrium. Sterile neutrinos are singlets on SM gauge groups, so the only important interactions
to consider are the weak interactions of quarks with active neutrinos. These processes are split
into neutral and charged channels and due to the conservation of the lepton number, each 2-to-2
scattering diagram contains only 2 quarks. One has to take into account the number of degrees of
freedom of quarks - in addition to the spin projections, each quark can carry one of 3 colors.
Color is conserved in weak interactions, so each interaction amplitude is simply multiplied by
$N_c = 3$.

As the first approximation, one can consider free quarks. The hadrons do not exist in the
quark-gluon plasma, so their interactions with the sterile neutrinos are not taken into account.

\begin{align}
    N + \overline{\nu}_\alpha &\to u_i + \overline{u}_i
    \\ N + \overline{\nu}_\alpha &\to d_i + \overline{d}_i
    \\ N + l^+_\alpha &\to u_i + \overline{d}_i
\end{align}

where $u_i$ and $d_i$ denote the up or down component of the $i$-th quark family.

##### Neutral channel reaction

\begin{align}
    N + \overline{\nu} &\to q + \overline{q}
\end{align}

\begin{align}
    \require{cancel}
    \imath \mathcal{M} = \left(\overline{\nu} \gamma^\mu \frac{-\imath g'}{2} P_L \frac{\imath M_D}{\cancel{N}} N \right)
        \frac{-\imath g_{\mu \nu}}{M_Z^2} (-\imath g') \left(\overline{q}_- \gamma^\nu (v_q - a_q \gamma_5) q_+ \right)
\end{align}

where the momenta of particles denoted by the particle symbols: $\nu, N, q_\pm$

\begin{align}
    \require{cancel}
    |\overline{\mathcal{M}}|^2 &= \left( \frac{g'}{M_Z} \right)^4 \frac{|\theta|^2}{2 m_N^2}
        Tr\left[ \cancel{N} (\cancel{N} \pm m_N) \cancel{N} P_R \gamma_\nu \cancel{\nu} \gamma^\mu P_L \right]
        Tr\left[ (\cancel{q_-} - m_q) \gamma_\mu (v_q -a_q \gamma_5)( \cancel{q_+} + m_q) (v_u + a_u \gamma_5) \gamma^\nu \right]
    \\ &= 8 \left( \frac{g'}{M_Z} \right)^4 |\theta|^2 \left(
        (a_q - v_q)^2 (N \cdot q_+) (\nu \cdot q_-)
        + (a_q + v_q)^2 (N \cdot q_-) (\nu \cdot q_+)
        - m_q^2 (a_q^2 - v_q^2) (N \cdot \nu)
    \right)
\end{align}

For the upper quarks:

\begin{align}
    |\overline{\mathcal{M_u}}(N \overline{\nu} \to q \overline{q})|^2 &=
        \frac29 \left( \frac{g'}{M_Z} \right)^4 |\theta|^2 \left(
        16 \sin^2{\theta_W} (N \cdot q_+) (\nu \cdot q_-)
        +(3-4 \sin^2{\theta_W})^2 (N \cdot q_-)(\nu \cdot q_+)
        -4 m_q^2 \sin^2{\theta_W} (-3+4 \sin^2{\theta_W}) (N \cdot \nu)
    \right)
\end{align}

For the lower quarks:

\begin{align}
    |\overline{\mathcal{M_d}}(N \overline{\nu} \to q \overline{q})|^2 &=
        \frac29 \left( \frac{g'}{M_Z} \right)^4 |\theta|^2 \left(
        4 \sin^2{\theta_W} (N \cdot q_+) (\nu \cdot q_-)
        +(3-2 \sin^2{\theta_W})^2 (N \cdot q_-)(\nu \cdot q_+)
        -2 m_q^2 \sin^2{\theta_W} (-3+2 \sin^2{\theta_W}) (N \cdot \nu)
    \right)
\end{align}

The charge-conjugated channels expressions can be obtained by simply swapping $q_\pm \to q_\mp$, so
the total interaction amplitudes are:

\begin{align}
    |\overline{\mathcal{M}}|^2 =
        |\overline{\mathcal{M}}(N \overline{\nu} \to q \overline{q})|^2
        + |\overline{\mathcal{M}}(N \nu \to \overline{q} q)|^2
\end{align}

$$
    |\overline{\mathcal{M_u}}|^2 = \frac29 \left( \frac{g'}{M_Z} \right)^4 |\theta|^2 \left(
        ((3-4 \sin^2{\theta_W})^2 + 16 \sin^2{\theta_W}) ((N \cdot q_+) (\nu \cdot q_-)+  (N \cdot q_-)(\nu \cdot q_+))
        + 8 m_q^2 \sin^2{\theta_W} (3-4 \sin^2{\theta_W}) (N \cdot \nu)
    \right)
$$

$$
    |\overline{\mathcal{M_d}}|^2 = \frac29 \left( \frac{g'}{M_Z} \right)^4 |\theta|^2 \left(
        ((3-2 \sin^2{\theta_W})^2 + 4 \sin^2{\theta_W}) ((N \cdot q_+) (\nu \cdot q_-)+  (N \cdot q_-)(\nu \cdot q_+))
        + 4 m_q^2 \sin^2{\theta_W} (3-2 \sin^2{\theta_W}) (N \cdot \nu)
    \right)
$$

##### Charged channel reaction

\begin{align}
    N + l^+ &\to u + \overline{d}
\end{align}

\begin{align}
    \require{cancel}
    \imath \mathcal{M} = \left(\overline{l} \gamma^\mu \frac{-\imath g}{2\sqrt{2}} P_L \frac{\imath M_D}{\cancel{N}} N \right)
        \frac{-\imath g_{\mu \nu}}{M_W^2} \frac{\imath g}{\sqrt{2}} \left(\overline{u} \gamma^\nu P_L  d \right)
\end{align}

where the momenta of particles denoted by the particle symbols: $l, N, q_\pm$

\begin{align}
    \require{cancel}
    |\overline{\mathcal{M}}(N l^+ \to u \overline{d})|^2 &= \left( \frac{g}{M_W} \right)^4 \frac{|\theta|^2}{32 m_N^2}
        Tr\left[ \cancel{N} (\cancel{N} \pm m_N) \cancel{N} P_R \gamma_\nu (\cancel{l}+m_l) \gamma^\mu P_L \right]
        Tr\left[ (\cancel{u} - m_u) \gamma_\mu P_L ( \cancel{d} + m_d) P_R \gamma^\nu \right]
    \\\\ &= \frac12 \left( \frac{g}{M_W} \right)^4 |\theta|^2 \left((d \cdot l) (N \cdot u) \right)
\end{align}

Together with the charge-conjugated channel this gives:

\begin{align}
    \require{cancel}
    |\overline{\mathcal{M}}|^2 &= \frac12 \left( \frac{g}{M_W} \right)^4 |\theta|^2 \left(
        (d \cdot l) (N \cdot u) + (u \cdot l) (N \cdot d)
    \right)
\end{align}