## Neutrinos, neutrons decoupling and BBN setup

Decoupling of the neutrinos and neutrons follows the similar pattern as the decoupling of the \
sterile neutrinos. Estimated decoupling temperature for both species is around $1.5 MeV$. However,\
it is important to discriminate them by more careful analysis that includes precise accounting for\
scattering cross-sections and regimes of the particles involved in reactions.

The most important reactions for neutrinos are

\begin{align}
    \nu + \nu &\to \nu + \nu \\\\
    \nu + e^- &\to \nu + e^-
\end{align}

(and crossing processes). Neutrons are governed by only 2 reactions:

\begin{align}
    p + e &\leftrightarrow n + \nu_e \\\\
    n &\to p + e + \overline{\nu}_e
\end{align}


The estimated temperature is relatively close to the electron mass, so the assumption of the \
all-relativistic particles in the plasma becomes approximate and the more precise treatment is \
required.

TODO: Derive neutron decoupling temperature using not only the $\Gamma \sim H$ trick, but also the\
    typical energy transfer timescale of the reactions.

After neutrinos decoupled, additionally to the contributions from the non-instantaneous decoupling,\
they gain corrections from the sterile neutrinos decay. One of the most important parts of the \
project is to calculate this correction spectrum. It depends both on the sterile neutrinos spectrum\
and decay channels.

Accordingly, corrections to the active neutrino spectrum influence the equilibrium of the neutrons\
reactions. For the scattering channel, additional active neutrinos lead to the shift in the balance\
to the production of $p$ and $e$, while in the decay channel they can increase the lifetime of the\
neutron:

\begin{align}
   \mathcal{F}\left[p+e\to n + \nu\right] &= f_p f_e (1-f_n)(1-f_\nu) - f_n f_\nu (1- f_p)(1-f_e) =
   \\\\ &= f_p f_e(1-f_n) - f_\nu \left[f_n(1-f_p-f_e) + f_p f_e \right]
\end{align}

\begin{align}
    \mathcal{F}\left[n \to p+e+\overline{\nu}\right] &= f_n (1- f_p)(1-f_e)(1 - f_\nu) - f_p f_e f_\nu (1-f_n) =
    \\\\ &= f_n (1-f_p)(1-f_e) - f_\nu \left[ f_n(1 - f_p - f_e) +f_p f_e \right]
\end{align}

$$ f_\nu \to f_\nu + \delta f_\nu $$

Other distributions are taken to be equilibrium, so the functional vanishes for them:


\begin{align}
    \mathcal{F}\left[p+e\to n + \nu\right] = \mathcal{F}\left[n \to p+e+\overline{\nu}\right] &= - \delta f_\nu \left[f_n(1-f_p-f_e) + f_p f_e \right]
\end{align}

The expression in the brackets is extremely close to $1$ when we consider non-relativistic protons\
and neutrons. So the Boltzmann integral is roughly proportional to the number of non-equiliubrium\
active neutrinos.

Physically, this effect contains 2 parts: firstly, more active neutrinos mean more interactions;\
secondly, active neutrinos phase space density repulses inverse reactions as they are fermions.\
One can imagine a situation when due to the densely packed neutrino phase space, neutrons are not \
allowed to decay.
