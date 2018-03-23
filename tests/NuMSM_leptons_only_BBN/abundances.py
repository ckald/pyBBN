import numpy as np
from matplotlib import pyplot as plt


########### 4He ABUNDANCE CHECKING ###########

Params = [
            (2e1, 1e-2, 0.24720),
            (2e1, 1e-4, 0.26618),
            (5e1, 1e-4, 0.24592),
            (5e1, 1e-6, 0.28641),
            (1e2, 1e-6, 0.24654),
            (1e2, 1e-8, 0.31388)
]

helium_abundance_lower = 0.2452
helium_abundance_upper = 0.2696


########### Plots of active neutrino non-equilibrium distribution functions ###########

pyBBN_mom_SBBN = np.arange(0, 51, 1)

HNL_active_e_SBBN = np.loadtxt("data_HNL_active/Electron_neutrino.distribution_SBBN.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Electron_neutrino.distribution_SBBN.txt", skiprows=2)[0][2:]

HNL_active_mu_SBBN = np.loadtxt("data_HNL_active/Muon_neutrino.distribution_SBBN.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Muon_neutrino.distribution_SBBN.txt", skiprows=2)[0][2:]


pyBBN_mom = np.arange(0, 101, 1)

HNL_active_e_20_2 = np.loadtxt("data_HNL_active/Electron_neutrino.distribution_20_1e-2.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Electron_neutrino.distribution_20_1e-2.txt", skiprows=2)[0][2:]

HNL_active_e_20_4 = np.loadtxt("data_HNL_active/Electron_neutrino.distribution_20_1e-4.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Electron_neutrino.distribution_20_1e-4.txt", skiprows=2)[0][2:]

HNL_active_e_50_4 = np.loadtxt("data_HNL_active/Electron_neutrino.distribution_50_1e-4.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Electron_neutrino.distribution_50_1e-4.txt", skiprows=2)[0][2:]

HNL_active_e_50_6 = np.loadtxt("data_HNL_active/Electron_neutrino.distribution_50_1e-6.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Electron_neutrino.distribution_50_1e-6.txt", skiprows=2)[0][2:]

HNL_active_e_100_6 = np.loadtxt("data_HNL_active/Electron_neutrino.distribution_100_1e-6.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Electron_neutrino.distribution_100_1e-6.txt", skiprows=2)[0][2:]

HNL_active_e_100_8 = np.loadtxt("data_HNL_active/Electron_neutrino.distribution_100_1e-8.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Electron_neutrino.distribution_100_1e-8.txt", skiprows=2)[0][2:]


HNL_active_mu_20_2 = np.loadtxt("data_HNL_active/Muon_neutrino.distribution_20_1e-2.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Muon_neutrino.distribution_20_1e-2.txt", skiprows=2)[0][2:]

HNL_active_mu_20_4 = np.loadtxt("data_HNL_active/Muon_neutrino.distribution_20_1e-4.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Muon_neutrino.distribution_20_1e-4.txt", skiprows=2)[0][2:]

HNL_active_mu_50_4 = np.loadtxt("data_HNL_active/Muon_neutrino.distribution_50_1e-4.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Muon_neutrino.distribution_50_1e-4.txt", skiprows=2)[0][2:]

HNL_active_mu_50_6 = np.loadtxt("data_HNL_active/Muon_neutrino.distribution_50_1e-6.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Muon_neutrino.distribution_50_1e-6.txt", skiprows=2)[0][2:]

HNL_active_mu_100_6 = np.loadtxt("data_HNL_active/Muon_neutrino.distribution_100_1e-6.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Muon_neutrino.distribution_100_1e-6.txt", skiprows=2)[0][2:]

HNL_active_mu_100_8 = np.loadtxt("data_HNL_active/Muon_neutrino.distribution_100_1e-8.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Muon_neutrino.distribution_100_1e-8.txt", skiprows=2)[0][2:]


HNL_active_tau_20_2 = np.loadtxt("data_HNL_active/Tau_neutrino.distribution_20_1e-2.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Tau_neutrino.distribution_20_1e-2.txt", skiprows=2)[0][2:]

HNL_active_tau_20_4 = np.loadtxt("data_HNL_active/Tau_neutrino.distribution_20_1e-4.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Tau_neutrino.distribution_20_1e-4.txt", skiprows=2)[0][2:]

HNL_active_tau_50_4 = np.loadtxt("data_HNL_active/Tau_neutrino.distribution_50_1e-4.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Tau_neutrino.distribution_50_1e-4.txt", skiprows=2)[0][2:]

HNL_active_tau_50_6 = np.loadtxt("data_HNL_active/Tau_neutrino.distribution_50_1e-6.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Tau_neutrino.distribution_50_1e-6.txt", skiprows=2)[0][2:]

HNL_active_tau_100_6 = np.loadtxt("data_HNL_active/Tau_neutrino.distribution_100_1e-6.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Tau_neutrino.distribution_100_1e-6.txt", skiprows=2)[0][2:]

HNL_active_tau_100_8 = np.loadtxt("data_HNL_active/Tau_neutrino.distribution_100_1e-8.txt", skiprows=2)[-1][2:] \
                / np.loadtxt("data_HNL_active/Tau_neutrino.distribution_100_1e-8.txt", skiprows=2)[0][2:]


plt.figure(figsize=(8,9))
plt.subplot(321)
plt.title(r"$M=20\ \mathrm{MeV},\ \theta^2=10^{-2}$", fontsize=8)
plt.plot(pyBBN_mom[:12], HNL_active_e_20_2[:12], color="orange", label=r"NuMSM$-e$")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_e_SBBN[:12], linestyle="--", color="orange", label=r"SBBN$-e$")

plt.plot(pyBBN_mom[:12], HNL_active_mu_20_2[:12], color="green", label=r"NuMSM$-\mu$")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_mu_SBBN[:12], linestyle="--", color="green", label=r"SBBN$-\mu$")

plt.plot(pyBBN_mom[:12], HNL_active_tau_20_2[:12], color="purple", label=r"NuMSM$-\tau$")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_mu_SBBN[:12], linestyle="--", color="purple", label=r"SBBN$-\tau$")

plt.subplots_adjust(hspace=0.4, wspace=0.3)
plt.legend(loc="upper left", frameon=False, prop={'size': 8})
plt.ylabel(r"$f_{\nu_\alpha} / f_\mathrm{eq}$", labelpad=8, fontsize=12)

plt.subplot(322)
plt.title(r"$M=20\ \mathrm{MeV},\ \theta^2=10^{-4}$", fontsize=8)
plt.plot(pyBBN_mom[:12], HNL_active_e_20_4[:12], color="orange")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_e_SBBN[:12], linestyle="--", color="orange")

plt.plot(pyBBN_mom[:12], HNL_active_mu_20_4[:12], color="green")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_mu_SBBN[:12], linestyle="--", color="green")

plt.plot(pyBBN_mom[:12], HNL_active_tau_20_4[:12], color="purple")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_mu_SBBN[:12], linestyle="--", color="purple")

plt.subplot(323)
plt.title(r"$M=50\ \mathrm{MeV},\ \theta^2=10^{-4}$", fontsize=8)
plt.plot(pyBBN_mom[:12], HNL_active_e_50_4[:12], color="orange")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_e_SBBN[:12], linestyle="--", color="orange")

plt.plot(pyBBN_mom[:12], HNL_active_mu_50_4[:12], color="green")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_mu_SBBN[:12], linestyle="--", color="green")

plt.plot(pyBBN_mom[:12], HNL_active_tau_50_4[:12], color="purple")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_mu_SBBN[:12], linestyle="--", color="purple")

plt.ylabel(r"$f_{\nu_\alpha} / f_\mathrm{eq}$", labelpad=8, fontsize=12)

plt.subplot(324)
plt.title(r"$M=50\ \mathrm{MeV},\ \theta^2=10^{-6}$", fontsize=8)
plt.plot(pyBBN_mom[:12], HNL_active_e_50_6[:12], color="orange")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_e_SBBN[:12], linestyle="--", color="orange")

plt.plot(pyBBN_mom[:12], HNL_active_mu_50_6[:12], color="green")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_mu_SBBN[:12], linestyle="--", color="green")

plt.plot(pyBBN_mom[:12], HNL_active_tau_50_6[:12], color="purple")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_mu_SBBN[:12], linestyle="--", color="purple")

plt.subplot(325)
plt.title(r"$M=100\ \mathrm{MeV},\ \theta^2=10^{-6}$", fontsize=8)
plt.plot(pyBBN_mom[:12], HNL_active_e_100_6[:12], color="orange")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_e_SBBN[:12], linestyle="--", color="orange")

plt.plot(pyBBN_mom[:12], HNL_active_mu_100_6[:12], color="green")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_mu_SBBN[:12], linestyle="--", color="green")

plt.plot(pyBBN_mom[:12], HNL_active_tau_100_6[:12], color="purple")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_mu_SBBN[:12], linestyle="--", color="purple")

plt.ylabel(r"$f_{\nu_\alpha} / f_\mathrm{eq}$", labelpad=8, fontsize=12)
plt.xlabel("$y=pa$ [MeV]", labelpad=10, fontsize=12)

plt.subplot(326)
plt.title(r"$M=100\ \mathrm{MeV},\ \theta^2=10^{-8}$", fontsize=8)
plt.plot(pyBBN_mom[:12], HNL_active_e_100_8[:12], color="orange")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_e_SBBN[:12], linestyle="--", color="orange")

plt.plot(pyBBN_mom[:12], HNL_active_mu_100_8[:12], color="green")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_mu_SBBN[:12], linestyle="--", color="green")

plt.plot(pyBBN_mom[:12], HNL_active_tau_100_8[:12], color="purple")
plt.plot(pyBBN_mom_SBBN[:12], HNL_active_mu_SBBN[:12], linestyle="--", color="purple")

plt.xlabel("$y=pa$ [MeV]", labelpad=10, fontsize=12)

plt.show(block=True)