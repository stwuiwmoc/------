# %%
import numpy as np
import matplotlib.pyplot as plt

import observation_simulation_myclass as osm

if __name__ == "__main__":
    t_obs_array = np.arange(1, 100)  # [s]

    R_3_0_params = osm.EmissionLineParameters(
        rambda=3.4128e-6,
        N_H3p=2.0e+16,
        g_ns=4,
        J_prime=4,
        A_if=177.6,
        E_prime=3382.9299,
        T_hypo=1200)

    print("I_obj = ", R_3_0_params.I_obj)

    TOPICS_params = osm.InstrumentParameters(
        is_ESPRIT=False,
        N_read=100,
        I_dark=20,
        G_Amp=9,
        has_fiber=False,
        l_fb=None,
        rambda_fi_center=3.414e-6,
        tau_fi_center=0.9,
        FWHM_fi=17e-9)

    T60_params = osm.TelescopeParameters(
        T_GBT=273,
        telescope_diameter=0.6)

    print(
        "I_GBT = ",
        T60_params.calc_I_GBT(
            rambda_=R_3_0_params.rambda,
            FWHM_fi=TOPICS_params.FWHM_fi) * 1e4,
        "e-04")

    obs_1_bin_params = osm.ObservationParameters(
        tau_alpha=0.812,
        T_sky=273,
        n_bin_spatial=1,
        t_obs=t_obs_array)

    obs_4_bin_params = osm.ObservationParameters(
        tau_alpha=0.812,
        T_sky=273,
        n_bin_spatial=4,
        t_obs=t_obs_array)

    print(
        "I_sky = ",
        obs_1_bin_params.calc_I_sky(
            rambda_=R_3_0_params.rambda,
            FWHM_fi=TOPICS_params.FWHM_fi) * 1e4,
        "e-04")

    R_3_0_obs_1_bin = osm.EmissionLineDisperse(
        emission_line_params=R_3_0_params,
        instrument_params=TOPICS_params,
        telescope_params=T60_params,
        observation_params=obs_1_bin_params)

    R_3_0_obs_4_bin = osm.EmissionLineDisperse(
        emission_line_params=R_3_0_params,
        instrument_params=TOPICS_params,
        telescope_params=T60_params,
        observation_params=obs_4_bin_params)

    # plot start
    fig1 = plt.figure(figsize=(7, 10))
    gs1 = fig1.add_gridspec(2, 1)

    # t_obs vs e- plot
    ax11 = fig1.add_subplot(gs1[0, 0])

    # Signal plot
    linestyle_signal = "-"
    linewidth_signal = 3

    ax11.plot(
        t_obs_array,
        R_3_0_obs_1_bin.S_all * TOPICS_params.G_sys,
        label="S_all",
        linestyle=linestyle_signal,
        linewidth=linewidth_signal)
    ax11.plot(
        t_obs_array,
        R_3_0_obs_1_bin.S_obj * TOPICS_params.G_sys,
        label="S_obj",
        linestyle=linestyle_signal,
        linewidth=linewidth_signal)

    # Noise plot
    linestyle_noise = "-."
    ax11.plot(
        t_obs_array,
        np.sqrt(R_3_0_obs_1_bin.S_all * TOPICS_params.G_sys),
        label="N_all",
        linestyle=linestyle_noise)
    ax11.plot(
        t_obs_array,
        np.sqrt(R_3_0_obs_1_bin.S_sky * TOPICS_params.G_sys),
        label="N_sky",
        linestyle=linestyle_noise)
    ax11.plot(
        t_obs_array,
        np.sqrt(R_3_0_obs_1_bin.S_GBT * TOPICS_params.G_sys),
        label="N_GBT",
        linestyle=linestyle_noise)
    ax11.plot(
        t_obs_array,
        np.sqrt(R_3_0_obs_1_bin.S_dark * TOPICS_params.G_sys),
        label="N_dark",
        linestyle=linestyle_noise)
    ax11.plot(
        t_obs_array,
        TOPICS_params.N_read * np.ones(len(t_obs_array)),
        label="N_read",
        linestyle=linestyle_noise)

    # FullWell plot
    ax11.plot(
        t_obs_array,
        TOPICS_params.S_FW_pix * TOPICS_params.G_sys * np.ones(len(t_obs_array)),
        label="FW limit",
        linestyle=":",
        linewidth=3)

    # axes decoration
    ax11.grid(axis="x", which="major")
    ax11.grid(axis="y", which="major")

    ax11.set_xlabel("Integration Time t_obs [s]")
    ax11.set_ylabel("Signal and Noise per 1 pixel [e-]")

    ax11.set_xscale("log")
    ax11.set_yscale("log")

    ax11.legend()

    # t_obs vs SNR plot
    ax12 = fig1.add_subplot(gs1[1, 0])

    # SNR plot
    ax12.plot(
        t_obs_array,
        R_3_0_obs_1_bin.SNR,
        label=str(R_3_0_obs_1_bin.n_bin) + " pix binning")
    ax12.plot(
        t_obs_array,
        R_3_0_obs_4_bin.SNR,
        label=str(R_3_0_obs_4_bin.n_bin) + " pix binning")

    ax12.grid()

    ax12.set_xlabel("Integration Time t_obs [s]")
    ax12.set_ylabel("SNR of each binning number")

    ax12.set_xscale("log")

    ax12.legend()
