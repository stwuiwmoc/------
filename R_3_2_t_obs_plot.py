# %%
import numpy as np
import matplotlib.pyplot as plt
import importlib

import observation_simulation_myclass as osm
import R_3_0_t_obs_plot


def mkfolder(suffix=""):
    import os
    """
    Parameters
    ----------
    suffix : str, optional
        The default is "".

    Returns
    -------
    str ( script name + suffix )
    """
    filename = os.path.basename(__file__)
    filename = filename.replace(".py", "") + suffix
    folder = "mkfolder/" + filename + "/"
    os.makedirs(folder, exist_ok=True)
    return folder


if __name__ == "__main__":
    importlib.reload(osm)
    importlib.reload(R_3_0_t_obs_plot)
    t_obs_array = np.arange(1, 100)  # [s]

    R_3_2 = osm.EmissionLineParameters(
        rambda=3.4207e-6,
        N_H3p=2.0e+16,
        g_ns=2,
        J_prime=4,
        A_if=86.91,
        E_prime=3287.2629,
        T_hypo=1200)

    TOPICS = osm.InstrumentParameters(
        is_ESPRIT=False,
        N_read=1200,
        I_dark=20,
        G_Amp=9,
        has_fiber=False,
        l_fb=None,
        rambda_fl_center=3.414e-6,
        tau_fl_center=0.9,
        FWHM_fl=17e-9)

    T60_PWV2000 = osm.TelescopeParameters(
        T_GBT=273,
        telescope_diameter=0.6,
        tau_GBT=0.66,
        T_sky=273,
        tau_sky=0.967)

    T60_PWV5000 = osm.TelescopeParameters(
        T_GBT=273,
        telescope_diameter=0.6,
        tau_GBT=0.66,
        T_sky=273,
        tau_sky=0.965)

    Pirika_Oct = osm.TelescopeParameters(
        T_GBT=273,
        telescope_diameter=0.6,  # T60とF=12が同じなので、イメージスケールの詳細を実装し終わるまではこの値で代入
        tau_GBT=0.66,
        T_sky=273,
        tau_sky=0.900)

    Pirika_Nov = osm.TelescopeParameters(
        T_GBT=273,
        telescope_diameter=0.6,  # T60とF=12が同じなので、イメージスケールの詳細を実装し終わるまではこの値で代入
        tau_GBT=0.66,
        T_sky=273,
        tau_sky=0.911)

    obs_1bin = osm.ObservationParameters(
        n_bin_spatial=1,
        t_obs=t_obs_array)

    obs_4bin = osm.ObservationParameters(
        n_bin_spatial=2 * 2,
        t_obs=t_obs_array)

    obs_16bin = osm.ObservationParameters(
        n_bin_spatial=4 * 4,
        t_obs=t_obs_array)

    result_1bin_T60_PWV2000 = osm.EmissionLineDisperse(
        emission_line_params=R_3_2,
        instrument_params=TOPICS,
        telescope_params=T60_PWV2000,
        observation_params=obs_1bin)

    result_1bin_T60_PWV5000 = osm.EmissionLineDisperse(
        emission_line_params=R_3_2,
        instrument_params=TOPICS,
        telescope_params=T60_PWV5000,
        observation_params=obs_1bin)

    result_4bin_T60_PWV2000 = osm.EmissionLineDisperse(
        emission_line_params=R_3_2,
        instrument_params=TOPICS,
        telescope_params=T60_PWV2000,
        observation_params=obs_4bin)

    result_4bin_T60_PWV5000 = osm.EmissionLineDisperse(
        emission_line_params=R_3_2,
        instrument_params=TOPICS,
        telescope_params=T60_PWV5000,
        observation_params=obs_4bin)

    result_1bin_Pirika_Oct = osm.EmissionLineDisperse(
        emission_line_params=R_3_2,
        instrument_params=TOPICS,
        telescope_params=Pirika_Oct,
        observation_params=obs_1bin)

    result_1bin_Pirika_Nov = osm.EmissionLineDisperse(
        emission_line_params=R_3_2,
        instrument_params=TOPICS,
        telescope_params=Pirika_Nov,
        observation_params=obs_1bin)

    result_16bin_Pirika_Oct = osm.EmissionLineDisperse(
        emission_line_params=R_3_2,
        instrument_params=TOPICS,
        telescope_params=Pirika_Oct,
        observation_params=obs_16bin)

    result_16bin_Pirika_Nov = osm.EmissionLineDisperse(
        emission_line_params=R_3_2,
        instrument_params=TOPICS,
        telescope_params=Pirika_Nov,
        observation_params=obs_16bin)

    # Fullwell Limit detection
    FW_limit_index_T60_PWV2000 = R_3_0_t_obs_plot.find_index_of_full_well_limit(
        S_FW_pix=TOPICS.S_FW_pix, S_all_pix_array=result_1bin_T60_PWV2000.S_all_pix)

    FW_limit_index_T60_PWV5000 = R_3_0_t_obs_plot.find_index_of_full_well_limit(
        S_FW_pix=TOPICS.S_FW_pix, S_all_pix_array=result_1bin_T60_PWV5000.S_all_pix)

    FW_limit_index_Pirika_Oct = R_3_0_t_obs_plot.find_index_of_full_well_limit(
        S_FW_pix=TOPICS.S_FW_pix, S_all_pix_array=result_1bin_Pirika_Oct.S_all_pix)

    FW_limit_index_Pirika_Nov = R_3_0_t_obs_plot.find_index_of_full_well_limit(
        S_FW_pix=TOPICS.S_FW_pix, S_all_pix_array=result_1bin_Pirika_Nov.S_all_pix)

    # plot start
    # plot T60 PWV=2000
    fig1 = plt.figure(figsize=(12, 10))
    gs1 = fig1.add_gridspec(2, 2)
    fig1.suptitle("T60, PWV=2000")

    ax13 = R_3_0_t_obs_plot.plot_input_data_table_plot(
        fig=fig1,
        position=gs1[0:2, 1],
        emission_line_params=R_3_2,
        telescope_params=T60_PWV2000,
        instrument_params=TOPICS,
        result_1bin=result_1bin_T60_PWV2000)

    ax11 = R_3_0_t_obs_plot.plot_t_obs_vs_Signal_and_Noise_per_1_pixel(
        fig=fig1,
        position=gs1[0, 0],
        t_obs_array_=t_obs_array,
        instrument_params=TOPICS,
        result_1bin=result_1bin_T60_PWV2000,
        FW_limit_index=FW_limit_index_T60_PWV2000)

    ax12 = R_3_0_t_obs_plot.plot_t_obs_vs_SNR(
        fig=fig1,
        position=gs1[1, 0],
        t_obs_array_=t_obs_array,
        result_1bin=result_1bin_T60_PWV2000,
        result_nbin=result_4bin_T60_PWV2000,
        FW_limit_index=FW_limit_index_T60_PWV2000)

    fig1.tight_layout()

    fig1.savefig(mkfolder() + "fig1.png")

    # plot T60 PWV=5000
    fig2 = plt.figure(figsize=(12, 10))
    gs2 = fig2.add_gridspec(2, 2)
    fig2.suptitle("T60, PWV=5000")

    ax23 = R_3_0_t_obs_plot.plot_input_data_table_plot(
        fig=fig2,
        position=gs2[0:2, 1],
        emission_line_params=R_3_2,
        telescope_params=T60_PWV5000,
        instrument_params=TOPICS,
        result_1bin=result_1bin_T60_PWV5000)

    ax21 = R_3_0_t_obs_plot.plot_t_obs_vs_Signal_and_Noise_per_1_pixel(
        fig=fig2,
        position=gs2[0, 0],
        t_obs_array_=t_obs_array,
        instrument_params=TOPICS,
        result_1bin=result_1bin_T60_PWV5000,
        FW_limit_index=FW_limit_index_T60_PWV5000)

    ax22 = R_3_0_t_obs_plot.plot_t_obs_vs_SNR(
        fig=fig2,
        position=gs2[1, 0],
        t_obs_array_=t_obs_array,
        result_1bin=result_1bin_T60_PWV5000,
        result_nbin=result_4bin_T60_PWV5000,
        FW_limit_index=FW_limit_index_T60_PWV5000)

    fig2.tight_layout()

    fig2.savefig(mkfolder() + "fig2.png")

    # plot Pirika October
    fig3 = plt.figure(figsize=(12, 10))
    gs3 = fig3.add_gridspec(2, 2)
    fig3.suptitle("Pirika, October")

    ax33 = R_3_0_t_obs_plot.plot_input_data_table_plot(
        fig=fig3,
        position=gs3[0:2, 1],
        emission_line_params=R_3_2,
        telescope_params=Pirika_Oct,
        instrument_params=TOPICS,
        result_1bin=result_1bin_Pirika_Oct)

    ax31 = R_3_0_t_obs_plot.plot_t_obs_vs_Signal_and_Noise_per_1_pixel(
        fig=fig3,
        position=gs3[0, 0],
        t_obs_array_=t_obs_array,
        instrument_params=TOPICS,
        result_1bin=result_1bin_Pirika_Oct,
        FW_limit_index=FW_limit_index_Pirika_Oct)

    ax32 = R_3_0_t_obs_plot.plot_t_obs_vs_SNR(
        fig=fig3,
        position=gs3[1, 0],
        t_obs_array_=t_obs_array,
        result_1bin=result_1bin_Pirika_Oct,
        result_nbin=result_16bin_Pirika_Oct,
        FW_limit_index=FW_limit_index_Pirika_Oct)

    fig3.tight_layout()

    fig3.savefig(mkfolder() + "fig3.png")

    # plot Pirika November
    fig4 = plt.figure(figsize=(12, 10))
    gs4 = fig4.add_gridspec(2, 2)
    fig4.suptitle("Pirika, November")

    ax43 = R_3_0_t_obs_plot.plot_input_data_table_plot(
        fig=fig4,
        position=gs4[0:2, 1],
        emission_line_params=R_3_2,
        telescope_params=Pirika_Nov,
        instrument_params=TOPICS,
        result_1bin=result_1bin_Pirika_Nov)

    ax41 = R_3_0_t_obs_plot.plot_t_obs_vs_Signal_and_Noise_per_1_pixel(
        fig=fig4,
        position=gs4[0, 0],
        t_obs_array_=t_obs_array,
        instrument_params=TOPICS,
        result_1bin=result_1bin_Pirika_Nov,
        FW_limit_index=FW_limit_index_Pirika_Nov)

    ax42 = R_3_0_t_obs_plot.plot_t_obs_vs_SNR(
        fig=fig4,
        position=gs4[1, 0],
        t_obs_array_=t_obs_array,
        result_1bin=result_1bin_Pirika_Nov,
        result_nbin=result_16bin_Pirika_Nov,
        FW_limit_index=FW_limit_index_Pirika_Nov)

    fig4.tight_layout()

    fig4.savefig(mkfolder() + "fig4.png")
