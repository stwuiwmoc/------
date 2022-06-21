# %%
import numpy as np
import matplotlib.pyplot as plt
import importlib

import observation_simulation_myclass as osm


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


def find_index_of_full_well_limit(
        S_FW_pix: float,
        S_all_pix_array: np.ndarray) -> int:
    """full well limitになる時のindexを見つける

    検出器のFW以上には電荷を蓄積できないため、S_all_pix =< S_FW_pix となる必要がある。
    Full well limit は、S_all_pix = S_FW_pixとなった時なので、
    S_all_pix と S_FW_pix の差分が最小になる時のindexを取得すればよい

    Parameters
    ----------
    S_FW_pix_array : float
        [DN / pix] FWに達した時のカウント値
    S_all_pix_array : np.ndarray
        [DN / pix] カウント値の合計

    Returns
    -------
    int
        full well limitになる時のindex
    """

    # FWとallの差分をとる
    diff_all_FW_array = np.abs(S_all_pix_array - S_FW_pix)

    # 差分が最小となる時のindexを取得する
    index = np.argmin(diff_all_FW_array)

    return index


def plot_t_obs_vs_Signal_and_Noise_per_1_pixel(
        fig, position, t_obs_array_, instrument_params, result_1bin, FW_limit_index):

    ax = fig.add_subplot(position)

    # Signal plot
    linestyle_signal = "-"
    linewidth_signal = 3

    ax.plot(
        t_obs_array_,
        result_1bin.S_all * instrument_params.G_sys,
        label="S_all",
        linestyle=linestyle_signal,
        linewidth=linewidth_signal)
    ax.plot(
        t_obs_array_,
        result_1bin.S_obj * instrument_params.G_sys,
        label="S_obj",
        linestyle=linestyle_signal,
        linewidth=linewidth_signal)

    # Noise plot
    linestyle_noise = "-."
    ax.plot(
        t_obs_array_,
        np.sqrt(result_1bin.S_all * instrument_params.G_sys),
        label="N_all",
        linestyle=linestyle_noise)
    ax.plot(
        t_obs_array_,
        np.sqrt(result_1bin.S_sky * instrument_params.G_sys),
        label="N_sky",
        linestyle=linestyle_noise)
    ax.plot(
        t_obs_array_,
        np.sqrt(result_1bin.S_GBT * instrument_params.G_sys),
        label="N_GBT",
        linestyle=linestyle_noise)
    ax.plot(
        t_obs_array_,
        np.sqrt(result_1bin.S_dark * instrument_params.G_sys),
        label="N_dark",
        linestyle=linestyle_noise)
    ax.plot(
        t_obs_array_,
        instrument_params.N_read * np.ones(len(t_obs_array_)),
        label="N_read",
        linestyle=linestyle_noise)

    # FullWell plot
    ax.plot(
        t_obs_array_,
        instrument_params.S_FW_pix * instrument_params.G_sys * np.ones(len(t_obs_array_)),
        label="FW limit",
        linestyle=":",
        linewidth=3,
        color="red")

    # FullWell limit plot
    ax.vlines(
        x=t_obs_array_[FW_limit_index],
        ymin=1,
        ymax=np.max(result_1bin.S_all * instrument_params.G_sys),
        linestyle=":",
        linewidth=3,
        color="red")

    # axes decoration
    ax.grid(axis="x", which="major")
    ax.grid(axis="y", which="major")

    ax.set_xlabel("Integration Time t_obs [s]")
    ax.set_ylabel("Signal and Noise per 1 pixel [e-]")

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.legend()

    return ax


def plot_t_obs_vs_SNR(
        fig, position, t_obs_array_, result_1bin, result_nbin, FW_limit_index):

    ax = fig.add_subplot(position)

    # SNR plot
    ax.plot(
        t_obs_array_,
        result_1bin.SNR,
        label=str(result_1bin.n_bin) + " pix binning")
    ax.plot(
        t_obs_array_,
        result_nbin.SNR,
        label=str(result_nbin.n_bin) + " pix binning")

    # FullWell limit plot
    ax.vlines(
        x=t_obs_array_[FW_limit_index],
        ymin=0,
        ymax=np.max(result_nbin.SNR),
        linestyle=":",
        linewidth=3,
        color="red")

    ax.grid()

    ax.set_xlabel("Integration Time t_obs [s]")
    ax.set_ylabel("SNR of each binning number")

    ax.set_xscale("log")

    ax.legend()

    return ax


def plot_input_data_table_plot(
        fig,
        position,
        emission_line_params,
        telescope_params,
        instrument_params,
        result_1bin):

    # data table plot
    table_text_list = [
        ["H3+ parameters", "", ""],
        ["rambda", emission_line_params.rambda, "m"],
        ["N_H3p", emission_line_params.N_H3p, "m^-2"],
        ["T_hypo", emission_line_params.T_hypo, "K"],
        ["", "", ""],
        ["Telescope parameters", "", ""],
        ["telescope_diameter", telescope_params.telescope_diameter, "m"],
        ["T_GBT", telescope_params.T_GBT, "K"],
        ["tau_GBT", telescope_params.tau_GBT, ""],
        ["T_sky", telescope_params.T_sky, "K"],
        ["tau_sky", telescope_params.tau_sky, ""],
        ["", "", ""],
        ["Observation_parameters", "", ""],
        ["t_obs", "", "s"],
        ["", "", ""],
        ["Intensity", "", ""],
        ["I_obj", result_1bin.I_obj, "W / m^2 / str"],
        ["I_sky", result_1bin.I_sky, "W / m^2 / str"],
        ["I_GBT", result_1bin.I_GBT, "W / m^2 / str"],
        ["", "", ""],
        ["Instrument parameters", "", ""],
        ["N_read", instrument_params.N_read, "e-rms"],
        ["I_dark", instrument_params.I_dark, "e-/s"],
        ["rambda_fl_center", instrument_params.rambda_fl_center, "m"],
        ["tau_fl_center", instrument_params.tau_fl_center, ""],
        ["FWHM_fl", instrument_params.FWHM_fl, "m"]]

    ax = osm.plot_parameter_table(
        fig=fig,
        position=position,
        parameter_table=table_text_list,
        fontsize=8)

    return ax


if __name__ == "__main__":
    importlib.reload(osm)
    t_obs_array = np.arange(1, 100)  # [s]

    R_3_0 = osm.EmissionLineParameters(
        rambda=3.4128e-6,
        N_H3p=2.0e+16,
        g_ns=4,
        J_prime=4,
        A_if=177.6,
        E_prime=3382.9299,
        T_hypo=1200)

    TOPICS = osm.InstrumentParameters(
        is_ESPRIT=False,
        N_read=600,
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
        tau_sky=0.744)

    T60_PWV5000 = osm.TelescopeParameters(
        T_GBT=273,
        telescope_diameter=0.6,
        tau_GBT=0.66,
        T_sky=273,
        tau_sky=0.564)

    Pirika_Oct = osm.TelescopeParameters(
        T_GBT=273,
        telescope_diameter=0.6,  # T60とF=12が同じなので、イメージスケールの詳細を実装し終わるまではこの値で代入
        tau_GBT=0.66,
        T_sky=273,
        tau_sky=0.131)

    Pirika_Nov = osm.TelescopeParameters(
        T_GBT=273,
        telescope_diameter=0.6,  # T60とF=12が同じなので、イメージスケールの詳細を実装し終わるまではこの値で代入
        tau_GBT=0.66,
        T_sky=273,
        tau_sky=0.286)

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
        emission_line_params=R_3_0,
        instrument_params=TOPICS,
        telescope_params=T60_PWV2000,
        observation_params=obs_1bin)

    result_1bin_T60_PWV5000 = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS,
        telescope_params=T60_PWV5000,
        observation_params=obs_1bin)

    result_4bin_T60_PWV2000 = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS,
        telescope_params=T60_PWV2000,
        observation_params=obs_4bin)

    result_4bin_T60_PWV5000 = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS,
        telescope_params=T60_PWV5000,
        observation_params=obs_4bin)

    result_1bin_Pirika_Oct = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS,
        telescope_params=Pirika_Oct,
        observation_params=obs_1bin)

    result_1bin_Pirika_Nov = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS,
        telescope_params=Pirika_Nov,
        observation_params=obs_1bin)

    result_16bin_Pirika_Oct = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS,
        telescope_params=Pirika_Oct,
        observation_params=obs_16bin)

    result_16bin_Pirika_Nov = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS,
        telescope_params=Pirika_Nov,
        observation_params=obs_16bin)

    # Fullwell Limit detection
    FW_limit_index_T60_PWV2000 = find_index_of_full_well_limit(
        S_FW_pix=TOPICS.S_FW_pix, S_all_pix_array=result_1bin_T60_PWV2000.S_all_pix)

    FW_limit_index_T60_PWV5000 = find_index_of_full_well_limit(
        S_FW_pix=TOPICS.S_FW_pix, S_all_pix_array=result_1bin_T60_PWV5000.S_all_pix)

    FW_limit_index_Pirika_Oct = find_index_of_full_well_limit(
        S_FW_pix=TOPICS.S_FW_pix, S_all_pix_array=result_1bin_Pirika_Oct.S_all_pix)

    FW_limit_index_Pirika_Nov = find_index_of_full_well_limit(
        S_FW_pix=TOPICS.S_FW_pix, S_all_pix_array=result_1bin_Pirika_Nov.S_all_pix)

    # plot start
    # plot T60 PWV=2000
    fig1 = plt.figure(figsize=(12, 10))
    gs1 = fig1.add_gridspec(2, 2)
    fig1.suptitle("T60, PWV=2000")

    ax13 = plot_input_data_table_plot(
        fig=fig1,
        position=gs1[0:2, 1],
        emission_line_params=R_3_0,
        telescope_params=T60_PWV2000,
        instrument_params=TOPICS,
        result_1bin=result_1bin_T60_PWV2000)

    ax11 = plot_t_obs_vs_Signal_and_Noise_per_1_pixel(
        fig=fig1,
        position=gs1[0, 0],
        t_obs_array_=t_obs_array,
        instrument_params=TOPICS,
        result_1bin=result_1bin_T60_PWV2000,
        FW_limit_index=FW_limit_index_T60_PWV2000)

    ax12 = plot_t_obs_vs_SNR(
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

    ax23 = plot_input_data_table_plot(
        fig=fig2,
        position=gs2[0:2, 1],
        emission_line_params=R_3_0,
        telescope_params=T60_PWV5000,
        instrument_params=TOPICS,
        result_1bin=result_1bin_T60_PWV5000)

    ax21 = plot_t_obs_vs_Signal_and_Noise_per_1_pixel(
        fig=fig2,
        position=gs2[0, 0],
        t_obs_array_=t_obs_array,
        instrument_params=TOPICS,
        result_1bin=result_1bin_T60_PWV5000,
        FW_limit_index=FW_limit_index_T60_PWV5000)

    ax22 = plot_t_obs_vs_SNR(
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

    ax33 = plot_input_data_table_plot(
        fig=fig3,
        position=gs3[0:2, 1],
        emission_line_params=R_3_0,
        telescope_params=Pirika_Oct,
        instrument_params=TOPICS,
        result_1bin=result_1bin_Pirika_Oct)

    ax31 = plot_t_obs_vs_Signal_and_Noise_per_1_pixel(
        fig=fig3,
        position=gs3[0, 0],
        t_obs_array_=t_obs_array,
        instrument_params=TOPICS,
        result_1bin=result_1bin_Pirika_Oct,
        FW_limit_index=FW_limit_index_Pirika_Oct)

    ax32 = plot_t_obs_vs_SNR(
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

    ax43 = plot_input_data_table_plot(
        fig=fig4,
        position=gs4[0:2, 1],
        emission_line_params=R_3_0,
        telescope_params=Pirika_Nov,
        instrument_params=TOPICS,
        result_1bin=result_1bin_Pirika_Nov)

    ax41 = plot_t_obs_vs_Signal_and_Noise_per_1_pixel(
        fig=fig4,
        position=gs4[0, 0],
        t_obs_array_=t_obs_array,
        instrument_params=TOPICS,
        result_1bin=result_1bin_Pirika_Nov,
        FW_limit_index=FW_limit_index_Pirika_Nov)

    ax42 = plot_t_obs_vs_SNR(
        fig=fig4,
        position=gs4[1, 0],
        t_obs_array_=t_obs_array,
        result_1bin=result_1bin_Pirika_Nov,
        result_nbin=result_16bin_Pirika_Nov,
        FW_limit_index=FW_limit_index_Pirika_Nov)

    fig4.tight_layout()

    fig4.savefig(mkfolder() + "fig4.png")
