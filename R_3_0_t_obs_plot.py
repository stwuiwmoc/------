# %%
import importlib

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

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
        fig: matplotlib.figure.Figure,
        position: matplotlib.gridspec.GridSpec,
        t_obs_array_: np.ndarray,
        instrument_params: osm.InstrumentParameters,
        result_1bin: osm.EmissionLineDisperse,
        FW_limit_index: int) -> matplotlib.axes._subplots.Axes:
    """1piexlでのSignalとNoiseをe-単位でプロット

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figureオブジェクト
    position : matplotlib.gridspec.GridSpec
        Gridspecオブジェクト
    t_obs_array_ : np.ndarray
        [s] 積分時間
    instrument_params : osm.InstrumentParameters
        自作インスタンス
    result_1bin : osm.EmissionLineDisperse
        自作インスタンス
        binning数が1のもの
    FW_limit_index : int
        FWに達する時のindex

    Returns
    -------
    matplotlib.axes._subplots.Axes
        Axesオブジェクト
    """

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
        fig: matplotlib.figure.Figure,
        position: matplotlib.gridspec.GridSpec,
        t_obs_array_: np.ndarray,
        result_1bin: osm.EmissionLineDisperse,
        result_nbin_small: osm.EmissionLineDisperse,
        result_nbin_large: osm.EmissionLineDisperse,
        FW_limit_index: int) -> matplotlib.axes._subplots.Axes:

    """binningなし、bininng数少なめ、binning数多めの各SNRをプロット

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figureオブジェクト
    position : matplotlib.gridspec.GridSpec
        Gridspecオブジェクト
    t_obs_array_ : np.ndarray
        [s] 積分時間
    result_1bin : osm.EmissionLineDisperse
        自作インスタンス
        binning数が1のもの
    result_nbin_small : osm.EmissionLineDisperse
        自作インスタンス
        binning数が少なめのもの
    result_nbin_large : osm.EmissionLineDisperse
        自作インスタンス
        binning数が多めのもの
    FW_limit_index : int
        FWに達する時のindex

    Returns
    -------
    matplotlib.axes._subplots.Axes
        Axesオブジェクト
    """

    ax = fig.add_subplot(position)

    # SNR plot
    ax.plot(
        t_obs_array_,
        result_1bin.SNR,
        label=str(result_1bin.n_bin) + " pix binning")
    ax.plot(
        t_obs_array_,
        result_nbin_small.SNR,
        label=str(result_nbin_small.n_bin) + " pix binning")
    ax.plot(
        t_obs_array_,
        result_nbin_large.SNR,
        label=str(result_nbin_large.n_bin) + " pix binning")

    # FullWell limit plot
    ax.vlines(
        x=t_obs_array_[FW_limit_index],
        ymin=0,
        ymax=np.max(result_nbin_large.SNR),
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
        fig: matplotlib.figure.Figure,
        position: matplotlib.gridspec.GridSpec,
        emission_line_params: osm.EmissionLineParameters,
        telescope_params: osm.TelescopeParameters,
        instrument_params: osm.InstrumentParameters,
        result_1bin: osm.EmissionLineDisperse) -> matplotlib.axes._subplots.Axes:
    """_summary_

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figureオブジェクト
    position : matplotlib.gridspec.GridSpec
        Gridspecオブジェクト
    emission_line_params : osm.EmissionLineParameters
        自作インスタンス
    telescope_params : osm.TelescopeParameters
        自作インスタンス
    instrument_params : osm.InstrumentParameters
        自作インスタンス
    result_1bin : osm.EmissionLineDisperse
        自作インスタンス
        binning数が1のもの

    Returns
    -------
    matplotlib.axes._subplots.Axes
        Axesオブジェクト
    """

    # data table plot
    table_text_list = [
        osm.get_latest_commit_datetime(),
        ["Have some change", "from above commit", osm.have_some_change_in_git_status()],
        ["", "", ""],
        ["H3+ parameters", "", ""],
        ["rambda", emission_line_params.rambda, "m"],
        ["N_H3p", emission_line_params.N_H3p, "m^-2"],
        ["T_hypo", emission_line_params.T_hypo, "K"],
        ["", "", ""],
        ["Telescope parameters", "", ""],
        ["D_t", telescope_params.D_t, "m"],
        ["T_GBT", telescope_params.T_GBT, "K"],
        ["tau_GBT", telescope_params.tau_GBT, ""],
        ["T_sky", telescope_params.T_sky, "K"],
        ["tau_sky", telescope_params.tau_sky, ""],
        ["", "", ""],
        ["Observation_parameters", "", ""],
        ["t_obs", "", "s"],
        ["", "", ""],
        ["Intensity", "", ""],
        ["I_obj", result_1bin.I_obj, "W / m^2 / sr"],
        ["I_sky", result_1bin.I_sky, "W / m^2 / sr"],
        ["I_GBT", result_1bin.I_GBT, "W / m^2 / sr"],
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

    T60_PWV2000 = osm.TelescopeParameters(
        T_GBT=273,
        D_t=0.6,
        FNO_t=12,
        tau_GBT=0.66,
        T_sky=273,
        tau_sky=0.744)

    T60_PWV5000 = osm.TelescopeParameters(
        T_GBT=273,
        D_t=0.6,
        FNO_t=12,
        tau_GBT=0.66,
        T_sky=273,
        tau_sky=0.564)

    Pirika_Oct = osm.TelescopeParameters(
        T_GBT=273,
        D_t=1.6,
        FNO_t=12,
        tau_GBT=0.66,
        T_sky=273,
        tau_sky=0.131)

    Pirika_Nov = osm.TelescopeParameters(
        T_GBT=273,
        D_t=1.6,
        FNO_t=12,
        tau_GBT=0.66,
        T_sky=273,
        tau_sky=0.286)

    TOPICS_in_T60 = osm.InstrumentParameters(
        telescope_params=T60_PWV2000,
        is_ESPRIT=False,
        N_read=600,
        I_dark=20,
        G_Amp=9,
        has_fiber=False,
        l_fb=None,
        rambda_fl_center=3.414e-6,
        tau_fl_center=0.9,
        FWHM_fl=17e-9)

    TOPICS_in_Pirika = osm.InstrumentParameters(
        telescope_params=Pirika_Oct,
        is_ESPRIT=False,
        N_read=600,
        I_dark=20,
        G_Amp=9,
        has_fiber=False,
        l_fb=None,
        rambda_fl_center=3.414e-6,
        tau_fl_center=0.9,
        FWHM_fl=17e-9)

    obs_1bin = osm.ObservationParameters(
        n_bin_spatial=1,
        t_obs=t_obs_array)

    obs_9bin = osm.ObservationParameters(
        n_bin_spatial=3 * 3,
        t_obs=t_obs_array)

    obs_25bin = osm.ObservationParameters(
        n_bin_spatial=5 * 5,
        t_obs=t_obs_array)

    obs_64bin = osm.ObservationParameters(
        n_bin_spatial=8 * 8,
        t_obs=t_obs_array)

    obs_100bin = osm.ObservationParameters(
        n_bin_spatial=10 * 10,
        t_obs=t_obs_array)

    result_1bin_T60_PWV2000 = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS_in_T60,
        telescope_params=T60_PWV2000,
        observation_params=obs_1bin)

    result_1bin_T60_PWV5000 = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS_in_T60,
        telescope_params=T60_PWV5000,
        observation_params=obs_1bin)

    result_9bin_T60_PWV2000 = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS_in_T60,
        telescope_params=T60_PWV2000,
        observation_params=obs_9bin)

    result_9bin_T60_PWV5000 = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS_in_T60,
        telescope_params=T60_PWV5000,
        observation_params=obs_9bin)

    result_25bin_T60_PWV2000 = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS_in_T60,
        telescope_params=T60_PWV2000,
        observation_params=obs_25bin)

    result_25bin_T60_PWV5000 = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS_in_T60,
        telescope_params=T60_PWV5000,
        observation_params=obs_25bin)

    result_1bin_Pirika_Oct = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS_in_Pirika,
        telescope_params=Pirika_Oct,
        observation_params=obs_1bin)

    result_1bin_Pirika_Nov = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS_in_Pirika,
        telescope_params=Pirika_Nov,
        observation_params=obs_1bin)

    result_64bin_Pirika_Oct = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS_in_Pirika,
        telescope_params=Pirika_Oct,
        observation_params=obs_64bin)

    result_64bin_Pirika_Nov = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS_in_Pirika,
        telescope_params=Pirika_Nov,
        observation_params=obs_64bin)

    result_100bin_Pirika_Oct = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS_in_Pirika,
        telescope_params=Pirika_Oct,
        observation_params=obs_100bin)

    result_100bin_Pirika_Nov = osm.EmissionLineDisperse(
        emission_line_params=R_3_0,
        instrument_params=TOPICS_in_Pirika,
        telescope_params=Pirika_Nov,
        observation_params=obs_100bin)

    # Fullwell Limit detection
    FW_limit_index_T60_PWV2000 = find_index_of_full_well_limit(
        S_FW_pix=TOPICS_in_T60.S_FW_pix, S_all_pix_array=result_1bin_T60_PWV2000.S_all_pix)

    FW_limit_index_T60_PWV5000 = find_index_of_full_well_limit(
        S_FW_pix=TOPICS_in_T60.S_FW_pix, S_all_pix_array=result_1bin_T60_PWV5000.S_all_pix)

    FW_limit_index_Pirika_Oct = find_index_of_full_well_limit(
        S_FW_pix=TOPICS_in_Pirika.S_FW_pix, S_all_pix_array=result_1bin_Pirika_Oct.S_all_pix)

    FW_limit_index_Pirika_Nov = find_index_of_full_well_limit(
        S_FW_pix=TOPICS_in_Pirika.S_FW_pix, S_all_pix_array=result_1bin_Pirika_Nov.S_all_pix)

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
        instrument_params=TOPICS_in_T60,
        result_1bin=result_1bin_T60_PWV2000)

    ax11 = plot_t_obs_vs_Signal_and_Noise_per_1_pixel(
        fig=fig1,
        position=gs1[0, 0],
        t_obs_array_=t_obs_array,
        instrument_params=TOPICS_in_T60,
        result_1bin=result_1bin_T60_PWV2000,
        FW_limit_index=FW_limit_index_T60_PWV2000)

    ax12 = plot_t_obs_vs_SNR(
        fig=fig1,
        position=gs1[1, 0],
        t_obs_array_=t_obs_array,
        result_1bin=result_1bin_T60_PWV2000,
        result_nbin_small=result_9bin_T60_PWV2000,
        result_nbin_large=result_25bin_T60_PWV2000,
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
        instrument_params=TOPICS_in_T60,
        result_1bin=result_1bin_T60_PWV5000)

    ax21 = plot_t_obs_vs_Signal_and_Noise_per_1_pixel(
        fig=fig2,
        position=gs2[0, 0],
        t_obs_array_=t_obs_array,
        instrument_params=TOPICS_in_T60,
        result_1bin=result_1bin_T60_PWV5000,
        FW_limit_index=FW_limit_index_T60_PWV5000)

    ax22 = plot_t_obs_vs_SNR(
        fig=fig2,
        position=gs2[1, 0],
        t_obs_array_=t_obs_array,
        result_1bin=result_1bin_T60_PWV5000,
        result_nbin_small=result_9bin_T60_PWV5000,
        result_nbin_large=result_25bin_T60_PWV5000,
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
        instrument_params=TOPICS_in_Pirika,
        result_1bin=result_1bin_Pirika_Oct)

    ax31 = plot_t_obs_vs_Signal_and_Noise_per_1_pixel(
        fig=fig3,
        position=gs3[0, 0],
        t_obs_array_=t_obs_array,
        instrument_params=TOPICS_in_Pirika,
        result_1bin=result_1bin_Pirika_Oct,
        FW_limit_index=FW_limit_index_Pirika_Oct)

    ax32 = plot_t_obs_vs_SNR(
        fig=fig3,
        position=gs3[1, 0],
        t_obs_array_=t_obs_array,
        result_1bin=result_1bin_Pirika_Oct,
        result_nbin_small=result_64bin_Pirika_Oct,
        result_nbin_large=result_100bin_Pirika_Oct,
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
        instrument_params=TOPICS_in_Pirika,
        result_1bin=result_1bin_Pirika_Nov)

    ax41 = plot_t_obs_vs_Signal_and_Noise_per_1_pixel(
        fig=fig4,
        position=gs4[0, 0],
        t_obs_array_=t_obs_array,
        instrument_params=TOPICS_in_Pirika,
        result_1bin=result_1bin_Pirika_Nov,
        FW_limit_index=FW_limit_index_Pirika_Nov)

    ax42 = plot_t_obs_vs_SNR(
        fig=fig4,
        position=gs4[1, 0],
        t_obs_array_=t_obs_array,
        result_1bin=result_1bin_Pirika_Nov,
        result_nbin_small=result_64bin_Pirika_Nov,
        result_nbin_large=result_100bin_Pirika_Nov,
        FW_limit_index=FW_limit_index_Pirika_Nov)

    fig4.tight_layout()

    fig4.savefig(mkfolder() + "fig4.png")
