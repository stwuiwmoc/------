# %%
import importlib

import matplotlib.pyplot as plt

import object_oriented_estimation_myclass as ooem


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
    importlib.reload(ooem)

    # グローバル変数の定義
    column_density_H3plus = 1.5e+16  # [/m^2] H3+カラム密度
    T_thermospheric_H3plus = 1200  # [K] H3+熱圏温度
    t_obs = 3  # [s] 積分時間
    n_bin_spatial_list = [2, 4, 8]

    # 各インスタンス生成
    light_all = ooem.LightGenenrator(
        rambda_division_width=0.1e-9,
        rambda_lower_limit=3.3e-6,
        rambda_upper_limit=3.5e-6)

    light_sky = ooem.LightGenenrator(
        rambda_division_width=0.1e-9,
        rambda_lower_limit=3.3e-6,
        rambda_upper_limit=3.5e-6)

    H3plus = ooem.H3ppyAuroralEmission(
        N_H3p=column_density_H3plus,
        T_hypo=T_thermospheric_H3plus,
        R_instrument=20000
    )

    Haleakala_Oct_good = ooem.EarthAtmosphere(
        T_ATM=273,
        ATRAN_result_filepath="raw_data/Ha_PWV2000_ZA22_Range3to4_R0.txt")

    PLANETS = ooem.GroundBasedTelescope(
        D_GBT=1.8,
        FNO_GBT=12,
        T_GBT=278,
        tau_GBT=0.8**7)

    ESPRIT = ooem.ImagingInstrument(
        is_TOPICS=False,
        rambda_BPF_center=3.413e-6,
        FWHM_BPF=17e-9,
        tau_BPF_center=0.88,
        tau_i_ND=1,
        G_Amp=9,
        I_dark=2.0,
        N_e_read=1500)

    fits_all = ooem.VirtualOutputFileGenerator()
    fits_sky = ooem.VirtualOutputFileGenerator()

    # plot作成の準備
    fig1 = plt.figure(figsize=(15, 15))
    gs1 = fig1.add_gridspec(4, 2)
    # M論用プロット
    fig2 = plt.figure(figsize=(12, 8))
    gs2 = fig2.add_gridspec(2, 2)

    # 輝線発光を加える
    H3plus.add_auroral_emission_to(light_instance=light_all)

    ax11 = light_all.show_rambda_vs_I_prime_plot(fig=fig1, position=gs1[0, 0])
    ax21 = light_all.show_rambda_vs_I_prime_plot(fig=fig2, position=gs2[0, 0])

    # 地球大気を通る
    Haleakala_Oct_good.pass_through(light_instance=light_all)
    Haleakala_Oct_good.pass_through(light_instance=light_sky)
    ax12 = light_all.show_rambda_vs_I_prime_plot(fig=fig1, position=gs1[1, 0])
    ax22 = light_all.show_rambda_vs_I_prime_plot(fig=fig2, position=gs2[0, 1])

    # 望遠鏡を通る
    PLANETS.pass_through(light_instance=light_all)
    PLANETS.pass_through(light_instance=light_sky)
    ax13 = light_all.show_rambda_vs_I_prime_plot(fig=fig1, position=gs1[2, 0])
    ax23 = light_all.show_rambda_vs_I_prime_plot(fig=fig2, position=gs2[1, 0])

    # 望遠鏡への撮像装置の設置
    ESPRIT.set_ImagingInstrument_to(GBT_instance=PLANETS)

    # 撮像してfitsに保存
    ESPRIT.shoot_light_and_save_to_fits(
        light_instance=light_all,
        virtual_output_file_instance=fits_all,
        t_obs=t_obs)
    ESPRIT.shoot_light_and_save_to_fits(
        light_instance=light_sky,
        virtual_output_file_instance=fits_sky,
        t_obs=t_obs)

    ax14 = light_all.show_rambda_vs_I_prime_plot(fig=fig1, position=gs1[3, 0])
    ax24 = light_all.show_rambda_vs_I_prime_plot(fig=fig2, position=gs2[1, 1])

    # SNRの計算
    SNRCalc = ooem.SNRCalculator(
        all_image_instance=fits_all,
        sky_image_instance=fits_sky,
        n_bin_spatial=n_bin_spatial_list[0])

    # plot
    ax11.set_title("H3+ emission lines")
    ax12.set_title("pass thruogh Earth Atmosphre")
    ax13.set_title("Pass through Ground-based-telescope")
    ax14.set_title("Pass through imaging instrument")

    parametar_table_list = [
        ooem.get_latest_commit_datetime(),
        ["Have some change", "from above commit", ooem.have_some_change_in_git_status()],
        ["", "", ""],
        ["H3plusEmission", "", ""],
        ["N_H3+", H3plus.get_N_H3p(), "/ m^2"],
        ["T_hypo", H3plus.get_T_hypo(), "K"],
        ["R_instrument", H3plus.get_R_instrument(), ""],
        ["", "", ""],
        ["EarthAtmosphre", "", ""],
        ["ATRAN result filename", Haleakala_Oct_good.get_ATRAN_result_filepath()[9:25], Haleakala_Oct_good.get_ATRAN_result_filepath()[25:]],
        ["T_ATM", Haleakala_Oct_good.get_T_ATM(), "K"],
        ["", "", ""],
        ["GroundBasedTelescope", "", ""],
        ["D_GBT", PLANETS.get_D_GBT(), "m"],
        ["FNO_GBT", PLANETS.get_FNO_GBT(), ""],
        ["tau_GBT", PLANETS.get_tau_GBT(), "K"],
        ["T_GBT", PLANETS.get_T_GBT(), "K"],
        ["", "", ""],
        ["ImagingInstrument", "", ""],
        ["rambda_BPF_center", ESPRIT.get_rambda_BPF_center(), "m"],
        ["FWHM_BPF", ESPRIT.get_FWHM_BPF(), "m"],
        ["tau_BPF_center", ESPRIT.get_tau_BPF_center(), ""],
        ["tau_i_ND", ESPRIT.get_tau_i_ND(), ""],
        ["I_dark_pix", ESPRIT.get_I_dark(), "e- / s"],
        ["N_e_read_pix", ESPRIT.get_N_e_read(), "e-rms"],
        ["", "", ""],
        ["Other parameters", "", ""],
        ["t_obs", t_obs, "s"],
        ["n_bin_spatial", n_bin_spatial_list[0], "pix"],
        ["theta_pix", ESPRIT.get_theta_pix(), "arcsec"],
        ["Field of View", ESPRIT.get_theta_pix() * 256, "arcsec"],
        ["spatial_resolution", SNRCalc.get_spatial_resolution(), "arcsec"],
        ["", "", ""],
        ["Results", "", ""],
        ["S_all_pix (obj image)", fits_all.get_S_all_pix(), "DN / pix"],
        ["S_all_pix (sky image)", fits_sky.get_S_all_pix(), "DN / pix"],
        ["S_dark_pix", fits_all.get_S_dark_pix(), "DN / pix"],
        ["S_read_pix", fits_all.get_S_read_pix(), "DN / pix"],
        ["R_electron/FW", fits_all.get_R_electron_FW(), ""],
        ["SNR", SNRCalc.get_SNR(), ""],
    ]

    ax15 = ooem.plot_parameter_table(
        fig=fig1, position=gs1[:, 1], parameter_table=parametar_table_list, fontsize=12)

    fig1.suptitle("H3+ 3.4um in Haleakala")
    fig1.tight_layout()
    fig1.savefig(mkfolder() + "fig1.png")

    # M論用プロット
    ax21.set_title("H3+ emission lines")
    ax22.set_title("pass thruogh Earth Atmosphre")
    ax23.set_title("Pass through Ground-based-telescope")
    ax24.set_title("Pass through imaging instrument")
    ax21.set_ylim(0, ax21.get_ylim()[1])
    ax22.set_ylim(0, ax21.get_ylim()[1])
    ax23.set_ylim(0, ax21.get_ylim()[1])
    ax24.set_ylim(0, ax21.get_ylim()[1])

    fig2.suptitle("H3+ 3.4um in Haleakala")
    fig2.tight_layout()

    # binning数を変えた時の空間分解能とSNRを表示
    for i in range(len(n_bin_spatial_list)):

        SNRCalc_for_nbin = ooem.SNRCalculator(
            all_image_instance=fits_all,
            sky_image_instance=fits_sky,
            n_bin_spatial=n_bin_spatial_list[i]
        )

        print(
            "n_bin_spatial =",
            n_bin_spatial_list[i],
            ", spatial_resolution =",
            SNRCalc_for_nbin.get_spatial_resolution(),
            ", SNR =",
            SNRCalc_for_nbin.get_SNR()
        )

        del SNRCalc_for_nbin
