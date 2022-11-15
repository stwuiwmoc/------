# %%
import importlib

import matplotlib.pyplot as plt
import numpy as np

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

    # ===========================================================================
    # グローバル変数の定義
    serial_name_1 = "1000KSun"  # アウトバースト（1000K）発生、太陽反射光あり
    serial_name_2 = "500KSun"  # 火山平穏状態（500K）、太陽反射光あり

    Io_input_filepath_1 = "mkfolder/convert_de_Kleer_etal_2014_fig1/" + serial_name_1 + "_rambda_vs_spectral_radiance.csv"
    Io_input_filepath_2 = "mkfolder/convert_de_Kleer_etal_2014_fig1/" + serial_name_2 + "_rambda_vs_spectral_radiance.csv"

    t_obs = 15  # [s] 積分時間
    n_bin_spatial_list = [4, 8, 1]

    # ===========================================================================
    # 各インスタンス生成
    light_all = ooem.LightGenenrator(
        rambda_division_width=0.1e-9,
        rambda_lower_limit=2.05e-6,
        rambda_upper_limit=2.45e-6)

    light_sky = ooem.LightGenenrator(
        rambda_division_width=0.1e-9,
        rambda_lower_limit=2.05e-6,
        rambda_upper_limit=2.45e-6)

    Io_continuum_1 = ooem.GenericEmissionFromCsv(
        csv_fpath=Io_input_filepath_1)

    Io_continuum_2 = ooem.GenericEmissionFromCsv(
        csv_fpath=Io_input_filepath_2)

    Nayoro_Nov = ooem.EarthAtmosphere(
        T_ATM=273,
        ATRAN_result_filepath="raw_data/Na_PWV8000_ZA46_Range2to2.5_R0.txt")

    Pirka = ooem.GroundBasedTelescope(
        D_GBT=1.6,
        FNO_GBT=12,
        T_GBT=268,
        tau_GBT=0.8**3)

    TOPICS = ooem.ImagingInstrument(
        is_TOPICS=True,
        rambda_BPF_center=2.295e-6,
        FWHM_BPF=10e-9,
        tau_BPF_center=0.7,
        tau_i_ND=1,
        G_Amp=9,
        I_dark=50,
        N_e_read=1200)

    fits_all_1 = ooem.VirtualOutputFileGenerator()
    fits_sky_1 = ooem.VirtualOutputFileGenerator()

    fits_all_2 = ooem.VirtualOutputFileGenerator()
    fits_sky_2 = ooem.VirtualOutputFileGenerator()

    # 望遠鏡への撮像装置の設置
    TOPICS.set_ImagingInstrument_to(GBT_instance=Pirka)

    # ===========================================================================
    # アウトバースト（1000K）発生、太陽反射光あり の観測

    # plot作成の準備
    fig1 = plt.figure(figsize=(15, 15))
    gs1 = fig1.add_gridspec(4, 2)

    # イオの発光を加える
    Io_continuum_1.add_spectral_radiance_to(light_instance=light_all)
    ax11 = light_all.show_rambda_vs_I_prime_plot(fig=fig1, position=gs1[0, 0])

    # 地球大気を通る
    Nayoro_Nov.pass_through(light_instance=light_all)
    ax12 = light_all.show_rambda_vs_I_prime_plot(fig=fig1, position=gs1[1, 0])

    # 望遠鏡を通る
    Pirka.pass_through(light_instance=light_all)
    ax13 = light_all.show_rambda_vs_I_prime_plot(fig=fig1, position=gs1[2, 0])

    # 撮像してfitsに保存
    TOPICS.shoot_light_and_save_to_fits(
        light_instance=light_all,
        virtual_output_file_instance=fits_all_1,
        t_obs=t_obs)
    ax14 = light_all.show_rambda_vs_I_prime_plot(fig=fig1, position=gs1[3, 0])

    # sky画像の撮像（観測対象の撮像と同じ手順）
    Nayoro_Nov.pass_through(light_instance=light_sky)
    Pirka.pass_through(light_instance=light_sky)
    TOPICS.shoot_light_and_save_to_fits(
        light_instance=light_sky,
        virtual_output_file_instance=fits_sky_1,
        t_obs=t_obs)

    # SNRの計算
    SNRCalc_1 = ooem.SNRCalculator(
        all_image_instance=fits_all_1,
        sky_image_instance=fits_sky_1,
        n_bin_spatial=n_bin_spatial_list[0])

    # ===========================================================================
    # 火山平穏（500K）、太陽反射光あり の観測

    # インスタンスのリセット
    light_all.multiply_I_prime_to(magnification=0)
    light_sky.multiply_I_prime_to(magnification=0)

    # plotの準備
    fig2 = plt.figure(figsize=(15, 15))
    gs2 = fig2.add_gridspec(4, 2)

    # アウトバーストありと同じ手順
    # 観測対象の撮像
    Io_continuum_2.add_spectral_radiance_to(light_instance=light_all)
    ax21 = light_all.show_rambda_vs_I_prime_plot(fig=fig2, position=gs2[0, 0])

    Nayoro_Nov.pass_through(light_instance=light_all)
    ax22 = light_all.show_rambda_vs_I_prime_plot(fig=fig2, position=gs2[1, 0])

    Pirka.pass_through(light_instance=light_all)
    ax23 = light_all.show_rambda_vs_I_prime_plot(fig=fig2, position=gs2[2, 0])

    TOPICS.shoot_light_and_save_to_fits(
        light_instance=light_all,
        virtual_output_file_instance=fits_all_2,
        t_obs=t_obs
    )
    ax24 = light_all.show_rambda_vs_I_prime_plot(fig=fig2, position=gs2[3, 0])

    # sky画像の撮像
    Nayoro_Nov.pass_through(light_instance=light_sky)
    Pirka.pass_through(light_instance=light_sky)
    TOPICS.shoot_light_and_save_to_fits(
        light_instance=light_sky,
        virtual_output_file_instance=fits_sky_2,
        t_obs=t_obs
    )

    # SNRの計算
    SNRCalc_2 = ooem.SNRCalculator(
        all_image_instance=fits_all_2,
        sky_image_instance=fits_sky_2,
        n_bin_spatial=n_bin_spatial_list[0])

    # ===========================================================================
    # アウトバーストあり/なしでの差分のSN計算

    # Signal の差分
    S_obj_diff = SNRCalc_1.get_S_obj() - SNRCalc_2.get_S_obj()
    print("S_obj_diff =", S_obj_diff)

    # 誤差伝搬（参考 http://www.tagen.tohoku.ac.jp/labo/ishijima/gosa-03.html）
    N_diff = np.sqrt(SNRCalc_1.get_N_all()**2 + SNRCalc_2.get_N_all())
    print("N_diff =", N_diff)

    # SNR の計算
    SNR_diff = S_obj_diff / N_diff
    print("SNR_diff =", SNR_diff)

    # ===========================================================================
    # アウトバーストありのplotの調整
    ax11.set_title("Io thermal continuum")
    ax12.set_title("pass thruogh Earth Atmosphere")
    ax13.set_title("Pass through Ground-based-telescope")
    ax14.set_title("Pass through imaging instrument")

    ax11.set_ylim(0, ax11.get_ylim()[1])
    ax12.set_ylim(0, ax11.get_ylim()[1])
    ax13.set_ylim(0, ax13.get_ylim()[1])
    ax14.set_ylim(0, ax13.get_ylim()[1])

    parametar_table_list_1 = [
        ooem.get_latest_commit_datetime(),
        ["Have some change", "from above commit", ooem.have_some_change_in_git_status()],
        ["", "", ""],
        ["GenericEmissionFromCsv", "", ""],
        ["serial_name", serial_name_1, ""],
        ["", "", ""],
        ["EarthAtmosphre", "", ""],
        ["ATRAN result filename", Nayoro_Nov.get_ATRAN_result_filepath()[9:25], Nayoro_Nov.get_ATRAN_result_filepath()[25:]],
        ["T_ATM", Nayoro_Nov.get_T_ATM(), "K"],
        ["", "", ""],
        ["GroundBasedTelescope", "", ""],
        ["D_GBT", Pirka.get_D_GBT(), "m"],
        ["FNO_GBT", Pirka.get_FNO_GBT(), ""],
        ["tau_GBT", Pirka.get_tau_GBT(), "K"],
        ["T_GBT", Pirka.get_T_GBT(), "K"],
        ["", "", ""],
        ["ImagingInstrument", "", ""],
        ["rambda_BPF_center", TOPICS.get_rambda_BPF_center(), "m"],
        ["FWHM_BPF", TOPICS.get_FWHM_BPF(), "m"],
        ["tau_BPF_center", TOPICS.get_tau_BPF_center(), ""],
        ["tau_i_ND", TOPICS.get_tau_i_ND(), ""],
        ["I_dark_pix", TOPICS.get_I_dark(), "e- / s"],
        ["N_e_read_pix", TOPICS.get_N_e_read(), "e-rms"],
        ["", "", ""],
        ["Other parameters", "", ""],
        ["t_obs", t_obs, "s"],
        ["n_bin_spatial", n_bin_spatial_list[0], "pix"],
        ["theta_pix", TOPICS.get_theta_pix(), "arcsec"],
        ["Field of View", TOPICS.get_theta_pix() * 256, "arcsec"],
        ["spatial_resolution", SNRCalc_1.get_spatial_resolution(), "arcsec"],
        ["", "", ""],
        ["Results", "", ""],
        ["S_all_pix (obj image)", fits_all_1.get_S_all_pix(), "DN / pix"],
        ["S_all_pix (sky image)", fits_sky_1.get_S_all_pix(), "DN / pix"],
        ["S_dark_pix", fits_all_1.get_S_dark_pix(), "DN / pix"],
        ["S_read_pix", fits_all_1.get_S_read_pix(), "DN / pix"],
        ["R_electron/FW", fits_all_1.get_R_electron_FW(), ""],
        ["SNR", SNRCalc_1.get_SNR(), ""],
    ]

    ax15 = ooem.plot_parameter_table(
        fig=fig1, position=gs1[:, 1], parameter_table=parametar_table_list_1, fontsize=12)

    fig1.suptitle("Io continuum in" + serial_name_1)
    fig1.tight_layout()
    fig1.savefig(mkfolder() + "fig1.png")

    # 火山平穏時のplotの調整
    ax21.set_title(ax11.get_title())
    ax22.set_title(ax12.get_title())
    ax23.set_title(ax13.get_title())
    ax24.set_title(ax14.get_title())

    ax21.set_ylim(0, ax11.get_ylim()[1])
    ax22.set_ylim(0, ax11.get_ylim()[1])
    ax23.set_ylim(0, ax13.get_ylim()[1])
    ax24.set_ylim(0, ax13.get_ylim()[1])

    parametar_table_list_2 = [
        ooem.get_latest_commit_datetime(),
        ["Have some change", "from above commit", ooem.have_some_change_in_git_status()],
        ["", "", ""],
        ["GenericEmissionFromCsv", "", ""],
        ["serial_name", serial_name_2, ""],
        ["", "", ""],
        ["Results", "", ""],
        ["S_all_pix (obj image)", fits_all_2.get_S_all_pix(), "DN / pix"],
        ["S_all_pix (sky image)", fits_sky_2.get_S_all_pix(), "DN / pix"],
        ["S_dark_pix", fits_all_2.get_S_dark_pix(), "DN / pix"],
        ["S_read_pix", fits_all_2.get_S_read_pix(), "DN / pix"],
        ["R_electron/FW", fits_all_2.get_R_electron_FW(), ""],
        ["SNR", SNRCalc_2.get_SNR(), ""],
    ]

    ax25 = ooem.plot_parameter_table(
        fig=fig2, position=gs2[:, 1], parameter_table=parametar_table_list_2, fontsize=12)

    fig2.suptitle("Io continuum in " + serial_name_2)
    fig2.tight_layout()
    fig2.savefig(mkfolder() + "fig2.png")

    # ===========================================================================
    # binning数を変えた時の空間分解能とSNRを表示
    for i in range(len(n_bin_spatial_list)):

        SNRCalc_1_for_nbin = ooem.SNRCalculator(
            all_image_instance=fits_all_1,
            sky_image_instance=fits_sky_1,
            n_bin_spatial=n_bin_spatial_list[i]
        )

        print(
            "n_bin_spatial =",
            n_bin_spatial_list[i],
            ", spatial_resolution =",
            SNRCalc_1_for_nbin.get_spatial_resolution(),
            ", SNR =",
            SNRCalc_1_for_nbin.get_SNR()
        )

        del SNRCalc_1_for_nbin
