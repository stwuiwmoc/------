# %%
import importlib
import numpy as np
import matplotlib.pyplot as plt

import object_oriented_estimation_myclass as ooem


if __name__ == "__main__":
    importlib.reload(ooem)

    # グローバル変数の定義
    column_density_H3plus = 2.0e+16  # [/m^2] H3+カラム密度
    T_thermospheric_H3plus = 1200  # [K] H3+熱圏温度

    # 各インスタンス生成
    light = ooem.LightGenenrator(
        rambda_division_width=0.1e-9,
        rambda_lower_limit=3.3e-6,
        rambda_upper_limit=3.5e-6)

    R_3_0 = ooem.H3plusAuroralEmission(
        rambda_obj=3.4128e-6,
        N_H3p=column_density_H3plus,
        g_ns=4,
        J_prime=4,
        A_if=177.6,
        E_prime=3382.9299,
        T_hypo=T_thermospheric_H3plus)

    R_3_1 = ooem.H3plusAuroralEmission(
        rambda_obj=3.4149e-6,
        N_H3p=column_density_H3plus,
        g_ns=2,
        J_prime=4,
        A_if=110.4,
        E_prime=3359.002,
        T_hypo=T_thermospheric_H3plus)

    R_3_2 = ooem.H3plusAuroralEmission(
        rambda_obj=3.4207e-6,
        N_H3p=column_density_H3plus,
        g_ns=2,
        J_prime=4,
        A_if=86.91,
        E_prime=3287.2629,
        T_hypo=T_thermospheric_H3plus)

    R_3_3 = ooem.H3plusAuroralEmission(
        rambda_obj=3.4270e-6,
        N_H3p=column_density_H3plus,
        g_ns=4,
        J_prime=4,
        A_if=42.94,
        E_prime=3169.252,
        T_hypo=T_thermospheric_H3plus)

    R_4_3 = ooem.H3plusAuroralEmission(
        rambda_obj=3.4547e-6,
        N_H3p=column_density_H3plus,
        g_ns=4,
        J_prime=5,
        A_if=62.91,
        E_prime=3489.2151,
        T_hypo=T_thermospheric_H3plus)

    R_4_4 = ooem.H3plusAuroralEmission(
        rambda_obj=3.4548e-6,
        N_H3p=column_density_H3plus,
        g_ns=2,
        J_prime=5,
        A_if=122.9,
        E_prime=3332.4121,
        T_hypo=T_thermospheric_H3plus)

    T60 = ooem.GroundBasedTelescope(
        D_GBT=0.6,
        FNO_GBT=12,
        T_GBT=280,
        tau_GBT=0.66)

    TOPICS = ooem.ImagingInstrument(
        rambda_fl_center=3.414e-6,
        FWHM_fl=17e-9,
        tau_fl_center=0.88,
        G_Amp=9,
        I_dark=50,
        N_read=1200)

    fits = ooem.VirtualOutputFileGenerator()

    # 輝線発光を加える
    R_3_0.add_auroral_emission_to(light_instance=light)
    print("Only R_3_0 emission, I = ", light.get_I())

    R_3_1.add_auroral_emission_to(light_instance=light)

    R_3_2.add_auroral_emission_to(light_instance=light)

    R_3_3.add_auroral_emission_to(light_instance=light)

    R_4_3.add_auroral_emission_to(light_instance=light)

    R_4_4.add_auroral_emission_to(light_instance=light)
    light.show_rambda_vs_I_prime_plot()

    # 望遠鏡を通る
    T60.pass_through(light_instance=light)
    light.show_rambda_vs_I_prime_plot()

    # 望遠鏡への撮像装置の設置
    TOPICS.set_ImagingInstrument_to(GBT_instance=T60)

    # 撮像してfitsに保存
    TOPICS.shoot_light_and_save_to_fits(
        light_instance=light,
        virtual_output_file_instance=fits,
        t_obs=30)

    light.show_rambda_vs_I_prime_plot()

    print("S_all_pix =", fits.get_S_all_pix())
    print("S_FW_pix =", fits.get_S_FW_pix())

    # osmでの結果と一致するかのテスト
    test_obj_light = ooem.LightGenenrator(
        rambda_lower_limit=3.3e-6,
        rambda_upper_limit=3.5e-6,
        rambda_division_width=1e-9)

    test_sky_light = ooem.LightGenenrator(
        rambda_lower_limit=3.3e-6,
        rambda_upper_limit=3.5e-6,
        rambda_division_width=1e-9)

    # 大気透過率を計算するクラスはまだ実装してないのでGBTのクラスで代用
    # osmでは大気透過率、大気の熱輻射もガウシアンとして近似しているのでこれでokのはず
    # 観測見積もり.mdから、名寄11月、PWV=8000umでの透過率を使用
    Nayoro_good = ooem.GroundBasedTelescope(
        D_GBT=1,
        FNO_GBT=1,
        T_GBT=273,
        tau_GBT=0.286)

    test_obj_image = ooem.VirtualOutputFileGenerator()
    test_sky_image = ooem.VirtualOutputFileGenerator()

    test_t_obs = np.arange(0, 100)

    # --- 観測対象の撮像 ---
    # osmと同じくR(3, 0)輝線発光だけを考える
    R_3_0.add_auroral_emission_to(light_instance=test_obj_light)
    test_obj_light.show_rambda_vs_I_prime_plot()

    Nayoro_good.pass_through(light_instance=test_obj_light)
    test_obj_light.show_rambda_vs_I_prime_plot()

    T60.pass_through(light_instance=test_obj_light)
    test_obj_light.show_rambda_vs_I_prime_plot()

    TOPICS.shoot_light_and_save_to_fits(
        light_instance=test_obj_light,
        virtual_output_file_instance=test_obj_image,
        t_obs=test_t_obs)

    test_obj_light.show_rambda_vs_I_prime_plot()

    # --- skyイメージの撮像 ---
    Nayoro_good.pass_through(light_instance=test_sky_light)

    T60.pass_through(light_instance=test_sky_light)

    TOPICS.shoot_light_and_save_to_fits(
        light_instance=test_sky_light,
        virtual_output_file_instance=test_sky_image,
        t_obs=test_t_obs)
    test_sky_light.show_rambda_vs_I_prime_plot()

    # --- plot ---
    fig1 = plt.figure(figsize=(11, 8))
    gs1 = fig1.add_gridspec(2, 2)

    ax11 = fig1.add_subplot(gs1[0, 0])
    ax11.plot(test_t_obs, test_obj_image.get_S_all_pix(), label="obj")
    ax11.plot(test_t_obs, test_sky_image.get_S_all_pix(), label="sky")
    ax11.hlines(
        y=test_obj_image.get_S_FW_pix(),
        xmin=test_t_obs.min(),
        xmax=test_t_obs.max())
    ax11.grid()
    ax11.set_ylabel("Signal [DN]")
    ax11.legend()

    ax12 = fig1.add_subplot(gs1[0, 1])
    ax12.plot(
        test_t_obs,
        test_obj_image.get_S_all_pix() * TOPICS.get_G_sys(),
        label="obj")
    ax12.plot(
        test_t_obs,
        np.sqrt(test_sky_image.get_S_all_pix()) * TOPICS.get_G_sys(),
        label="sky")
    ax12.hlines(
        y=test_obj_image.get_S_FW_pix() * TOPICS.get_G_sys(),
        xmin=test_t_obs.min(),
        xmax=test_t_obs.max())
    ax12.set_xscale("log")
    ax12.set_yscale("log")
    ax12.grid(which="both")
    ax12.set_ylabel("Signal [e-]")
    ax12.legend()

    ax13 = fig1.add_subplot(gs1[1, 1])
    ax13.plot(
        test_t_obs,
        (test_obj_image.get_S_all_pix() - test_sky_image.get_S_all_pix()) / np.sqrt(test_sky_image.get_S_all_pix()))
    ax13.set_xscale("log")
    ax13.grid(which="both")
    ax13.set_ylabel("S/N")
