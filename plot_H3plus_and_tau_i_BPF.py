# %%
import importlib

import matplotlib.pyplot as plt

import object_oriented_estimation_myclass as ooem

if __name__ == "__main__":
    importlib.reload(ooem)

    # グローバル変数の定義
    column_density_H3plus = 1.5e+16  # [/m^2] H3+カラム密度
    T_thermospheric_H3plus = 1200  # [K] H3+熱圏温度

    # 各インスタンス生成
    light_H3plus = ooem.LightGenenrator(
        rambda_division_width=0.1e-9,
        rambda_lower_limit=3.3e-6,
        rambda_upper_limit=3.5e-6)

    light_BPF = ooem.LightGenenrator(
        rambda_division_width=0.1e-9,
        rambda_lower_limit=3.3e-6,
        rambda_upper_limit=3.5e-6)

    H3plus = ooem.H3ppyAuroralEmission(
        N_H3p=column_density_H3plus,
        T_hypo=T_thermospheric_H3plus,
        R_instrument=20000
    )

    Pirka = ooem.GroundBasedTelescope(
        D_GBT=1.6,
        FNO_GBT=12,
        T_GBT=268,
        tau_GBT=0.8**3)

    TOPICS = ooem.ImagingInstrument(
        is_TOPICS=True,
        rambda_BPF_center=3.413e-6,
        FWHM_BPF=17e-9,
        tau_BPF_center=0.88,
        tau_i_ND=1 / 0.9**6,  # tau_iをBPFの透過率だけにしたいので
        G_Amp=9,
        FW=168800,
        eta=0.889,
        I_dark=1.86,
        N_e_read=900)

    fits = ooem.VirtualOutputFileGenerator()

    # plot作成の準備
    fig1 = plt.figure(figsize=(8, 6))
    gs1 = fig1.add_gridspec(1, 1)

    # 輝線発光のプロット
    H3plus.add_auroral_emission_to(light_instance=light_H3plus)
    ax11 = light_H3plus.show_rambda_vs_I_prime_plot(fig=fig1, position=gs1[0, 0])

    # BPFのプロット
    # 透過率をプロットしたいので発光強度1の光とする
    light_BPF.add_I_prime_to(I_prime_xx=1)

    # 望遠鏡への設置
    TOPICS.set_ImagingInstrument_to(GBT_instance=Pirka)

    TOPICS.shoot_light_and_save_to_fits(
        light_instance=light_BPF,
        virtual_output_file_instance=fits,
        t_obs=1
    )

    # 透過率のプロット
    ax11_2 = ax11.twinx()
    ax11_2.plot(
        light_BPF.get_rambda(),
        light_BPF.get_I_prime(),
        color="red",
        linewidth=3
    )

    # グラフの調整
    ax11.set_xlabel(ax11.get_xlabel(), fontsize=20)
    ax11.set_ylabel(ax11.get_ylabel(), fontsize=20)
    ax11.set_xticklabels(ax11.get_xticklabels(), fontsize=15)
    ax11.set_yticklabels(ax11.get_yticklabels(), fontsize=15)

    ax11_2.set_ylabel("tau_BPF", fontsize=20)
    ax11_2.set_ylim(ax11_2.get_ylim()[0], 1)
