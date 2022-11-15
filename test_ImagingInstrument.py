# %%
import importlib

import matplotlib.pyplot as plt

import object_oriented_estimation_myclass as ooem

if __name__ == "__main__":
    importlib.reload(ooem)

    # インスタンス生成

    light = ooem.LightGenenrator(
        rambda_lower_limit=1e-6,
        rambda_upper_limit=8e-6,
        rambda_division_width=1e-9
    )

    Pirka = ooem.GroundBasedTelescope(
        D_GBT=1.6,
        FNO_GBT=12,
        T_GBT=273,
        tau_GBT=0.8**3
    )

    ESPRIT = ooem.ImagingInstrument(
        is_TOPICS=False,
        rambda_BPF_center=2e-6,
        FWHM_BPF=1e-1,
        tau_BPF_center=1,
        tau_i_ND=1,
        G_Amp=9,
        I_dark=50,
        N_e_read=1200
    )

    fits = ooem.VirtualOutputFileGenerator()

    # 望遠鏡に観測装置を設置
    ESPRIT.set_ImagingInstrument_to(GBT_instance=Pirka)

    # 図の準備
    fig1 = plt.figure()
    gs1 = fig1.add_gridspec(2, 1)

    # ESPRIT透過率のテスト
    # 分かりやすいように元の光の発光強度は全波長で１とする
    light.add_I_prime_to(I_prime_xx=1)
    ax11 = light.show_rambda_vs_I_prime_plot(fig=fig1, position=gs1[0, 0])

    ESPRIT.shoot_light_and_save_to_fits(
        light_instance=light,
        virtual_output_file_instance=fits,
        t_obs=1
    )

    ax12 = light.show_rambda_vs_I_prime_plot(fig=fig1, position=gs1[1, 0])
