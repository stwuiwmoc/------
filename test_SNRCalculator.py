# %%
import importlib

import object_oriented_estimation_myclass as ooem

if __name__ == "__main__":
    importlib.reload(ooem)

    # グローバル変数の定義
    t_obs = 1  # [s]
    n_bin_spatial = 2

    # インスタンス生成
    light_all = ooem.LightGenenrator(
        rambda_lower_limit=2e-6,
        rambda_upper_limit=4e-6,
        rambda_division_width=1e-9
    )

    light_sky = ooem.LightGenenrator(
        rambda_lower_limit=2e-6,
        rambda_upper_limit=4e-6,
        rambda_division_width=1e-9
    )

    Pirka = ooem.GroundBasedTelescope(
        D_GBT=1.6,
        FNO_GBT=12,
        T_GBT=273,
        tau_GBT=1
    )

    TOPICS = ooem.ImagingInstrument(
        is_TOPICS=True,
        rambda_BPF_center=3e-6,
        FWHM_BPF=1,
        tau_BPF_center=1,
        tau_i_ND=1,
        G_Amp=9,
        I_dark=50,
        N_e_read=1200
    )

    fits_all = ooem.VirtualOutputFileGenerator()
    fits_sky = ooem.VirtualOutputFileGenerator()

    TOPICS.set_ImagingInstrument_to(GBT_instance=Pirka)

    # 観測の模擬
    light_all.add_I_prime_to(I_prime_xx=1.0e3)
    TOPICS.shoot_light_and_save_to_fits(
        light_instance=light_all,
        virtual_output_file_instance=fits_all,
        t_obs=t_obs
    )

    light_sky.add_I_prime_to(I_prime_xx=0.5e3)
    TOPICS.shoot_light_and_save_to_fits(
        light_instance=light_sky,
        virtual_output_file_instance=fits_sky,
        t_obs=t_obs
    )

    SNRCalc = ooem.SNRCalculator(
        all_image_instance=fits_all,
        sky_image_instance=fits_sky,
        n_bin_spatial=n_bin_spatial
    )
