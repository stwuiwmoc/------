# %%
import importlib

import matplotlib.pyplot as plt

import object_oriented_estimation_myclass as ooem

if __name__ == "__main__":
    importlib.reload(ooem)

    # 各インスタンス生成
    light_all = ooem.LightGenenrator(
        rambda_division_width=0.1e-9,
        rambda_lower_limit=2.6e-6,
        rambda_upper_limit=5.0e-6)

    Jupiter_surface = ooem.GenericEmissionFromCsv(
        csv_fpath="mkfolder/convert_Norwood_etal_2015_fig2/rambda_vs_spectral_radiance.csv")

    # plot作成の準備
    fig1 = plt.figure(figsize=(5, 10))
    gs1 = fig1.add_gridspec(2, 1)

    # 木星表面の赤外発光を加える
    Jupiter_surface.add_spectral_radiance_to(light_instance=light_all)
    ax11 = light_all.show_rambda_vs_I_prime_plot(fig=fig1, position=gs1[0, 0])
    ax12 = light_all.show_rambda_vs_I_prime_plot(fig=fig1, position=gs1[1, 0])

    # plotの調整
    ax12.set_yscale("log")

    fig1.tight_layout()
