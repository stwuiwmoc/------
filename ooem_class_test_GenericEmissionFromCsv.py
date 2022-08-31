# %%
import importlib

import matplotlib.pyplot as plt

import object_oriented_estimation_myclass as ooem

if __name__ == "__main__":
    importlib.reload(ooem)

    # 各インスタンス生成
    light_all = ooem.LightGenenrator(
        rambda_division_width=0.1e-9,
        rambda_lower_limit=3.3e-6,
        rambda_upper_limit=3.5e-6)

    Jupiter_surface = ooem.GenericEmissionFromCsv(
        csv_fpath="mkfolder/convert_Norwood_etal_2015_fig2/rambda_vs_spectral_radiance.csv")

    # plot
    fig1 = plt.figure()
