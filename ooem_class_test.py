# %%
import importlib
import matplotlib.pyplot as plt

import object_oriented_estimation_myclass as ooem


if __name__ == "__main__":
    importlib.reload(ooem)

    light = ooem.LightGenenrator(
        rambda_division_width=0.1e-9,
        rambda_lower_limit=3.4e-6,
        rambda_upper_limit=3.5e-6)

    # plot
    fig1 = plt.figure()
    gs1 = fig1.add_gridspec(1, 1)
    ax11 = fig1.add_subplot(gs1[0, 0])
