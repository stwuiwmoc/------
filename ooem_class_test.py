# %%
import importlib
import numpy as np

import object_oriented_estimation_myclass as ooem


if __name__ == "__main__":
    importlib.reload(ooem)

    light = ooem.LightGenenrator(
        rambda_division_width=0.1e-9,
        rambda_lower_limit=3.3e-6,
        rambda_upper_limit=3.5e-6)

    light.show_rambda_vs_I_prime_plot()

    light.add_I_prime_to(
        I_prime_xx=np.arange(light.get_len())
    )

    light.show_rambda_vs_I_prime_plot()

    light.multiply_I_prime_to(
        magnification=np.arange(light.get_len())
    )

    light.show_rambda_vs_I_prime_plot()
