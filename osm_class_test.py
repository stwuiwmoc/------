# %%
import observation_simulation_myclass as osm
import importlib

if __name__ == "__main__":
    importlib.reload(osm)

    N_H3p_ion = 5.0e+15
    T_hypothesis = 600

    instrument_params_for_Q_1_0 = osm.InstrumentParameters(
        is_ESPRIT=True,
        N_read=100,
        I_dark=20,
        G_Amp=9,
        has_fiber=True,
        l_fb=15,
        rambda_fi_center=3.953e-6,
        tau_fi_center=0.9,
        FWHM_fi=266e-9)

    instrument_params_for_R_3_0 = osm.InstrumentParameters(
        is_ESPRIT=True,
        N_read=100,
        I_dark=20,
        G_Amp=9,
        has_fiber=True,
        l_fb=15,
        rambda_fi_center=3.414e-6,
        tau_fi_center=0.9,
        FWHM_fi=17e-9)

    telescope_params = osm.TelescopeParameters(
        T_GBT=273,
        telescope_diameter=0.6)

    observation_params = osm.ObservationParameters(
        tau_alpha=0.812,
        t_obs=30 * 60,
        T_sky=273)

    Q_1_0_params = osm.EmissionLineParameters(
        rambda=3.9530e-6,
        N_H3p=N_H3p_ion,
        g_ns=4,
        J_prime=1,
        A_if=128.7,
        E_prime=2552.5691,
        T_hypo=T_hypothesis)

    R_3_0_params = osm.EmissionLineParameters(
        rambda=3.4128e-6,
        N_H3p=N_H3p_ion,
        g_ns=4,
        J_prime=4,
        A_if=177.6,
        E_prime=3382.9299,
        T_hypo=T_hypothesis)

    Q_1_0_obs = osm.EmissionLineDisperse(
        emission_line_params=Q_1_0_params,
        instrument_params=instrument_params_for_Q_1_0,
        telescope_params=telescope_params,
        observation_params=observation_params)

    R_3_0_obs = osm.EmissionLineDisperse(
        emission_line_params=R_3_0_params,
        instrument_params=instrument_params_for_R_3_0,
        telescope_params=telescope_params,
        observation_params=observation_params)

    temp_decision = osm.TemperatureFromSpectroscopy(
        emission_disperse_FD=Q_1_0_obs,
        emission_disperse_HB=R_3_0_obs)
