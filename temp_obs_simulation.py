# %%
import numpy as np
import matplotlib.pyplot as plt

import observation_simulation_myclass as osm


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


def calc_I_obj(physical_consts, line_params, N_H3p_: float, T_vib_: int):
    """calc_I_obj H3+の赤外輝線発光強度を計算

    Parameters
    ----------
    physical_consts : _type_
        PhysicalConstants
    line_params : _type_
        EmissionLineParameters
    N_H3p_ : float
        H3+カラム密度 [/m^2]
    T_vib_ : int
        H3+振動温度 [K]
        100 -> 1000まで100, 200, ... と100K刻みの温度か、1200K のみ入力可

    Returns
    -------
    float
        発光輝線強度 [W/m^2/str]
    """
    h = physical_consts.h
    c = physical_consts.c
    k_b = physical_consts.k_b
    Q_T = physical_consts.Q_T_dict[T_vib_]

    g_ns = line_params.g_ns
    J_prime = line_params.J_prime
    omega_if = line_params.omega_if
    A_if = line_params.A_if
    E_prime = line_params.E_prime

    I_obj_ = N_H3p_ * g_ns * (2 * J_prime + 1) * h * c * (omega_if * 1e2) * A_if\
        * np.exp(- h * c * (E_prime * 1e2) / (k_b * T_vib_)) \
        / (4 * np.pi * Q_T)
    return I_obj_


def calc_S_obj(physical_consts, line_params, I_obj_, A_t_, Omega_, tau_alpha_, tau_e_, eta_, G_sys_, t_obs_, n_pix_):
    h = physical_consts.h
    c = physical_consts.c
    lambda_ = line_params.lambda_um
    S_obj_ = ((I_obj_ * A_t_ * Omega_ * tau_alpha_ * tau_e_) / (h * c / (lambda_ * 1e-6))) * eta_ * (1 / G_sys_) * t_obs_ * n_pix_
    return S_obj_


def calc_del_T_del_R(physical_consts, line_params_fd, line_params_hb, beta_, R_s_):
    h = physical_consts.h
    c = physical_consts.c
    k_b = physical_consts.k_b
    E_prime_fd = line_params_fd.E_prime
    E_prime_hb = line_params_hb.E_prime
    del_T_del_R_ = - h * c / k_b * ((E_prime_hb - E_prime_fd) * 1e2) * (np.log(beta_) - np.log(R_s_)) ** -2 / R_s_
    return del_T_del_R_


def calc_SNR(S_obj_, S_GBT_sky_, S_dark_, N_read_, G_sys_, n_pix_):
    SNR_ = S_obj_ / np.sqrt(S_obj_ + S_GBT_sky_ + S_dark_ + (N_read_ / G_sys_)**2 * n_pix_)
    return SNR_


if __name__ == "__main__":

    # variables
    diamiter_list = [0.6, 1.8]
    delta_T_list = []
    SNR_fd_list = []
    SNR_hb_list = []

    N_read = np.arange(0, 1200)  # e-rms/pix
    I_dark = 50  # e-/s/pix
    G_sys = 10.9  # [e-/DN]
    t_obs = 60 * 15  # [s]
    T_vib = 700  # [K]
    N_H3p = 5.0e15  # [/m^2]

    # physical constants
    PHYSICAL_CONSTS = osm.PhysicalConstants()

    # instrument constants
    tau_alpha = 0.9  # 大気透過率
    tau_e = 0.1063  # 光学系全体の透過率
    eta = 0.889  # 検出器量子効率
    plate_scale = 0.3  # [arcsec/pixel]
    slit_width_arcsec = 0.7  # arcsec
    n_pix = slit_width_arcsec / plate_scale
    Omega = 2.35e-11 * plate_scale  # 画素が見込む立体角

    # input parameter

    beta = 4.795144  # Tvib導出の式のβ

    # Q(1, 0-)輝線
    line_fd = osm.EmissionLineParameters(
        lambda_um=3.9530,
        g_ns=4,
        J_prime=1,
        A_if=128.7,
        E_prime=2552.5691)

    # R(3, 0-)輝線
    line_hb = osm.EmissionLineParameters(
        lambda_um=3.4128,
        g_ns=4,
        J_prime=4,
        A_if=177.6,
        E_prime=3382.9299)

    I_obj_fd = calc_I_obj(PHYSICAL_CONSTS, line_fd, N_H3p, T_vib)
    I_obj_hb = calc_I_obj(PHYSICAL_CONSTS, line_hb, N_H3p, T_vib)
    I_GBT_sky = 2.99e-6  # GBT発光とSky発光（3.4um, ESPRITスリット幅）

    for diamiter in diamiter_list:

        # instrument constants
        A_t = (diamiter / 2) ** 2 * np.pi  # 望遠鏡開口面積

        S_obj_fd = calc_S_obj(PHYSICAL_CONSTS, line_fd, I_obj_fd, A_t, Omega, tau_alpha, tau_e, eta, G_sys, t_obs, n_pix)
        S_obj_hb = calc_S_obj(PHYSICAL_CONSTS, line_hb, I_obj_hb, A_t, Omega, tau_alpha, tau_e, eta, G_sys, t_obs, n_pix)

        S_GBT_sky_fd = calc_S_obj(PHYSICAL_CONSTS, line_fd, I_GBT_sky, A_t, Omega, 1, tau_e, eta, G_sys, t_obs, n_pix)  # I_GBT + I_skyは大気透過率を考慮しないため1としている
        S_GBT_sky_hb = calc_S_obj(PHYSICAL_CONSTS, line_hb, I_GBT_sky, A_t, Omega, 1, tau_e, eta, G_sys, t_obs, n_pix)  # I_GBT + I_skyは大気透過率を考慮しないため1としている
        S_dark = I_dark / G_sys * t_obs * n_pix

        SNR_fd = calc_SNR(S_obj_fd, S_GBT_sky_fd, S_dark, N_read, G_sys, n_pix)
        SNR_hb = calc_SNR(S_obj_hb, S_GBT_sky_hb, S_dark, N_read, G_sys, n_pix)

        R_s = S_obj_hb / S_obj_fd

        delta_S_obj_fd = S_obj_fd / SNR_fd
        delta_S_obj_hb = S_obj_hb / SNR_hb

        del_T_del_R = calc_del_T_del_R(PHYSICAL_CONSTS, line_fd, line_hb, beta, R_s)
        delta_R_s = np.sqrt((1 / S_obj_fd * delta_S_obj_hb) ** 2 + (S_obj_hb / (S_obj_fd ** 2) * delta_S_obj_fd) ** 2)
        delta_T = abs(del_T_del_R) * delta_R_s

        SNR_fd_list.append(SNR_fd)
        SNR_hb_list.append(SNR_hb)
        delta_T_list.append(delta_T)

    fig1 = plt.figure()
    gs1 = fig1.add_gridspec(1, 1)

    ax11 = fig1.add_subplot(gs1[0, 0])
    ax11.plot(N_read, delta_T_list[0], label="T60")
    ax11.plot(N_read, delta_T_list[1], label="PLANETS")
    ax11.grid()

    ax11.set_xlabel("N_read [e-rms]")
    ax11.set_ylabel("abs( T_error ) [K]")
    ax11.set_title(
        "Temperature dicision error\n" +
        "T_vib ~ " + str(T_vib) + " [K] and I_dark = " + str(I_dark) + " [e-/s]")
    ax11.legend()

    fig1.savefig(mkfolder() + "fig1.png")

    fig2 = plt.figure(figsize=(7, 7))
    gs2 = fig2.add_gridspec(2, 1)

    ax21 = fig2.add_subplot(gs2[0, 0])
    ax21.plot(N_read, SNR_fd_list[0], label="T60")
    ax21.plot(N_read, SNR_fd_list[1], label="PLANETS")
    ax21.grid()

    ax21.set_xlabel("N_read [e-rms]")
    ax21.set_ylabel("Line Emission S/N")
    ax21.set_title(
        "Fundamental Band Line Emission\n" +
        "T_vib ~ " + str(T_vib) + " [K] and I_dark = " + str(I_dark) + " [e-/s]")
    ax21.legend()

    ax22 = fig2.add_subplot(gs2[1, 0])
    ax22.plot(N_read, SNR_hb_list[0], label="T60")
    ax22.plot(N_read, SNR_hb_list[1], label="PLANETS")
    ax22.grid()

    ax22.set_xlabel("N_read [e-rms]")
    ax22.set_ylabel("Line Emission S/N")
    ax22.set_title(
        "Hot Band Line Emission\n" +
        "T_vib ~ " + str(T_vib) + " [K] and I_dark = " + str(I_dark) + " [e-/s]")
    ax22.legend()

    fig2.tight_layout()
