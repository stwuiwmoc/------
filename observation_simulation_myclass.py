import numpy as np
from scipy import constants as phys_consts


def mkhelp(instance):
    import inspect
    attr_list = list(instance.__dict__.keys())
    for attr in attr_list:
        if attr.startswith("_"):
            continue
        print(attr)
    for method in inspect.getmembers(instance, inspect.ismethod):
        if method[0].startswith("_"):
            continue
        print(method[0] + "()")


class PhysicalConstants:

    def __init__(self) -> None:
        self.Q_T_dict = {  # [無次元] partition function
            100: 7.36, 200: 20.726, 300: 37.608, 400: 57.649, 500: 80.579,
            600: 106.393, 700: 135.33, 800: 167.812, 900: 204.388, 1000: 245.683, 1200: 345.179}

    def h(self):
        mkhelp(self)


class InstrumentParameters:

    def __init__(
            self, N_read, I_dark, G_Amp, l_f, telescope_diameter) -> None:

        # 入力されたパラメーターの代入
        self.N_read = N_read
        self.I_dark = I_dark
        self.G_Amp = G_Amp
        self.l_f = l_f
        self.telescope_diameter = telescope_diameter

        # システムゲインの導出
        self.G_sys = self.__calc_G_sys()

        # 装置透過率の導出
        self.tau_t = 0.66
        self.tau_f = self.__calc_tau_f()
        self.tau_s = self.__calc_tau_s()
        self.tau_e = self.tau_t * self.tau_f * self.tau_s

        # ピクセル数関連の導出
        s_plate = 0.3  # <-ESPRITの値 プレートスケールは検出器までの光学系依存なのでTOPICSでは値が変わることに注意
        self.Omega = 2.35e-11 * s_plate
        w_slit = 0.7
        self.n_pix = w_slit / s_plate

        # その他の文字の定義
        self.A_t = np.pi * (self.telescope_diameter / 2) ** 2
        self.eta = 0.889

    def h(self):
        mkhelp(self)

    def __calc_G_sys(self):
        """システムゲインG_sysの計算

        Returns
        -------
        float
            システムゲイン
        """
        e = phys_consts.e
        G_Amp = self.G_Amp
        C_PD = 7.20e-14
        G_SF = 0.699
        ADU_ADC = 10 / 2**16

        G_sys = C_PD / (e * G_SF) * ADU_ADC / G_Amp
        return G_sys

    def __calc_tau_f(self):
        """InF3ファイバーの透過率計算

        Returns
        -------
        float
            ファイバー部分での透過率
        """
        l_f = self.l_f
        tau_f_unit = 0.98
        tau_f_coupling = 0.5

        tau_f = tau_f_coupling * tau_f_unit ** l_f
        return tau_f

    def __calc_tau_s(self):
        """ESPRITの装置透過率の計算

        Returns
        -------
        float
            ESPRIT全体での装置透過率
        """
        tau_s_lens = 0.66
        tau_s_mirror = 0.86

        # 実際は回折効率は波長依存性がかなりある（宇野2012D論p106）が、ひとまず固定値として計算
        tau_s_grating = 0.66

        tau_s = tau_s_lens * tau_s_mirror * tau_s_grating
        return tau_s


class ObservationParameters:

    def __init__(self, t_obs) -> None:

        # 入力パラメータの代入
        self.t_obs = t_obs
        self.tau_alpha = 0.9

        # 参照元に I_GBT + I_sky のみの合算値しかないので、実装では md上の $I_{GBT} + I_{sky}$ を I_GBT_sky として記述
        self.I_GBT_sky = 2.99e-6  # ESPRIT 3.4umでの値を仮置きした、今後もう少し複雑な機能を実装するかも

    def h(self):
        mkhelp(self)


class EmissionLineParameters:

    def __init__(
            self,
            lambda_um: float,
            g_ns: int,
            J_prime: int,
            A_if: float,
            E_prime: float) -> None:
        """__init__ 各輝線のパラメータ

        Parameters
        ----------
        lambda_um : float
            輝線の中心波長 [um]
        g_ns : int
            nuclear spin weight, 2 or 4 [無次元]
        J_prime : int
            回転準位 [無次元]
        A_if : float
            アインシュタインのA係数 [/s]
        E_prime : float
            energy of upper statement [/cm]
        """
        self.lambda_um = lambda_um
        self.omega_if = 1 / (self.lambda_um * 1e-6) * 1e-2  # 波数 [/cm]
        self.g_ns = g_ns
        self.J_prime = J_prime
        self.A_if = A_if
        self.E_prime = E_prime

    def h(self):
        mkhelp(self)
